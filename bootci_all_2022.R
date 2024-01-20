


cifmax=function(x, Year, cif){if(min(Year)>x){yy=0}else{
  yy=max(cif[Year<=x])}
  return(yy)}

# ACEI then ARB
#get KM estimate given weight
wtkm.b=function(dat, wt="fwt"){
  des=svydesign(ids=~1, weights=~eval(parse(text=wt)), data=dat) 
  s1<-svykm(Surv(death_days, AnyDeathEvent) ~ group, design=des)
  (trt0=s1[[1]])
  Year=trt0$time/365
  cif=1-trt0$surv
  CIF1=c(0,sapply(0.1*c(1:70),cifmax, Year=Year, cif=cif))
  trt1=s1[[2]]
  Year=trt1$time/365
  cif=1-trt1$surv
  CIF2=c(0,sapply(0.1*c(1:70),cifmax, Year=Year, cif=cif))
  return(c(CIF1, CIF2))}



#indices=sample(1:nrow(dboot), replace=TRUE)
#n=nrow(dboot)
#subgroup=c("Age75", "BLFemale", "Black","IntensiveTrt", "BLMCI")
#ncate=2
#dic<-c("ACEI","ARB")


hrcibootall = function(data,indices) { # bootstrap MI + ATE - all main tables 
dboot1=data[indices, ]
#dboot1=dboot[sample(1:nrow(dboot), replace=TRUE), ]

n=nrow(dboot1)
imp=mice(dboot1, m=10, predictorMatrix=pred1, meth=meth1, seed=3654, maxit=10)
z=dboot1$trt
## propensity score logistic regression model 
fit=with(imp, glm(eval(parse(text=ps.formula)), family=binomial))
## p scores... 
ps.list=lapply(fit$analyses, predict, type="response") ## ps for each imputed data;
#psmx=list.cbind(ps.list)
#ps.pool=apply(psmx, 1, mean)
                       
 
### Add ps and ATE weights to "imp"     
imp.long=complete(imp, action="long", include=TRUE)
imp.long$ps=c(rep(NA, n),unlist(ps.list))

imp.long=mutate(imp.long, ATE=case_when(trt==1~1/ps, trt==0~1/(1-ps)), 
                OW=case_when(trt==1~1-ps, trt==0~ps),
                minps=pmin(ps, 1-ps),
                mwt=case_when(trt==1~minps/ps, trt==0~minps/(1-ps)),
                BLMCI2=case_when(
                  is.na(BLMCI) & race %in% c("NH White", "Other") & BLeducation %in% c("High school diploma or GED","Less than high school diploma")~BLMoCA<=22,
                  is.na(BLMCI) & race %in% c("NH White", "Other") & BLeducation %in% c("Some college (No degree)")~BLMoCA<=24,
                  is.na(BLMCI) & race %in% c("NH White", "Other") & BLeducation %in% c("College graduate")~BLMoCA<=25,
                  is.na(BLMCI) & BLeducation %in% c("High school diploma or GED","Less than high school diploma")~BLMoCA<=23,
                  is.na(BLMCI) & BLeducation %in% c("Some college (No degree)")~BLMoCA<=24,
                  is.na(BLMCI) & BLeducation %in% c("College graduate")~BLMoCA<=24)%>% as.numeric,
                Black=as.numeric(race=="NH Black"),
                Age75=as.numeric(BLage>=75)) %>% mutate_at(
                  c(evtvars,"BLlivewithothers", "atrialFib", "BLdepress"), function(x){as.numeric(as.character(x))})

imp.long$BLMCI[is.na(imp.long$BLMCI)]=imp.long$BLMCI2[is.na(imp.long$BLMCI)]
imp.long$ATE=unlist(tapply(imp.long$ATE, imp.long$.imp, function(x){y=x; if(sum(!is.na(x))>0){cut=quantile(x, 0.99); y[x>cut]=cut};return(y)}))
imp1=as.mids(imp.long) ### imp1 has ATE added!


### now for each imputed dataset
ASMD1=NULL;
ASMD2=NULL;
ASMDb=NULL;
for (jj in 1:imp$m){
  db=complete(imp1,jj)
  wt1=db$ATE
  wt2=db$OW
  
  #also generate the baseline weight
  wt_bs<-rep(1,n)
  # the design matrix for balance check
  covM<-as.data.frame(model.matrix(formula(ps.formula),db))
  covM<-as.matrix(covM)
  if (ncol(covM)>1){
    #drop intercept
    if(unique(covM[,1])==1){
      covM<-covM[,-1,drop=F]
    }
  }
  ncov <- dim(covM)[2]+1   #number of covariate, plus baseline
  
  ######## Prepare subgroup matrix for balance check#############################################################################
  submatrix<-as.matrix(db[,subgroup,drop=FALSE])
  
  # Access the overall and subgroup asd, all level all_sub_level, by expanding subgroup by their levels
  all_sub_level <- list()
  for (k in 1:ncol(submatrix)){
    all_sub_level[[k]]<- table(submatrix[,k])
  }
  # number of subgroup levels by each subgrouping variable
  
  # subgroup levels, how namy level each subgroup
  nlevel <- c(unlist(lapply(all_sub_level,length)))
  nsubg <- c(unlist(all_sub_level))
  
  nv <- length(nsubg)   # number of subgroups
  subg.name_raw <- names(nsubg) #raw subgroup levels
  subcovnames <- colnames(submatrix) #raw subgroup names
  
  subg.name <- paste0(rep(subcovnames,nlevel),'=',subg.name_raw)#combine variable name to levels
  names(all_sub_level) <- subcovnames
  names(nsubg)<-subg.name #assign new subgroup names to summary table
  
  colnames(covM)
  overall_asd1 <- apply(covM, 2, abs_stand_diff, z, wt1) #ATE weighted
  overall_asd2 <- apply(covM, 2, abs_stand_diff, z, wt2) #ATE weighted
  overall_asd_bs <- apply(covM, 2, abs_stand_diff, z, wt_bs) #baseline
  
  #Calculate balance per subgroup across covariates
  groups_asd1 <- c() #weighted
  groups_asd2 <- c() #weighted
  groups_asd_bs <- c() #baseline
  names_col <- c()
  
  for(r in 1:ncol(submatrix)){
    level.name <- names(all_sub_level[[r]])
    for(g in 1:length(level.name)){
      find_g <- which(submatrix[,r]==level.name[g])
      g_asd1 <- apply(covM[find_g, ], 2, abs_stand_diff, z[find_g], wt1[find_g])
      groups_asd1 <- cbind(groups_asd1, g_asd1)
      g_asd2 <- apply(covM[find_g, ], 2, abs_stand_diff, z[find_g], wt2[find_g])
      groups_asd2 <- cbind(groups_asd2, g_asd2)
      g_asd_bs<-apply(covM[find_g, ], 2, abs_stand_diff, z[find_g], wt_bs[find_g])
      groups_asd_bs<- cbind(groups_asd_bs, g_asd_bs)
    }
  }
  
  colnames(groups_asd2) <-subg.name
  colnames(groups_asd1)  <-subg.name
  colnames(groups_asd_bs) <-subg.name
  
  ASD1 <- cbind(Overall=overall_asd1,groups_asd1) #weighted
  ASD2 <- cbind(Overall=overall_asd2,groups_asd2) #weighted
  ASD_bs <- cbind(Overall=overall_asd_bs,groups_asd_bs) #baseline
  
  nsubg <- c(length(z), nsubg)   # subgroup sample size including overall
  names(nsubg)[1]<-'Overall'
  
  ASMD1 = cbind(ASMD1, c(ASD1))
  ASMD2 =cbind(ASMD2, c(ASD2))
  ASMDb = cbind(ASMDb, c(ASD_bs))
}

asmdb=apply(ASMDb,1,mean)
asmd1=apply(ASMD1,1,mean)
asmd2=apply(ASMD2,1,mean)
asmd=c(Base=asmdb, WT.ate=asmd1, WT.OW=asmd2)

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##         CIFs 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cumcif.b=function(evtvar, dayvar, wt.type="ATE"){
  cif1=NULL
  cif0=NULL
  for (i in 1:10){
    dat = imp.long[imp.long$.imp==i, ]
    cif0= cbind(cif0, getcif(dat=dat, tx=0, evtvar=evtvar, dayvar=dayvar, wt.type=wt.type))
    cif1 =cbind(cif1, getcif(dat=dat,tx=1, evtvar=evtvar, dayva=dayvar, wt.type=wt.type ))
  }
  mean0= data.frame(Year=cif0[,1]/365, cif=apply(cif0[, 2*(1:10)],1,mean))
  mean1= data.frame(Year=cif1[,1]/365, cif=apply(cif1[, 2*(1:10)],1,mean))
  c(c(0,sapply(0.1*c(1:70),cifmax, Year=mean0$Year, cif=mean0$cif)),c(0,sapply(0.1*c(1:70),cifmax, Year=mean1$Year, cif=mean1$cif)))
}

## Death ###
dboot1$fwt=1#- unweighted
cd0=wtkm.b(dboot1, "fwt")
cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
cd1=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
cd2=apply(cd, 1, mean)
CIFd=c(death0=cd0, death1=cd1, death2=cd2)

### All primary/secondarey outcomes 
###  - pd_mci_amnestic- unweighted
rc=getcif(dat=dboot1, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=dboot1, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1, 0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
CIFe=c(cife0=ce0, cife1=ce1, cife2=ce2)

## PD MCI only
rc=getcif(dat=dboot1, tx=0, evtvar="mci_protocolcomp", dayvar="mci_protocol_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=dboot1, tx=1, evtvar="mci_protocolcomp", dayvar="mci_protocol_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
cse10=c(0,e1, 0, e2)
##Event wighted
cse11=cumcif.b("mci_protocolcomp", "mci_protocol_death_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
cse12=cumcif.b("mci_protocolcomp", "mci_protocol_death_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

SCIF1=c(pmci.un=cse10, pmci.ate=cse11, pmci.ow=cse12)

### Amnestic MCI alone
## PD MCI only
rc=getcif(dat=dboot1, tx=0, evtvar="mci_amnesticcomp", dayvar="mci_amnestic_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=dboot1, tx=1, evtvar="mci_amnesticcomp", dayvar="mci_amnestic_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
cse20=c(0,e1, 0, e2)
##Event wighted
cse21=cumcif.b("mci_amnesticcomp", "mci_amnestic_death_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
cse22=cumcif.b("mci_amnesticcomp", "mci_amnestic_death_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
SCIF2=c(amci.un=cse20, amci.ate=cse21, amci.ow=cse22)

### PD alone
rc=getcif(dat=dboot1, tx=0, evtvar="pdcomp", dayvar="pd_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=dboot1, tx=1, evtvar="pdcomp", dayvar="pd_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
cse30=c(0,e1, 0, e2)

cse31=cumcif.b("pdcomp", "pd_death_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
cse32=cumcif.b("pdcomp", "pd_death_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
SCIF3=c(pd.un=cse30, pd.ate=cse31, pd.ow=cse32)

## PD and MCI
rc=getcif(dat=dboot1, tx=0, evtvar="pd_mci_protocolcomp", dayvar="pd_mci_protocol_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=dboot1, tx=1, evtvar="pd_mci_protocolcomp", dayvar="pd_mci_protocol_death_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
cse40=c(0,e1,0, e2)
##Event wighted
cse41=cumcif.b("pd_mci_protocolcomp", "pd_mci_protocol_death_days") 
cse42=cumcif.b("pd_mci_protocolcomp", "pd_mci_protocol_death_days", "OW") 
SCIF4=c(pdmci.un=cse40, pdmci.ate=cse41, pdmci.ow=cse42)


## Female  Death ###
imps=filter(imp1, BLFemale==1);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)

sgd=c(female.d0=ssa, female.d1=ssb, female.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1, 0,e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(female.e0=ce0, female.e1=ce1, female.e2=ce2)

sg.female=c(sgd, sge)

## male  Death ###
imps=filter(imp1, BLFemale==0);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(male.d0=ssa, male.d1=ssb, male.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1, 0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(male.e0=ce0, male.e1=ce1, male.e2=ce2)

sg.male=c(sgd, sge)

##Age <75###
imps=filter(imp1, Age75==0);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(young.d0=ssa, young.d1=ssb, young.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1,0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(young.e0=ce0, young.e1=ce1, young.e2=ce2)

sg.young=c(sgd, sge)


##Age >=75###
imps=filter(imp1, Age75==1);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(old.d0=ssa, old.d1=ssb, old.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1, 0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(old.e0=ce0, old.e1=ce1, old.e2=ce2)
sg.old=c(sgd, sge)


###%%%%%%%%%%%%%%%%%%%%%%%%  Black
imps=filter(imp1, Black==1);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(black.d0=ssa, black.d1=ssb, black.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1,0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(black.e0=ce0, black.e1=ce1, black.e2=ce2)

sg.black=c(sgd, sge)


imps=filter(imp1, Black==0);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(nblack.d0=ssa, nblack.d1=ssb, nblack.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1,0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(nblack.e0=ce0, nblack.e1=ce1, nblack.e2=ce2)

sg.nonblack=c(sgd, sge)


##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Intensive
imps=filter(imp1, IntensiveTrt==1);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(Int.d0=ssa, Int.d1=ssb, Int.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1, 0,e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(Int.e0=ce0, Int.e1=ce1, Int.e2=ce2)

sg.Intensive=c(sgd, sge)


imps=filter(imp1, IntensiveTrt==0);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(Std.d0=ssa, Std.d1=ssb, Std.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1, 0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(Std.e0=ce0, Std.e1=ce1, Std.e2=ce2)

sg.Std=c(sgd, sge)

###BLMCI
imps=filter(imp1, BLMCI==1);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(MCI.d0=ssa, MCI.d1=ssb, MCI.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1,0, e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(MCI.e0=ce0, MCI.e1=ce1, MCI.e2=ce2)

sg.MCI=c(sgd, sge)


imps=filter(imp1, BLMCI==0);
cimp=complete(imps) %>% mutate(fwt=1)
imp.long=complete(imps, action="long");dim(imp.long)
ssa=wtkm.b(dat=cimp, wt="fwt")

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="ATE")
cd=cbind(cd, cifd)}
ssb=apply(cd, 1, mean)

cd=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifd=wtkm.b(dat=dat, wt="OW")
cd=cbind(cd, cifd)}
ssc=apply(cd, 1, mean)
sgd=c(nMCI.d0=ssa, nMCI.d1=ssb, nMCI.d2=ssc)

###  - pd_mci_amnestic- unweighted
rc=getcif(dat=cimp, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e1=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
rc=getcif(dat=cimp, tx=1, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)
e2=sapply(0.1*c(1:70),cifmax, Year=rc$Year, cif=rc$cif)
ce0=c(0,e1, 0,e2)
##Event wighted
ce1=cumcif.b("evtcomp", "comp_days") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()
ce2=cumcif.b("evtcomp", "comp_days", "OW") #%>% mutate(Year=round(Year, 1)%>% as.character())  %>% group_by(group, Year) %>% summarise(CIF=max(cif)) %>% as.data.frame()

sge=c(nMCI.e0=ce0, nMCI.e1=ce1, nMCI.e2=ce2)

sg.nMCI=c(sgd, sge)


##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
##                 Table 3 - weighted
##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## ATE
## pd_mci_amnestic
fit=with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE, weights=ATE))
est1=summary(pool(fit))$estimate 

# pd_mci_protocol
fit=with(imp1, coxph(Surv(pd_mci_protocol_days, pd_mci_protocol) ~ trt, robust=TRUE, weights=ATE))
est2=summary(pool(fit))$estimate

# pd
fit=with(imp1, coxph(Surv(pd_days, pd) ~ trt, robust=TRUE, weights=ATE))
est3=summary(pool(fit))$estimate
# mci_amnestic
fit=with(imp1, coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE, weights=ATE))
est4=summary(pool(fit))$estimate
# mci_protocol
fit=with(imp1, coxph(Surv(mci_protocol_days, mci_protocol) ~ trt, robust=TRUE, weights=ATE))
est5=summary(pool(fit))$estimate

# pd_mci_amnestic_death
fit=with(imp1, coxph(Surv(pd_mci_amnestic_death_days, pd_mci_amnestic_death) ~ trt, robust=TRUE, weights=ATE))
est6=summary(pool(fit))$estimate

# pd_mci_protocol_death
fit=with(imp1, coxph(Surv(pd_mci_protocol_death_days, pd_mci_protocol_death) ~ trt, robust=TRUE, weights=ATE))
est7=summary(pool(fit))$estimate

# pd_death
fit=with(imp1, coxph(Surv(pd_death_days, pd_death) ~ trt, robust=TRUE, weights=ATE))
est8=summary(pool(fit))$estimate

# mci_amnestic_death
fit=with(imp1, coxph(Surv(mci_amnestic_death_days, mci_amnestic_death) ~ trt, robust=TRUE, weights=ATE))
est9=summary(pool(fit))$estimate

# mci_protocol_death
fit=with(imp1, coxph(Surv(mci_protocol_death_days, mci_protocol_death) ~ trt, robust=TRUE, weights=ATE))
est10=summary(pool(fit))$estimate

# AnyDeathEvent
fit=with(imp1, coxph(Surv(death_days, AnyDeathEvent) ~ trt, robust=TRUE, weights=ATE))
est11=summary(pool(fit))$estimate
T3.ate=c(est1, est2, est3,est4, est5, est6, est7, est8, est9, est10, est11)

## OW
## pd_mci_amnestic
fit=with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE, weights=OW))
esto1=summary(pool(fit))$estimate 

# pd_mci_protocol
fit=with(imp1, coxph(Surv(pd_mci_protocol_days, pd_mci_protocol) ~ trt, robust=TRUE, weights=OW))
esto2=summary(pool(fit))$estimate

# pd
fit=with(imp1, coxph(Surv(pd_days, pd) ~ trt, robust=TRUE, weights=OW))
esto3=summary(pool(fit))$estimate
# mci_amnestic
fit=with(imp1, coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE, weights=OW))
esto4=summary(pool(fit))$estimate
# mci_protocol
fit=with(imp1, coxph(Surv(mci_protocol_days, mci_protocol) ~ trt, robust=TRUE, weights=OW))
esto5=summary(pool(fit))$estimate

# pd_mci_amnestic_death
fit=with(imp1, coxph(Surv(pd_mci_amnestic_death_days, pd_mci_amnestic_death) ~ trt, robust=TRUE, weights=OW))
esto6=summary(pool(fit))$estimate

# pd_mci_protocol_death
fit=with(imp1, coxph(Surv(pd_mci_protocol_death_days, pd_mci_protocol_death) ~ trt, robust=TRUE, weights=OW))
esto7=summary(pool(fit))$estimate

# pd_death
fit=with(imp1, coxph(Surv(pd_death_days, pd_death) ~ trt, robust=TRUE, weights=OW))
esto8=summary(pool(fit))$estimate

# mci_amnestic_death
fit=with(imp1, coxph(Surv(mci_amnestic_death_days, mci_amnestic_death) ~ trt, robust=TRUE, weights=OW))
esto9=summary(pool(fit))$estimate

# mci_protocol_death
fit=with(imp1, coxph(Surv(mci_protocol_death_days, mci_protocol_death) ~ trt, robust=TRUE, weights=OW))
esto10=summary(pool(fit))$estimate

# AnyDeathEvent
fit=with(imp1, coxph(Surv(death_days, AnyDeathEvent) ~ trt, robust=TRUE, weights=OW))
esto11=summary(pool(fit))$estimate
T3.ow=c(esto1, esto2, esto3,esto4, esto5, esto6, esto7, esto8, esto9, esto10, esto11)


# Unweighted
## pd_mci_amnestic
fit=with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE))
estu1=summary(pool(fit))$estimate

# pd_mci_protocol
fit=with(imp1, coxph(Surv(pd_mci_protocol_days, pd_mci_protocol) ~ trt, robust=TRUE))
estu2=summary(pool(fit))$estimate
# pd
fit=with(imp1, coxph(Surv(pd_days, pd) ~ trt, robust=TRUE))
estu3=summary(pool(fit))$estimate

# mci_amnestic
fit=with(imp1, coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE))
estu4=summary(pool(fit))$estimate

# mci_protocol
fit=with(imp1, coxph(Surv(mci_protocol_days, mci_protocol) ~ trt, robust=TRUE))
estu5=summary(pool(fit))$estimate
# pd_mci_amnestic_death
fit=with(imp1, coxph(Surv(pd_mci_amnestic_death_days, pd_mci_amnestic_death) ~ trt, robust=TRUE))
estu6=summary(pool(fit))$estimate
# pd_mci_protocol_death
fit=with(imp1, coxph(Surv(pd_mci_protocol_death_days, pd_mci_protocol_death) ~ trt, robust=TRUE))
estu7=summary(pool(fit))$estimate
# pd_death
fit=with(imp1, coxph(Surv(pd_death_days, pd_death) ~ trt, robust=TRUE))
estu8=summary(pool(fit))$estimate
# mci_amnestic_death
fit=with(imp1, coxph(Surv(mci_amnestic_death_days, mci_amnestic_death) ~ trt, robust=TRUE))
estu9=summary(pool(fit))$estimate
# mci_protocol_death
fit=with(imp1, coxph(Surv(mci_protocol_death_days, mci_protocol_death) ~ trt, robust=TRUE))
estu10=summary(pool(fit))$estimate
# AnyDeathEvent
fit=with(imp1, coxph(Surv(death_days, AnyDeathEvent) ~ trt, robust=TRUE))
estu11=summary(pool(fit))$estimate
T3.un=c(estu1, estu2, estu3,estu4, estu5, estu6, estu7, estu8, estu9, estu10, estu11)


#### Subgroup - primary outcomes. 
so2=function(imps,wt.type=NULL){ if(is.null(wt.type)) {fit=with(imps, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt))}else{
                                       fit=with(imps, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE, weights=eval(parse(text=wt.type))))}
  return(est=summary(pool(fit))$estimate)}

sgest.ate=c(so2(filter(imp1,BLage<75),wt.type="ATE"),
so2(filter(imp1,BLage>=75),wt.type="ATE"),

so2(filter(imp1,BLFemale==0),wt.type="ATE"),
so2(filter(imp1,BLFemale==1),wt.type="ATE"),

so2(filter(imp1,Black==1),wt.type="ATE"),
so2(filter(imp1,Black==0),wt.type="ATE"),

so2(filter(imp1,IntensiveTrt==0),wt.type="ATE"),
so2(filter(imp1,IntensiveTrt==1),wt.type="ATE"),

so2(filter(imp1,BLMCI==0),wt.type="ATE"),
so2(filter(imp1,BLMCI==1),wt.type="ATE"))


sgest.ow=c(so2(filter(imp1,BLage<75),wt.type="OW"),
  so2(filter(imp1,BLage>=75),wt.type="OW"),
  
  so2(filter(imp1,BLFemale==0),wt.type="OW"),
  so2(filter(imp1,BLFemale==1),wt.type="OW"),
  
  so2(filter(imp1,Black==1),wt.type="OW"),
  so2(filter(imp1,Black==0),wt.type="OW"),

  
  so2(filter(imp1,IntensiveTrt==0),wt.type="OW"),
  so2(filter(imp1,IntensiveTrt==1),wt.type="OW"),
  
  so2(filter(imp1,BLMCI==0),wt.type="OW"),
  so2(filter(imp1,BLMCI==1),wt.type="OW"))


sgest.un=c(so2(filter(imp1,BLage<75)),
           so2(filter(imp1,BLage>=75)),
           
           so2(filter(imp1,BLFemale==0)),
           so2(filter(imp1,BLFemale==1)),
           
           so2(filter(imp1,Black==1)),
           so2(filter(imp1,Black==0)),

           so2(filter(imp1,IntensiveTrt==0)),
           so2(filter(imp1,IntensiveTrt==1)),
           
           so2(filter(imp1,BLMCI==0)),
           so2(filter(imp1,BLMCI==1)))

           

# Unweighted
# minimally adjusted
hrb=NULL
fit=with(imp1,coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic)~ trt+BLage+BLFemale+BLrace))
summary(pool(fit))
est=summary(pool(fit))$estimate[1]
hrb=c(hrb, exp(est))

## Multivariate adjusted
fm2=paste("Surv(pd_mci_amnestic_days, pd_mci_amnestic)~trt+", paste(predvars, collapse="+"))
fit <- with(imp1,coxph(eval(parse(text=fm2))))
est=summary(pool(fit))$estimate[1]
hrb=c(hrb, exp(est))

## propensity score as covariate
fit <- with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic)~trt+ps ))
est=summary(pool(fit))$estimate[1]
hrb=c(hrb, exp(est))

### propensity score stratification
#https://rdrr.io/cran/MatchIt/man/method_subclass.html
imp.data <- complete(imp, "long")
fit.list <- setNames(vector("list", length(unique(imp.data$.imp))),  unique(imp.data$.imp))
for (i in unique(imp.data$.imp)) {
  m.out <- matchit(eval(parse(text=ps.formula)), method = "subclass", subclass=10, estimand = "ATE", data = imp.data[imp.data$.imp == i,])
  dat_m= match.data(m.out, data = imp.data[imp.data$.imp == i,])
fit.list[[i]] <- coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE, weights=weights, data=dat_m)
}
fit.list.mira <- as.mira(fit.list) #combines into mira object for pool()
hr=exp(summary(pool(fit.list.mira))$estimate[1]);names(hr)="Stratification"
hrb=c(hrb, hr)

## Propensity score matching
fit.list <- setNames(vector("list", length(unique(imp.data$.imp))),  unique(imp.data$.imp))
for (i in unique(imp.data$.imp)) {
  m.out <- matchit(eval(parse(text=ps.formula)), method = "nearest", caliper=0.2, caliper.sd=TRUE,data = imp.data[imp.data$.imp == i,])
  dat_m= match.data(m.out, data = imp.data[imp.data$.imp == i,])
fit.list[[i]] <- coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE, weights=weights, data=dat_m)
}
fit.list.mira <- as.mira(fit.list) #combines into mira object for pool()
hr=exp(summary(pool(fit.list.mira))$estimate[1]); names(hr)="Matching"
hrb=c(hrb, hr)

## Matching weight adjusted
fit=with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE, weights=mwt))
est=summary(pool(fit))$estimate
hrb=c(hrb, exp(est))
hrs.s2=hrb


### AE?
aehr=function(evt, timevar, wt.type=NULL){
fm=paste("Surv(", timevar, ",", evt, ")~trt", sep="")
if(is.null(wt.type)){fit=with(imp1, coxph(eval(parse(text=fm))))}else{
fit=with(imp1, coxph(eval(parse(text=fm)), robust=TRUE, weights=eval(parse(text=wt.type))))}
est=summary(pool(fit))$estimate
hr=exp(est);
names(hr)=evt
return(hr)}


## unweighted
ae=aehr(evt="AnySAEEvent", timevar="t_AnySAEEvent");
ae1=aehr("hypotenSAEEvent" , timevar="t_hypotenSAEEvent") ## warnings
ae2=aehr(evt="SyncopeSAEEvent" , timevar="t_SyncopeSAEEvent") ## warnings
ae3=aehr(evt="BradySAEEvent" , timevar="t_BradySAEEvent") # warnings
ae4=aehr(evt="electroSAEEvent" , timevar="t_electroSAEEvent")
ae6=aehr(evt="AKISAEEvent" , timevar="t_AKISAEEvent" )
ae7=aehr(evt="hypotenallEvent" , timevar="t_hypotenallEvent") 
ae8=aehr(evt="SyncopeallEvent" , timevar="t_SyncopeallEvent") 
ae9=aehr(evt="BradyallEvent" , timevar="t_BradyallEvent")
ae10=aehr(evt="electroallEvent" , timevar="t_electroallEvent")
ae12=aehr(evt="AKIallEvent" , timevar="t_AKIallEvent")
ae13=aehr(evt="LowNaEvent" , timevar="t_LowNaEvent")
ae15=aehr(evt="LowKEvent" , timevar="t_LowKEvent")
ae16=aehr(evt="HighKEvent" , timevar="t_HighKEvent")
ae17=aehr(evt="OHmeasuredEvent" , timevar="t_LowKEvent" )
ae18=aehr(evt="OHbothEvent" , timevar="t_OHbothEvent")
aes.un=c(ae,  ae1, ae2, ae3, ae4, ae6,  ae7,ae8, ae9, ae10,  ae12, ae13, ae15, ae16, ae17, ae18)
aes.un

## ATE
ae=aehr(evt="AnySAEEvent", timevar="t_AnySAEEvent", wt.type = "ATE");
ae1=aehr("hypotenSAEEvent" , timevar="t_hypotenSAEEvent",wt.type = "ATE") ## warnings
ae2=aehr(evt="SyncopeSAEEvent" , timevar="t_SyncopeSAEEvent" , wt.type = "ATE" ) ## warnings
ae3=aehr(evt="BradySAEEvent" , timevar="t_BradySAEEvent", wt.type = "ATE"  ) # warnings
ae4=aehr(evt="electroSAEEvent" , timevar="t_electroSAEEvent", wt.type = "ATE")
ae6=aehr(evt="AKISAEEvent" , timevar="t_AKISAEEvent" , wt.type = "ATE")
ae7=aehr(evt="hypotenallEvent" , timevar="t_hypotenallEvent", wt.type = "ATE") 
ae8=aehr(evt="SyncopeallEvent" , timevar="t_SyncopeallEvent", wt.type = "ATE") 
ae9=aehr(evt="BradyallEvent" , timevar="t_BradyallEvent", wt.type = "ATE")
ae10=aehr(evt="electroallEvent" , timevar="t_electroallEvent", wt.type = "ATE")
ae12=aehr(evt="AKIallEvent" , timevar="t_AKIallEvent", wt.type = "ATE")
ae13=aehr(evt="LowNaEvent" , timevar="t_LowNaEvent", wt.type = "ATE")
ae15=aehr(evt="LowKEvent" , timevar="t_LowKEvent", wt.type = "ATE")
ae16=aehr(evt="HighKEvent" , timevar="t_HighKEvent", wt.type = "ATE")
ae17=aehr(evt="OHmeasuredEvent" , timevar="t_LowKEvent" , wt.type = "ATE")
ae18=aehr(evt="OHbothEvent" , timevar="t_OHbothEvent", wt.type = "ATE")
aes.ate=c(ae,  ae1, ae2, ae3, ae4, ae6,  ae7,ae8, ae9, ae10,  ae12, ae13, ae15, ae16, ae17, ae18)
aes.ate


## OW
ae=aehr(evt="AnySAEEvent", timevar="t_AnySAEEvent", wt.type = "OW");
ae1=aehr("hypotenSAEEvent" , timevar="t_hypotenSAEEvent",wt.type = "OW") ## warnings
ae2=aehr(evt="SyncopeSAEEvent" , timevar="t_SyncopeSAEEvent" , wt.type = "OW" ) ## warnings
ae3=aehr(evt="BradySAEEvent" , timevar="t_BradySAEEvent", wt.type = "OW"  ) # warnings
ae4=aehr(evt="electroSAEEvent" , timevar="t_electroSAEEvent", wt.type = "OW")
ae6=aehr(evt="AKISAEEvent" , timevar="t_AKISAEEvent" , wt.type = "OW")
ae7=aehr(evt="hypotenallEvent" , timevar="t_hypotenallEvent", wt.type = "OW") 
ae8=aehr(evt="SyncopeallEvent" , timevar="t_SyncopeallEvent", wt.type = "OW") 
ae9=aehr(evt="BradyallEvent" , timevar="t_BradyallEvent", wt.type = "OW")
ae10=aehr(evt="electroallEvent" , timevar="t_electroallEvent", wt.type = "OW")
ae12=aehr(evt="AKIallEvent" , timevar="t_AKIallEvent", wt.type = "OW")
ae13=aehr(evt="LowNaEvent" , timevar="t_LowNaEvent", wt.type = "OW")
ae15=aehr(evt="LowKEvent" , timevar="t_LowKEvent", wt.type = "OW")
ae16=aehr(evt="HighKEvent" , timevar="t_HighKEvent", wt.type = "OW")
ae17=aehr(evt="OHmeasuredEvent" , timevar="t_LowKEvent" , wt.type = "OW")
ae18=aehr(evt="OHbothEvent" , timevar="t_OHbothEvent", wt.type = "OW")
aes.ow=c(ae,  ae1, ae2, ae3, ae4, ae6,  ae7,ae8, ae9, ae10,  ae12, ae13, ae15, ae16, ae17, ae18)
aes.ow

### Negative control
so2=function(imps,wt.type=NULL){ if(is.null(wt.type)) {fit=with(imps, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt))}else{
  fit=with(imps, coxph(Surv(t_CompositeEvent, CompositeEvent) ~ trt, robust=TRUE, weights=eval(parse(text=wt.type))))}
  return(est=summary(pool(fit))$estimate)}


nsgest.ate=c(so2(imp1,wt.type="ATE"),
  so2(filter(imp1,BLage<75),wt.type="ATE"),
            so2(filter(imp1,BLage>=75),wt.type="ATE"),
            
            so2(filter(imp1,BLFemale==0),wt.type="ATE"),
            so2(filter(imp1,BLFemale==1),wt.type="ATE"),
            
            so2(filter(imp1,Black==1),wt.type="ATE"),
            so2(filter(imp1,Black==0),wt.type="ATE"),

            
            so2(filter(imp1,IntensiveTrt==0),wt.type="ATE"),
            so2(filter(imp1,IntensiveTrt==1),wt.type="ATE"),
  
  so2(filter(imp1,BLMCI==0),wt.type="ATE"),
  so2(filter(imp1,BLMCI==1),wt.type="ATE"))

nsgest.ow=c(so2(imp1,wt.type="OW"),
  so2(filter(imp1,BLage<75),wt.type="OW"),
           so2(filter(imp1,BLage>=75),wt.type="OW"),
           
           so2(filter(imp1,BLFemale==0),wt.type="OW"),
           so2(filter(imp1,BLFemale==1),wt.type="OW"),
           
           so2(filter(imp1,Black==1),wt.type="OW"),
           so2(filter(imp1,Black==0),wt.type="OW"),
           
           so2(filter(imp1,IntensiveTrt==0),wt.type="OW"),
           so2(filter(imp1,IntensiveTrt==1),wt.type="OW"),
  
  so2(filter(imp1,BLMCI==0),wt.type="OW"),
  so2(filter(imp1,BLMCI==1),wt.type="OW") )

nsgest.un=c(so2(imp1),
  so2(filter(imp1,BLage<75)),
           so2(filter(imp1,BLage>=75)),
           
           so2(filter(imp1,BLFemale==0)),
           so2(filter(imp1,BLFemale==1)),
           
           so2(filter(imp1,Black==1)),
           so2(filter(imp1,Black==0)),
           
           
           so2(filter(imp1,IntensiveTrt==0)),
           so2(filter(imp1,IntensiveTrt==1)),
  
  so2(filter(imp1,BLMCI==0)),
  so2(filter(imp1,BLMCI==1)) )



out=c(asmd, CIFd, CIFe, SCIF1, SCIF2, SCIF3, SCIF4, sg.young, sg.old,sg.male, sg.female, sg.black, sg.nonblack, sg.Intensive, sg.Std, sg.MCI, sg.nMCI, 
  t3ate=T3.ate, t3ow=T3.ow, t3un=T3.un,sgest.un=sgest.un, sgest.ate=sgest.ate, sgest.ow=sgest.ow, hrs.adj=hrs.s2, AEs=c(aes.un,aes.ate,aes.ow), 
  nsg.un=nsgest.un, nsg.ate=nsgest.ate, nsg.ow=nsgest.ow )

return(out)}









