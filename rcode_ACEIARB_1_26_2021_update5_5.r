setwd("C:/Users/u0942431/Box/Zhang, Chong - SDBC/macro");
source("functions_2019.r");source("functions2.r");
require(survey);require(MatchIt); require(mice);library(boot);require(tidyverse);require(plyr);library(viridis); require(car);require(rlist);require(twang)
require(broom.mixed);require(ggpubr);require(survminer);require(glmnet)

setwd("C:/Users/u0942431/Box/Zhang, Chong - SDBC/SPRINTMIND234cz/MIND new users ACEI or ARB/Analysis")
source("support_functions_2022.r")

d=read.csv("C:/Users/u0942431/Box/SPRINT 2018/MIND new users ACEI or ARB/Data/cogdata12.csv", na.strings=c(""," ","Null","NULL"))
#d=read.csv("C:/Users/u0942431/Box/SPRINT 2018/MIND new users ACEI or ARB/Data/cogdata24.csv", na.strings=c(""," ","Null","NULL"))
nuni(d$maskid)
#d$maskid
tablem(d$BLMCI)
glimpse(d);dim(d) # 2040  119
nuni(d$maskid) #2040  /  2291
table(d$ACEIvsARB)#727 1313

## Add histogram
#function(data, nbins=50){ggplot(data,aes(ps))+geom_histogram(alpha=0.5,bins =nbins)+ylab(" ")+facet_wrap(~group,ncol=1)+xlab("Propensity score")}
dp=d %>% mutate(grp=factor(ACEIvsARB))
levels(dp$grp)=c('ARB','ACEI')
ggplot(dp, aes(indexdate))+geom_histogram(alpha=0.5, bins=30)+ylab('')+facet_wrap(~grp,ncol=1)+xlab("Days from Initiation")+scale_x_continuous(breaks=c(0,60, 120,180, 240,300,360))
ggsave('Days_from.pdf' , width=7, height=6)


###%%%%%%%%%%%%%%%%%%%%      %%%    Cleaning of data, Recoding variables for later use        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d$BLeducation=factor(d$BLeducation)
table(d$BLeducation)
levels(d$BLeducation)=c("Less than high school diploma", "High school diploma or GED", "Some college (No degree)", "College graduate")
d$BLsmoking=factor(d$BLsmoking)
table(d$BLsmoking)
levels(d$BLsmoking)=c("Current", "Former", "Never")
tablem(d$BLcaresource)
d$BLcaresource=factor(d$BLcaresource)
levels(d$BLcaresource)=c("Doctors office/outpatient clinic", "Community healthcare facility/other", "No usual source of care")
d$BLrace=factor(d$BLrace)
table(d$BLrace)
levels(d$BLrace)= c("NH White", "NH Black", "Hispanic", "Other")
anyins=d[,c("privOther","medicare","medicaid","va")] %>% rowSums()
table(anyins, d$uninsured)

## %%%%%%%%%%%%%%%%%%%%%%%%%       &&&    Descriptive Summary -- just raw data          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
covars=c("BLage","BLFemale", "BLrace","BLlivewithothers","BLsmoking","BLeducation","uninsured","medicare","medicaid","va","privOther","BLcaresource",
         "BLClinicalCVD" ,"atrialFib",  "BLdepress",  "BLMoCA","BLLogicalMemoryII","BLDigitSymbolCoding", "BLSBP","BLDBP","BLseatHeartRate","BLpotassium","BLCreatr","BLUMALCR", 
         "BLtotChol","BLHDL","BLTriglycerides","BLBMI", "BLGlucose","BLeGFR","BLAspirinYN","BLStatinYN","BLNsaidYN", "BLNonAntiHyperMedNum","BLAntiHyperMedNum",  "BLCCB", 
         "BLThiazideDiuretic" ,"BLLoopDiuretic", "BLBetaBlocker","BLAlphaBlocker","BLOtherAntiHyperMed","IntensiveTrt")

predvars=setdiff(covars, c("BLeGFR", "uninsured"))  ## variables involvedl in ps model
length(predvars)

varlabs=c("Age","Female","Race","live with others","Smoking status","Education","Insurance status- Uninsured","Having Medicare","Having Medicaid","Having VA insurance","Having Private insurance",
          "Source of care", "Medical history-Clinical CVD","-AFib","-Depression",   "Montreal Cognitive","Logical Memory Form II","Digit Symbol Coding","Systolic BP","Diastolic BP",
          "Resting heart rate","Serum potassium","Serum creatinine","Albumin to creatinine ratio","Total cholesterol","HDL cholesterol","Triglycerides", "BMI", "Serum glucose", "eGFR",  
          "Baseline medication - Aspirin","-Statin","-NSAID","-Number of non antihypertensive medications", "-Number of antihypertensive medications", "CCB","Thiazide diuretic","Loop diuretic",  
          "Beta blocker","Alpha blocker","Other antihypertensive medication","Intensive arm")
glimpse(d[, covars])
tp=guess(d[,covars], 6)
tp[tp==1]=0

dig=rep(0, length(tp))
tp[covars %in% c("BLNonAntiHyperMedNum","BLAntiHyperMedNum","BLpotassium","BLCreatr","BLUMALCR", "BLtotChol","BLHDL","BLTriglycerides","BLBMI", "BLGlucose", "BLMoCA","BLLogicalMemoryII","BLDigitSymbolCoding")]=1
cbind(covars, varlabs, tp)
classvars=covars[tp>1]
classvars

t1=my2grtab(xdat = d[,covars], gvar=d$ACEIvsARB, gnames=c("ACEI", "ARB"), typev=tp, dnames=varlabs, combinegroups = T,footMissing = FALSE, allSums = FALSE, onelinebinary = TRUE, showbinarylevel = FALSE, output="Table0.doc")# 
Table0=prfmt(t1[[1]])

dim(d)
#AnyDeathEvent death_days
#Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt,
#pd_mci_amnestic_death_days, pd_mci_amnestic_death
d$status=1-d$pd_mci_amnestic
d$year=d$pd_mci_amnestic_days/365
survfit(Surv(year, status)~1, data=d)
table(d$ACEIvsARB)


x=d$pd_mci_amnestic_days/365
miqr(x)


tp1=tp[covars %in% predvars]
varlabs1=varlabs[covars %in% predvars]

## Our main events and time to events
daysvars=findchar("_days", names(d))
evtvars=sub("_days","", daysvars)
evtvars[daysvars=="death_days"]="AnyDeathEvent"  
cbind(evtvars, daysvars)

dim(d)
d[, daysvars] %>% summary

##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###       Defining outcome for the competing risk CIF 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tablem(d$AnyDeathEvent)
tablem(d$pd_mci_amnestic)
tablem(d$pd_mci_amnestic_death)
d$evtcomp=d$pd_mci_amnestic_death;
d$evtcomp[!d$pd_mci_amnestic %in% 1 & d$evtcomp==1]=2
tablem(d$evtcomp)
tablem(d$evtcomp, d$pd_mci_amnestic)
d$comp_days=d$pd_mci_amnestic_death_days

d$pdcomp=d$pd_death;
d$pdcomp[!d$pd %in% 1 & d$pdcomp==1]=2
tablem(d$pdcomp)

d$mci_protocolcomp=d$mci_protocol_death;
d$mci_protocolcomp[!d$mci_protocol %in% 1 & d$mci_protocolcomp==1]=2
table(d$mci_protocolcomp, d$mci_protocol)
d$pd_mci_protocolcomp=d$pd_mci_protocol_death;
d$pd_mci_protocolcomp[!d$pd_mci_protocol %in% 1 & d$pd_mci_protocolcomp==1]=2
table(d$pd_mci_protocolcomp, d$pd_mci_protocol)
table(d$pd_mci_protocol, d$AnyDeathEvent)
d$mci_amnesticcomp=d$mci_amnestic_death;
d$mci_amnesticcomp[!d$mci_amnestic %in% 1 & d$mci_amnesticcomp==1]=2
table(d$mci_amnesticcomp, d$mci_amnestic)

###%%%%%%%%%%%%%%       &&&&&&&&&&&&&             Missing Table 2                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outmis=names(findmiss(d[, evtvars]) )#%>% as.data.fr
outlabs=c("Probable Dimentia", "Protocol defined MCI","Amnestic MCI", "Probable Dimentia or Amnestic MCI", "Probable Dimentia or Amnestic MCI or Death")
cbind(outlabs, outmis)
#printtab(outmis, "Outcomemiss.doc")
findmiss(d %>% select(ends_with("comp")))

covarsmis=names(findmiss(d[, predvars]))
lab=varlabs[covars %in% covarsmis]
laball=c(lab, outlabs)

dmis=is.na(d[, c(covarsmis, outmis)]) %>% as.data.frame
head(dmis)
t1m=my2grtab(xdat = dmis,gvar=d$ACEIvsARB, gnames=c("ACEI","ARB"), typev=rep(2, ncol(dmis)), dnames=laball,onelinebinary=TRUE, showbinarylevel=FALSE, combinegroups = TRUE, output="table2.doc")
(Table2=t1m[[1]] %>% select(-Levels))

## follow up time KM estimation
d$fucensor=1-d$pd_mci_amnestic
survfit(Surv(pd_mci_amnestic_days/365, fucensor)~1, data=d)
#      n  events  median 0.95LCL 0.95UCL 
# 2039.00 1661.00    4.89    4.59    5.05 

dm1=d %>% mutate(trt=1-ACEIvsARB) %>% select(-BLLeftVentricalHypertrophyany3, -BLcollege_grad, -BLMI)  
dm1$group=factor(dm1$trt)
levels(dm1$group)=c("ACEI", "ARB")
table(dm1$trt)
table(dm1$ACEIvsARB)
table(dm1$group);dim(dm1)
tablem(dm1$BLMCI)

dm1$race=dm1$BLrace
tablem(dm1$BLrace)
tablem(dm1$BLeducation)
levels(dm1$BLrace )=c("NH White" ,"NH Black" ,"Other" , "Other") 
tablem(dm1$race)

setdiff(names(dm1),predvars)
setdiff(predvars,names(dm1))
nuni(dm1$maskid)

## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###       MI, Weighting, and all analysis of effect - point estimate using original data (no resampling..)                 ##############
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dboot=dm1  %>% mutate_at(c("BLlivewithothers", "atrialFib", "BLdepress"), as.factor) ## do we impute outcomes? - no
findmiss(dboot)
nuni(dboot$maskid)
mvars=names(findmiss(dboot));length(mvars) # 29
glimpse(dboot[,mvars])

outcomemiss=intersect(mvars, c(evtvars, daysvars)) 
outcomemiss

## MI now
ini=mice(dboot, maxit=0)
(meth1=ini$meth)
findchar("comp", names(d))
meth1[mvars] 
#intersect(mvars, findchar("comp", names(d)))
meth1[c(outcomemiss, "fucensor", "evtcomp", "comp_days", "BLMCI") ]="" # do not impute these

puse=c("trt", predvars) ## only use these
pdrop=setdiff(names(dboot), puse)
pred1=quickpred(dboot, mincor=0.1, exclude=pdrop)

dim(pred1)
pred1[mvars,] %>% rowSums ## total predictors per missing variable
colSums(pred1)[colSums(pred1)>0] ## all variables ever used as predictor

mvars1=setdiff(mvars, outcomemiss)
pred1[mvars1, ] %>% rowSums %>% sort
meth1[meth1!=""] 


imp=mice(dboot, m=10, predictorMatrix=pred1, meth=meth1, seed=3654, maxit=10)
findmiss(complete(imp,2))
#save(imp, file="Imputed24.RData")
#names(dboot)
#findchar("comp", names(dboot))
#load("Imputed24.RData")

subgroup=c("Age75", "BLFemale", "Black","IntensiveTrt", "BLMCI")

require(boot);require(dplyr);require(doParallel); require(doRNG);

ncore=detectCores()
cl <- makeCluster(ncore-1) # don't use all your cores 
registerDoParallel(cl)

## bootstrap - takes days 
#source("bootci_all_2022short.r")

#t1=Sys.time()
#registerDoRNG(12345)
#bt=boot(dboot, hrcibootCIFselect, R=12,parallel = 'snow',ncpus=ncore-1  ) # 1 hour per 50 runs
#t2=Sys.time()
#t2-t1
#save(bt, file="boot2022.RData")


###  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###          This section does the check for important interactions ###  
###  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#setwd("C:/Users/u0942431/Box/Zhang, Chong - SDBC/SPRINTMIND234cz/MIND new users ACEI or ARB/Analysis/SGAcode/Design")
#source("SumStat_sub.R")
#source("print_Sumstat_sub.R")
#source("Plot_SumStat_sub.R")
#source('PSmethod_sub.R')
#setwd("C:/Users/u0942431/Box/Zhang, Chong - SDBC/SPRINTMIND234cz/MIND new users ACEI or ARB/Analysis")

## We modify the predictor list by relacing BLrace with 'Black' because we don't want to have them both
predvars1=c("BLage","BLFemale", "Black","BLlivewithothers","BLsmoking","BLeducation","medicare","medicaid","va","privOther","BLcaresource",
           "BLClinicalCVD" ,"atrialFib",  "BLdepress",  "BLMoCA","BLLogicalMemoryII","BLDigitSymbolCoding", "BLSBP","BLDBP","BLseatHeartRate","BLpotassium","BLCreatr","BLUMALCR", 
           "BLtotChol","BLHDL","BLTriglycerides","BLBMI", "BLGlucose","BLAspirinYN","BLStatinYN","BLNsaidYN", "BLNonAntiHyperMedNum","BLAntiHyperMedNum",  "BLCCB", 
           "BLThiazideDiuretic" ,"BLLoopDiuretic", "BLBetaBlocker","BLAlphaBlocker","BLOtherAntiHyperMed","IntensiveTrt")
length(predvars1) #40 the reason is we use 'Black' instead of BLrace

ps.form_m <- paste("trt~",paste(predvars1,collapse = "+")) # main effect model 
ps.form_f = ps.form_m 
for(i in 1:length(subgroup)){
  predvars2=setdiff(predvars1, subgroup[i])
  for (j in 1:length(predvars2)){
    ps.form_f=paste(ps.form_f, paste(subgroup[i],"*",predvars2[j]), sep=" + ")
  }}
fm=as.formula(ps.form_f)  ## use this formula to check interactions

## Cycle through 10 imputed data and save the variables selected.
Kept=vector("list", 10)
for (m in 1:10){
  dd=complete(imp, m) # one of the 10 imputed dataset - for now just use this.
  dd$Age75=as.numeric(dd$BLage>=75)
  dd$Black=(dd$BLrace=="NH Black")
  dd1=filter(dd, is.na(BLMCI)) %>%mutate (
  BLMCI=case_when(BLrace %in% c("NH White", "Other") & BLeducation %in% c("High school diploma or GED","Less than high school diploma")~BLMoCA<=22,
                    BLrace%in% c("NH White", "Other") & BLeducation %in% c("Some college (No degree)")~BLMoCA<=24,
                    BLrace%in% c("NH White", "Other") & BLeducation %in% c("College graduate")~BLMoCA<=25,
                    BLeducation %in% c("High school diploma or GED","Less than high school diploma")~BLMoCA<=23,
                    BLeducation %in% c("Some college (No degree)")~BLMoCA<=24,
                    BLeducation %in% c("College graduate")~BLMoCA<=24)%>% as.numeric)
  dd2=filter(dd, !is.na(BLMCI))
  dd=rbind(dd1, dd2)
  summary(dd$BLUMALCR)
  fullmatrix <- model.matrix(fm, dd) ## update fullmatrix-- use same fm but different data
  
  # create an idx for penalty.factor 
  idx <- rep(0, dim(fullmatrix)[2])
  idx[grep('\\:',colnames(fullmatrix))] <- 1
  idx[colnames(fullmatrix) %in% "BLMCI"] <- 1
  cbind(colnames(fullmatrix),idx) 
  set.seed(100+m)
  fitLASSO <- cv.glmnet(y=dd$trt, x=fullmatrix, penalty.factor=idx, family="binomial", maxit=50000)
  nonzero_coef <- rownames(coef(fitLASSO, s='lambda.min'))[which(coef(fitLASSO, s='lambda.min')!=0)][-1]
  nonzero_coef
  Kept[[m]]=nonzero_coef}

unlist(lapply(Kept, length))
Kept 
# For most imputed datasets, none of the interactions were selected, so from now on we just use the main effects model for propensity score.

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###     Checking is done ## we now move on and use main effect model.
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
##            propensity score logistic regression model 
##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

ps.formula="trt ~ BLage + BLFemale + BLrace + BLlivewithothers + BLsmoking + 
  BLeducation + medicare + medicaid + va + privOther + BLcaresource + 
  BLClinicalCVD + atrialFib + BLdepress + BLMoCA + BLLogicalMemoryII + 
  BLDigitSymbolCoding + BLSBP + BLDBP + BLseatHeartRate + BLpotassium + 
  BLCreatr + BLUMALCR + BLtotChol + BLHDL + BLTriglycerides + 
  BLBMI + BLGlucose + BLAspirinYN + BLStatinYN + BLNsaidYN + 
  BLNonAntiHyperMedNum + BLAntiHyperMedNum + BLCCB + BLThiazideDiuretic + 
  BLLoopDiuretic + BLBetaBlocker + BLAlphaBlocker + BLOtherAntiHyperMed + 
  IntensiveTrt"  #### If we will eventually use the main effect model, this is it.

fit=with(imp, glm(eval(parse(text=ps.formula)), family=binomial))
ps.list=lapply(fit$analyses, predict, type="response") ## ps for each imputed data;
psmx=list.cbind(ps.list)
ps.pool=apply(psmx, 1, mean)
summary(ps.pool)
dm1$ps=ps.pool
fig1=hist.ps(dm1)


n=nrow(dboot)
ncate=2
z=dboot$trt

### Add ps and ATE weights to "imp"     
imp.long=complete(imp, action="long", include=TRUE)
imp.long$ps=c(rep(NA, n),unlist(ps.list))

#findchar("comp",names(imp.long))

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
                Age75=as.numeric(BLage>=75)) %>% mutate_at( c(evtvars,"BLlivewithothers", "atrialFib", "BLdepress"), function(x){as.numeric(as.character(x))})
				

imp.long$BLMCI[is.na(imp.long$BLMCI)]=imp.long$BLMCI2[is.na(imp.long$BLMCI)]
imp.long$BLMCI %>% summary

max(imp.long$ATE, na.rm=TRUE) # 11.83679
imp.long$ATE=unlist(tapply(imp.long$ATE, imp.long$.imp, function(x){y=x; if(sum(!is.na(x))>0){cut=quantile(x, 0.99); y[x>cut]=cut};return(y)}))
# 5.515691
imp1=as.mids(imp.long) ### imp1 has ATE added!

##averaged:
adm1=dm1 %>% mutate (ATE=case_when(trt==1~1/ps, trt==0~1/(1-ps)))
max(adm1$ATE)

source("ASMDCal.r") # This produces the connect-S plots.
#tavdiff=mydata

b1=b2=b1sd=b2sd=bdiff=NULL
a1=a2=a1sd=a2sd=adiff=NULL
### now for each imputed dataset
t1=Sys.time()
for (i in 1:imp$m){
  dm=complete(imp,i)
  dm$ps=ps.list[[i]]
  
  #### calculating  ATE 
  db=dm %>% mutate(ATE=case_when(trt==1~1/ps, trt==0~1/(1-ps)),
                   minps=pmin(ps, 1-ps),
                   mwt=case_when(trt==1~minps/ps, trt==0~minps/(1-ps))) %>% mutate_at(c(evtvars,"BLlivewithothers", "atrialFib", "BLdepress"),
                                                                                      function(x){as.numeric(as.character(x))})
  cut=quantile(db$ATE, 0.99)
  db$ATE[db$ATE>cut]=cut #r(dm1, ATE<cut)
  
  before=bal.stat(data = db,  vars = predvars, treat.var = "trt",   w.all = 1, sampw = 1,  estimand = "ATE",
                  get.means = TRUE, get.ks = TRUE,na.action = "level",multinom = FALSE)$results # unadjusted
  after=bal.stat(data = db,  vars = predvars, treat.var = "trt",   w.all = db$ATE, sampw = 1,  estimand = "ATE",
                 get.means = TRUE, get.ks = TRUE,na.action = "level",multinom = FALSE)$results # ate adjusted
  
  b1=cbind(b1, before$tx.mn)
  b2=cbind(b2, before$ct.mn)
  b1sd=cbind(b1sd, before$tx.sd)
  b2sd=cbind(b2sd, before$ct.sd)
  bdiff=cbind(bdiff, before$std.eff.sz)
  a1=cbind(a1, after$tx.mn)
  a2=cbind(a2, after$ct.mn)
  a1sd=cbind(a1sd, after$tx.sd)
  a2sd=cbind(a2sd, after$ct.sd)
  adiff=cbind(adiff, after$std.eff.sz)
}
t2=Sys.time()
t2-t1 # 1.1min

uwn=data.frame(arb=apply(b1,1,mean), arb.sd=apply(b1sd,1,mean), acei=apply(b2,1,mean), acei.sd=apply(b2sd,1,mean), smd=apply(bdiff,1,mean))
ate=data.frame(arb=apply(a1,1,mean),arb.sd=apply(a1sd, 1, mean), acei=apply(a2,1,mean),acei.sd=apply(a2sd, 1, mean), smd=apply(adiff,1,mean))
var=rownames(after)

## only doing this to get nicer labels
tp1=tp[covars %in% predvars]
varlabs1=varlabs[covars %in% predvars]
tt=myNOgrtab(xdat = db[,predvars],  typev=tp1, dnames=varlabs1, allSums = FALSE, annot = FALSE,  onelinebinary = TRUE, showbinarylevel = FALSE)
tlabs=tt[[1]][, c("Variable", "Levels")]
labs=prfmt(tt[[1]])$Variable
cbind(labs, var)
classvars=predvars[tp1>1];classvars

## replace raw var names with labeled ones
head(uwn)
tab0=cbind(tlabs,uwn) %>% mutate( prt1=case_when(arb<1 & acei<1 ~ paste(round(100*arb,1), "%", sep="" ),T~paste(round(arb,1),"(", round(arb.sd,1),")", sep="")),
                                  prt2=case_when(arb<1 & acei<1 ~ paste(round(100*acei,1), "%", sep="" ),T~paste(round(acei,1),"(", round(acei.sd,1),")", sep="")))%>% select(Variable, Levels, prt1, prt2, smd)

tab1=ate %>% mutate(prt1=case_when(arb<1 & acei<1 ~ paste(round(100*arb,1), "%", sep="" ),T~paste(round(arb,1),"(", round(arb.sd,1),")", sep="")),
                    prt2=case_when(arb<1 & acei<1 ~ paste(round(100*acei,1), "%", sep="" ),T~paste(round(acei,1),"(", round(acei.sd,1),")", sep=""))) %>% select(prt11=prt1, prt21=prt2, smd1=smd)

Table1=cbind(tab0, tab1) %>% mutate( std.eff.sz=abs(round(smd, 2)),std.eff.sz1=abs(round(smd1, 2)) ) %>% select(Variable, Levels, prt1, prt2, std.eff.sz, prt11, prt21, std.eff.sz1)
head(Table1) 
Table1=prfmt(Table1)
names(Table1)=c("Variable",  "ARB1", "ACEI1", "ASMD1",  "ARB2", "ACEI2", "ASMD2")
labelsp=paste(tlabs$Variable, tlabs$Levels)
labelsp
pos=length(var)+1-1:length(var)             
dplot=rbind(data.frame(vars1=labelsp, diff=abs(uwn$smd), stop.methods="Unadjusted", pos),
            data.frame(vars1=labelsp, diff=abs(ate$smd), stop.methods="ATE adjusted", pos))  %>% arrange(pos)
dplot$vars1=fct_inorder(dplot$vars1)

dplot$stop.methods=factor(dplot$stop.methods) %>% reord(c(2,1))
levels(dplot$stop.methods) =c("Before weighting", "After IP weighting")
head(dplot)      

bal1=ggplot(data=dplot, aes(vars1, diff, group="stop.methods"))+geom_point(aes(color=stop.methods))+geom_hline(yintercept = c(0.1, 0.2), color="red", lty=2)+xlab(" ")+coord_flip()
bal.plot1=bal1+scale_color_manual(values=c("black","red"),name="")+theme(legend.position = "bottom")+ylab("Absolute Standardized differences")
(fig2=bal.plot1)

png("Figure2_balance_plot_2022.png", res=300, width=2000, height=2200)
fig2
dev.off()

head(Table1)
printtab(Table1, "Table1stddiff2022.doc")


##&&&&&&&&&&&&&&&&&&&&&&fffff&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
##                 Table 3 - weighted
##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## ATE
## pd_mci_amnestic
fit=with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate 
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
## pool event rate ?
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_amnestic", timevar="pd_mci_amnestic_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_amnestic", datin=dm1) 
outp=data.frame(evt="pd_mci_amnestic", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr, HR)
outp
# pd_mci_protocol
fit=with(imp1, coxph(Surv(pd_mci_protocol_days, pd_mci_protocol) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_protocol", timevar="pd_mci_protocol_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_protocol", datin=dm1) 
outs1=data.frame(evt="pd_mci_protocol", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr, HR)
# pd
fit=with(imp1, coxph(Surv(pd_days, pd) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd", timevar="pd_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd", datin=dm1) 
outs2=data.frame(evt="pd", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# mci_amnestic
fit=with(imp1, coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_amnestic", timevar="mci_amnestic_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_amnestic", datin=dm1) 
outs3=data.frame(evt="mci_amnestic", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# mci_protocol
fit=with(imp1, coxph(Surv(mci_protocol_days, mci_protocol) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_protocol", timevar="mci_protocol_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_protocol", datin=dm1) 
outs4=data.frame(evt="mci_protocol", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# pd_mci_amnestic_death
fit=with(imp1, coxph(Surv(pd_mci_amnestic_death_days, pd_mci_amnestic_death) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_amnestic_death", timevar="pd_mci_amnestic_death_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_amnestic_death", datin=dm1) 
outs5=data.frame(evt="pd_mci_amnestic_death", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# pd_mci_protocol_death
fit=with(imp1, coxph(Surv(pd_mci_protocol_death_days, pd_mci_protocol_death) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_protocol_death", timevar="pd_mci_protocol_death_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_protocol_death", datin=dm1) 
outs6=data.frame(evt="pd_mci_protocol_death", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# pd_death
fit=with(imp1, coxph(Surv(pd_death_days, pd_death) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_death", timevar="pd_death_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_death", datin=dm1) 
outs7=data.frame(evt="pd_death", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# mci_amnestic_death
fit=with(imp1, coxph(Surv(mci_amnestic_death_days, mci_amnestic_death) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_amnestic_death", timevar="mci_amnestic_death_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_amnestic_death", datin=dm1) 
outs8=data.frame(evt="mci_amnestic_death", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# mci_protocol_death
fit=with(imp1, coxph(Surv(mci_protocol_death_days, mci_protocol_death) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_protocol_death", timevar="mci_protocol_death_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_protocol_death", datin=dm1) 
outs9=data.frame(evt="mci_protocol_death", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
# AnyDeathEvent
fit=with(imp1, coxph(Surv(death_days, AnyDeathEvent) ~ trt, robust=TRUE, weights=ATE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="AnyDeathEvent", timevar="death_days",  wt=db$ATE, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="AnyDeathEvent", datin=dm1) 
outs10=data.frame(evt="AnyDeathEvent", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr,HR)
out.ate=rbind(outp, outs1, outs2, outs3, outs4, outs5, outs6, outs7, outs8, outs9, outs10)

## OW
## pd_mci_amnestic
fit=with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate 
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
## pool event rate ?
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_amnestic", timevar="pd_mci_amnestic_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_amnestic", datin=dm1) 
outp=data.frame(evt="pd_mci_amnestic", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), hr, HR)
outp
# pd_mci_protocol
fit=with(imp1, coxph(Surv(pd_mci_protocol_days, pd_mci_protocol) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_protocol", timevar="pd_mci_protocol_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_protocol", datin=dm1) 
outs1=data.frame(evt="pd_mci_protocol", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr, HR)
# pd
fit=with(imp1, coxph(Surv(pd_days, pd) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd", timevar="pd_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd", datin=dm1) 
outs2=data.frame(evt="pd", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_amnestic
fit=with(imp1, coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_amnestic", timevar="mci_amnestic_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_amnestic", datin=dm1) 
outs3=data.frame(evt="mci_amnestic", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_protocol
fit=with(imp1, coxph(Surv(mci_protocol_days, mci_protocol) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_protocol", timevar="mci_protocol_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_protocol", datin=dm1) 
outs4=data.frame(evt="mci_protocol", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# pd_mci_amnestic_death
fit=with(imp1, coxph(Surv(pd_mci_amnestic_death_days, pd_mci_amnestic_death) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_amnestic_death", timevar="pd_mci_amnestic_death_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_amnestic_death", datin=dm1) 
outs5=data.frame(evt="pd_mci_amnestic_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# pd_mci_protocol_death
fit=with(imp1, coxph(Surv(pd_mci_protocol_death_days, pd_mci_protocol_death) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_protocol_death", timevar="pd_mci_protocol_death_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_protocol_death", datin=dm1) 
outs6=data.frame(evt="pd_mci_protocol_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# pd_death
fit=with(imp1, coxph(Surv(pd_death_days, pd_death) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_death", timevar="pd_death_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_death", datin=dm1) 
outs7=data.frame(evt="pd_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_amnestic_death
fit=with(imp1, coxph(Surv(mci_amnestic_death_days, mci_amnestic_death) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_amnestic_death", timevar="mci_amnestic_death_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_amnestic_death", datin=dm1) 
outs8=data.frame(evt="mci_amnestic_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_protocol_death
fit=with(imp1, coxph(Surv(mci_protocol_death_days, mci_protocol_death) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_protocol_death", timevar="mci_protocol_death_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_protocol_death", datin=dm1) 
outs9=data.frame(evt="mci_protocol_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# AnyDeathEvent
fit=with(imp1, coxph(Surv(death_days, AnyDeathEvent) ~ trt, robust=TRUE, weights=OW))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="AnyDeathEvent", timevar="death_days",  wt=db$OW, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="AnyDeathEvent", datin=dm1) 
outs10=data.frame(evt="AnyDeathEvent", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
out.ow=rbind(outp, outs1, outs2, outs3, outs4, outs5, outs6, outs7, outs8, outs9, outs10)

# Unweighted
## pd_mci_amnestic
fit=with(imp1, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)

## pool event rate ?
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_amnestic", timevar="pd_mci_amnestic_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_amnestic", datin=dm1) 
outp=data.frame(evt="pd_mci_amnestic", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
outp
# pd_mci_protocol
fit=with(imp1, coxph(Surv(pd_mci_protocol_days, pd_mci_protocol) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_protocol", timevar="pd_mci_protocol_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_protocol", datin=dm1) 
outs1=data.frame(evt="pd_mci_protocol", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# pd
fit=with(imp1, coxph(Surv(pd_days, pd) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd", timevar="pd_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd", datin=dm1) 
outs2=data.frame(evt="pd", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_amnestic
fit=with(imp1, coxph(Surv(mci_amnestic_days, mci_amnestic) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_amnestic", timevar="mci_amnestic_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_amnestic", datin=dm1) 
outs3=data.frame(evt="mci_amnestic", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_protocol
fit=with(imp1, coxph(Surv(mci_protocol_days, mci_protocol) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_protocol", timevar="mci_protocol_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_protocol", datin=dm1) 
outs4=data.frame(evt="mci_protocol", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# pd_mci_amnestic_death
fit=with(imp1, coxph(Surv(pd_mci_amnestic_death_days, pd_mci_amnestic_death) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_amnestic_death", timevar="pd_mci_amnestic_death_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_amnestic_death", datin=dm1) 
outs5=data.frame(evt="pd_mci_amnestic_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# pd_mci_protocol_death
fit=with(imp1, coxph(Surv(pd_mci_protocol_death_days, pd_mci_protocol_death) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_mci_protocol_death", timevar="pd_mci_protocol_death_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_mci_protocol_death", datin=dm1) 
outs6=data.frame(evt="pd_mci_protocol_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# pd_death
fit=with(imp1, coxph(Surv(pd_death_days, pd_death) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="pd_death", timevar="pd_death_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="pd_death", datin=dm1) 
outs7=data.frame(evt="pd_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_amnestic_death
fit=with(imp1, coxph(Surv(mci_amnestic_death_days, mci_amnestic_death) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_amnestic_death", timevar="mci_amnestic_death_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_amnestic_death", datin=dm1) 
outs8=data.frame(evt="mci_amnestic_death", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)
# mci_protocol_death
fit=with(imp1, coxph(Surv(mci_protocol_death_days, mci_protocol_death) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="mci_protocol_death", timevar="mci_protocol_death_days",  wt=1, datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="mci_protocol_death", datin=dm1) 
outs9=data.frame(evt="mci_protocol_death", n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2]%>% round(1), rate.arb=rates[2]%>% round(1), hr,HR)
# AnyDeathEvent
fit=with(imp1, coxph(Surv(death_days, AnyDeathEvent) ~ trt, robust=TRUE))
est=summary(pool(fit))$estimate
hr=exp(est);hr
se=summary(pool(fit))$std.error
HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt="AnyDeathEvent", timevar="death_days",  wt=1,  datin=db)})
rates=apply(rts, 1, mean)
num=numraw(evt="AnyDeathEvent", datin=dm1) 
outs10=data.frame(evt="AnyDeathEvent", n.acei=num[1], rate.acei=rates[1]%>% round(1), n.arb=num[2], rate.arb=rates[2]%>% round(1), hr,HR)

outuwn=rbind(outp, outs1, outs2, outs3, outs4, outs5, outs6, outs7, outs8, outs9, outs10)
outuwn.use=select(outuwn, rate.acei , rate.arb, hr,HR)

outcome.labs=c("Probable Dementia or amnestic MCI (censoring death)", 
               "Probable Dementia or protocol-defined MCI",
               "Probable Dementia alone",
               "Amnestic MCI alone",
               "Protocol-defined MCI alone",
               ## combining death
               "Probable Dementia or amnestic MCI or death",
               "Probable Dementia or protocol-defined MCI or death",
               "Probable Dementia or death",
               "Amnestic MCI or death",
               "Protocol-defined MCI or death",
               "Death")
Table3.point=cbind(Variables=outcome.labs, out.ate[,-1],  out.ow[,-c(1,2,4)], outuwn.use)
names(Table3.point)=c("Variables", "Num.ACEI", "Rate.ACEI.ate", "Num.ARB",  "Rate.ARB.ate", "hr.ate","HR.ate", 
                      "Rate.ACEI.ow", "Rate.ARB.ow", "hr.ow","HR.ow",
                      "Rate.ACEI.un",  "Rate.ARB.un","hr.un", "HR.un")
Table3.point

#dlab=Table3.point %>% select(Variables, HR.ate,HR.ow, HR.un)
#dlab$HR.ate=sub(",","-",dlab$HR.ate)
#dlab$HR.un=sub(",","-",dlab$HR.un)
#dlab$HR.ow=sub(",","-",dlab$HR.ow)
#dlab

#setdiff(ls(), obs1)

load("bt2022.RData") 
bt=bt2022
t0=bt$t0
sb=summary(bt)
pov=unique(substr(names(t0),1,4))
head(sb)
rownames(sb)=names(na.omit(t0))

## confint for Table 3:
dcis=sb[substr(rownames(sb),1,4)=="t3at",]
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
hrat=fmtci(exp(est), exp(lcl), exp(ucl), digits=2)

dcis=sb[substr(rownames(sb),1,4)=="t3ow",]
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
hrow=fmtci(exp(est), exp(lcl), exp(ucl), digits=2)

dcis=sb[substr(rownames(sb),1,4)=="t3un",]
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
hrun=fmtci(exp(est), exp(lcl), exp(ucl), digits=2)

Table3=Table3.point[,c("Variables", "Num.ACEI",  "Num.ARB", "Rate.ACEI.ate", "Rate.ARB.ate","HR.ate", 
                             "Rate.ACEI.ow", "Rate.ARB.ow", "HR.ow",
                             "Rate.ACEI.un",  "Rate.ARB.un", "HR.un")]
cbind(Table3$HR.ate, hrat)
Table3$HR.ate=hrat
Table3$HR.un=hrun
Table3$HR.ow=hrow
Table3
#printtab(Table3, "Table3new24.doc")
#printtab(data.frame(hrat, hrow, hrun), "BootCI_table3.doc")

table(dm1$evtcomp, dm1$trt)

## Source this for CIF plots ##
#source("getCIF.r")
## %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
## Subgroup Analysis ##

gett4=function(){rbind(
  cbind(group="Age<75",so2(filter(imp1,BLage<75),wt.type=wtuse)),
  cbind(group="Age>=75", so2(filter(imp1,BLage>=75),wt.type=wtuse)),
  
  cbind(group="Male", so2(filter(imp1,BLFemale==0),wt.type=wtuse)),
  cbind(group="Female",so2(filter(imp1,BLFemale==1),wt.type=wtuse)),
  
  cbind(group="Black",so2(filter(imp1,Black==1),wt.type=wtuse)),
  cbind(group="Not Black",so2(filter(imp1,Black==0),wt.type=wtuse)),
  
  cbind(group="Standard trt", so2(filter(imp1,IntensiveTrt==0),wt.type=wtuse)),
  cbind(group="Intensive trt", so2(filter(imp1,IntensiveTrt==1),wt.type=wtuse)),
  
  cbind(group="No MCI", so2(filter(imp1,BLMCI==0),wt.type=wtuse)),
  cbind(group="MCI", so2(filter(imp1,BLMCI==1),wt.type=wtuse)))}

wtuse=NULL
Table4.un=gett4()

  
wtuse="ATE"
(Table4.ate=gett4())
wtuse="OW"
Table4.ow=gett4()
Table4.ow


#Table5.un=
 gett5=function(){rbind(cbind(group="Overall",nso2(imp1,wt.type=wtuse)),
  cbind(group="Age<75",nso2(filter(imp1,BLage<75),wt.type=wtuse)),
  cbind(group="Age>=75", nso2(filter(imp1,BLage>=75),wt.type=wtuse)),
  
  cbind(group="Male", nso2(filter(imp1,BLFemale==0),wt.type=wtuse)),
  cbind(group="Female",nso2(filter(imp1,BLFemale==1),wt.type=wtuse)),
  
  cbind(group="Black",nso2(filter(imp1,Black==1),wt.type=wtuse)),
  cbind(group="Not Black",nso2(filter(imp1,Black==0),wt.type=wtuse)),
  
  cbind(group="Standard trt", nso2(filter(imp1,IntensiveTrt==0),wt.type=wtuse)),
  cbind(group="Intensive trt", nso2(filter(imp1,IntensiveTrt==1),wt.type=wtuse)),
  
  cbind(group="No MCI", nso2(filter(imp1,BLMCI==0),wt.type=wtuse)),
  cbind(group="MCI", nso2(filter(imp1,BLMCI==1),wt.type=wtuse)))}
 
wtuse=NULL
Table5.un=gett5()
wtuse="ATE"
Table5.ate=gett5()
wtuse="OW"
Table5.ow=gett5()


#load("boot2022suoershort.RData")
#length(btss$t0)
load("boott2022short.RData")
sbs=summary(bts);dim(sbs)
t0=bts$t0; length(t0)
names(t0)
nuni(names(t0))
dupout(names(t0))
rname=names(t0)
rname[substr(rname,1,4)=="AEs."]=paste(names(t0)[substr(rname,1,4)=="AEs."], 1:length(dupout(names(t0))))
rownames(sbs)=rname


#sbs=sb
## confint for Table 4:
dcis=sbs[substr(rownames(sbs),1,9) =="sgest.ate",]
sghrt=function(){
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
hrsg=fmtci(exp(est), exp(lcl), exp(ucl), digits=2)

#the Wald test will compare (log(HR2) - log(HR1))/sqrt(SE1^2 + SE2^2) to the standard normal distribution to obtain a pvalue. 
z1=abs(est[2]-est[1])/sqrt(se[1]**2+se[[2]]**2);
p1=(1-pnorm(z1))*2
z2=abs(est[4]-est[3])/sqrt(se[4]**2+se[[3]]**2);
p2=(1-pnorm(z2))*2
z3=abs(est[6]-est[5])/sqrt(se[6]**2+se[[5]]**2);
p3=(1-pnorm(z3))*2
z4=abs(est[8]-est[7])/sqrt(se[8]**2+se[[7]]**2);
p4=(1-pnorm(z4))*2
z5=abs(est[10]-est[9])/sqrt(se[10]**2+se[[9]]**2);
p5=(1-pnorm(z5))*2
p=sapply(c(p1,NA, p2,NA, p3,NA,p4,NA,p5,NA),fpvalue)
return(data.frame(HR=hrsg, p))}

new=sghrt()
cbind(Table4.ate, new)
Table4.ate=cbind(Table4.ate%>% select(-HR), new)
Table4.ate

dp=Table4.ate %>% mutate (est=as.numeric(substr(HR, 1,4)), 
                           lcl=as.numeric(substr(HR, 7,10)), ucl=as.numeric(substr(HR, 12,15)))
#write.csv(dp, "dpnew.csv")

### Want forest plot for this:
#dp=read.csv("dpforest.csv", na.strings=c("", " "))
dp=read.csv("dpnew.csv", na.strings=c("", " "))
head(dp)

unique(dp$group) %>% nchar

dplot=dp %>% mutate(aceirate=case_when(is.na(n.acei)~"", T~paste(n.acei, " (", round(rate.acei,1),")", sep="")),pval=sapply(as.numeric(dp$p.int), fpvalue),
                    arbrate=case_when(is.na(n.arb)~ " ", T~paste(n.arb,  " (", round(rate.arb,1), ")", sep="")),precesion=1/(ucl-lcl),
                    group=case_when(is.na(n.arb)|group=="Overall"~str_pad(group, 48, "right"), !is.na(n.arb)~str_pad(group, 60, "left"))) %>% select(group, aceirate, arbrate, est, lcl, ucl, hr, pval, precesion)


dplot$hr[is.na(dplot$hr)]=""
dplot$pval[dplot$pval=="NA"]=""
dplot$group[is.na(dplot$group)]=""
dplot$group
head(dplot);
head(dp)
myfunc=function(x){gsub( ","," - ", x)}
myfunc.vector=function(xv){sapply(xv, myfunc)}

dplot=dplot %>% mutate_at(vars(starts_with("hr")), myfunc.vector)

####   plot only   ###############
mytheme1=theme(panel.background = element_blank(),axis.text= element_blank(),axis.ticks = element_blank(), plot.margin = unit(c(0.2,0,0,0), "cm") , axis.line=element_blank(), legend.position="none")
#mytheme1=theme(panel.background = element_blank(), plot.margin = unit(c(0,0,0,0), "cm") , axis.line=element_blank(), legend.position="none")

dplot$ord=nrow(dplot)-1:nrow(dplot)+1
dplot$group
dplot
grlab=dplot$group
(index1=c(1, which(!is.na(grlab) & is.na(dplot$est))))
grlab[index1]
grlab[-index1]=""
grlab

dplot$ord
grlablevel=trim(dplot$group)
grlablevel[index1]=""

#expression(grlablevel)
cbind(dplot$ord, grlablevel)
dplot$ord[which(grepl(">=", grlablevel))] # 13
grlablevel[grepl(">=", grlablevel)]=""; #c(expression("">=75), expression("">=30), expression(">=4"))
grlablevel
#grlablevel
grlab_short=str_wrap(grlab, width=25)
dplot$ord
pl.lab=ggplot(data=dplot, aes(ord, y=0))+geom_text(aes(x=ord, y=0, label=grlab_short), hjust=0)+ylim(-0.53,1.5)+xlab("")+ylab(" ")+coord_flip() +ggtitle("Subgroup")
pl.lab
plab1=pl.lab+geom_text(aes(x=ord, y=1,label=grlablevel), hjust=0)+theme(plot.title =element_text(hjust=0.5))+xlim(-4,18)+geom_text(x=13,y=1,label=expression("">=75),hjust=0.1,size=3.6)
plab1=plab1+mytheme1 

phr=ggplot(data=dplot, aes(ord, est))+geom_point(aes(size=precesion), shape=15)+geom_errorbar(aes(ymin=lcl, ymax=ucl),width=0.4)+ylab("")+xlab("");phr
phr1=phr+scale_y_continuous(limits = c(0.2, 3), breaks=c(0.25, 0.5, 1, 2),trans = 'log10',expand = c(0, 0))+geom_segment(x=0, y=0, xend=25, yend=0,linetype=2)+coord_flip();
phr1
phr2=phr1+geom_segment(aes(x=-2.5, y=1.1, xend=-2.5, yend=1.75), arrow=arrow(length = unit(0.3,"cm")), size=0.5)+geom_vline(xintercept=0);phr2
phr3=phr2+geom_segment(aes(x=-2.5, y=0.9, xend=-2.5, yend=0.25), arrow=arrow(length = unit(0.3,"cm")), size=0.5)+mytheme1+ggtitle("Hazard Ratio")
phr4=phr3+geom_text(aes(x=-1, y=0.5, label="0.5"))+geom_text(aes(x=-1, y=1, label="1"))+geom_text(aes(x=-1, y=1.5, label="1.5"))+geom_text(aes(x=-1, y=2, label="2"))
phr4
(phr5=phr4+geom_text(aes(x=-4, y=0.4, label="Favors ARB Initiators"))+geom_text(aes(x=-4, y=1.6, label="Favors ACEI Initiators")))#+geom_text(aes(x=-1, y=0, label="0")))
phr6=phr5+theme(plot.title =element_text(hjust=0.9))+xlim(-4,18);phr6

phr2=phr1+theme( legend.position="none", axis.text.y = element_blank(), plot.title =element_text(hjust=1), axis.text.x= element_text(size=12))
phr2=phr2+ggtitle("Hazard Ratio")
phr1+geom_segment(aes(x=-4, y=0, xend=-4, yend=1))
hrlab=dplot$hr
pl.lab=ggplot(data=dplot, aes(ord, y=0))+geom_text(aes(x=ord, y=0,label=hrlab), hjust=0)+ylim(0,0.5)+xlab("")+ylab(" ")+coord_flip() +ggtitle("(95% CI)")
plab2=pl.lab+mytheme1 +theme(plot.title =element_text(hjust=-0.25),  axis.ticks = element_blank())+xlim(-4,18)
plab2

lab=dplot$pval
pr=ggplot(data=dplot, aes(ord, y=0))+geom_text(aes(x=ord, y=0.15,label=lab), hjust=1)+ylim(0,0.5)+xlab("")+ylab(" ")+coord_flip() 
pr3=pr+theme(panel.background = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.2,0,0,0),"cm"), plot.title =element_text(hjust=-0.5))
pr3=pr3+ggtitle("    P value", subtitle="for interaction                      ")+xlim(-4,18) +theme(plot.subtitle =element_text(hjust=0))
pr3

lab1=dplot$aceirate
pr1=ggplot(data=dplot, aes(ord, y=0))+geom_text(aes(x=ord, y=0.015,label=lab1), hjust=0.5)+ylim(0.01,0.02)+xlab("")+ylab(" ")+coord_flip()+theme(panel.background = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.2,0,0,0), "cm") , )
pr1=pr1+ggtitle("ACEI ", subtitle="per 100 person years)      ")+theme(plot.title =element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.5))+xlim(-4,18)
pr1

lab2=dplot$arbrate
pr=ggplot(data=dplot, aes(ord, y=0))+geom_text(aes(x=ord, y=0.015,label=lab2), hjust=0.5)+ylim(0.01,0.02)+xlab("")+ylab(" ")+coord_flip() 
pr2=pr+theme(panel.background = element_blank(),axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.2,0,0,0), "cm") )
pr2=pr2+ggtitle("ARB", subtitle="#Events (#events ")+theme(plot.title =element_text(hjust=0.5), plot.subtitle=element_text(hjust=0.4))+xlim(-4,18)
pr2

#png("Forest_interation_2022_.png", width=2500, height=1100, res=160)
ggarrange(plab1, pr2, pr1, phr6, plab2, pr3, ncol=6, align="hv", widths = c(3.3,1.3,1.35,3.4,2,1.3)  )
#dev.off()

ggsave('forest_log.pdf', width=14, height=7)

dcis=sbs[substr(rownames(sbs),1,8) =="sgest.ow",]
new=sghrt()
Table4.ow=cbind(Table4.ow%>% select(-HR), new)
Table4.ow

dcis=sbs[substr(rownames(sbs),1,8) =="sgest.un",]
new=sghrt()
Table4.un=cbind(Table4.un%>% select(-HR), new)
Table4.un


dcis=sbs[substr(rownames(sbs),1,6) =="nsg.un",]
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
(hrsg=fmtci(exp(est), exp(lcl), exp(ucl), digits=2))
cbind(hrsg, Table5.un$HR)
est=est[-1]
#the Wald test will compare (log(HR2) - log(HR1))/sqrt(SE1^2 + SE2^2) to the standard normal distribution to obtain a pvalue. 
z1=abs(est[2]-est[1])/sqrt(se[1]**2+se[[2]]**2);
p1=(1-pnorm(z1))*2
z2=abs(est[4]-est[3])/sqrt(se[4]**2+se[[3]]**2);
p2=(1-pnorm(z2))*2
z3=abs(est[6]-est[5])/sqrt(se[6]**2+se[[5]]**2);
p3=(1-pnorm(z3))*2
z4=abs(est[8]-est[7])/sqrt(se[8]**2+se[[7]]**2);
p4=(1-pnorm(z4))*2
z5=abs(est[10]-est[9])/sqrt(se[10]**2+se[[9]]**2);
p5=(1-pnorm(z5))*2
p=sapply(c(NA, p1,NA, p2,NA, p3,NA,p4,NA,p5,NA),fpvalue);
#printtab(data.frame(HR=hrsg, p),"Tab5un.doc")
Table5.un$HR=hrsg
Table5.un$p=p
Table5.un

dcis=sbs[substr(rownames(sbs),1,7) =="nsg.ate",]
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
(hrsg=fmtci(exp(est), exp(lcl), exp(ucl), digits=2))
cbind(hrsg, Table5.ate$HR)
est=est[-1]
#the Wald test will compare (log(HR2) - log(HR1))/sqrt(SE1^2 + SE2^2) to the standard normal distribution to obtain a pvalue. 
z1=abs(est[2]-est[1])/sqrt(se[1]**2+se[[2]]**2);
p1=(1-pnorm(z1))*2
z2=abs(est[4]-est[3])/sqrt(se[4]**2+se[[3]]**2);
p2=(1-pnorm(z2))*2
z3=abs(est[6]-est[5])/sqrt(se[6]**2+se[[5]]**2);
p3=(1-pnorm(z3))*2
z4=abs(est[8]-est[7])/sqrt(se[8]**2+se[[7]]**2);
p4=(1-pnorm(z4))*2
z5=abs(est[10]-est[9])/sqrt(se[10]**2+se[[9]]**2);
p5=(1-pnorm(z5))*2
p=sapply(c(NA, p1,NA, p2,NA, p3,NA,p4,NA,p5,NA),fpvalue);
Table5.ate$HR=hrsg
Table5.ate$p=p
Table5.ate
printtab(Table5.ate,"Tab5ate.doc")


dcis=sbs[substr(rownames(sbs),1,6) =="nsg.ow",]
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
(hrsg=fmtci(exp(est), exp(lcl), exp(ucl), digits=2))
cbind(hrsg, Table5.ow$HR)
est=est[-1]
#the Wald test will compare (log(HR2) - log(HR1))/sqrt(SE1^2 + SE2^2) to the standard normal distribution to obtain a pvalue. 
z1=abs(est[2]-est[1])/sqrt(se[1]**2+se[[2]]**2);
p1=(1-pnorm(z1))*2
z2=abs(est[4]-est[3])/sqrt(se[4]**2+se[[3]]**2);
p2=(1-pnorm(z2))*2
z3=abs(est[6]-est[5])/sqrt(se[6]**2+se[[5]]**2);
p3=(1-pnorm(z3))*2
z4=abs(est[8]-est[7])/sqrt(se[8]**2+se[[7]]**2);
p4=(1-pnorm(z4))*2
z5=abs(est[10]-est[9])/sqrt(se[10]**2+se[[9]]**2);
p5=(1-pnorm(z5))*2
p=sapply(c(NA, p1,NA, p2,NA, p3,NA,p4,NA,p5,NA),fpvalue);

Table5.ow$HR=hrsg
Table5.ow$p=p
Table5.ow

printtab(Table5.ow, "Tab5ow.doc")

Tab4names=c("Subgroup","#Events ACEI","Event rate ACEI","#Events ARB","Event rate ARB","Hazard rario", "p")
names(Table4.un)=names(Table4.ow)=names(Table4.ate)=Tab4names


###  AE
###   Table 3 - weighted
## pd_mci_amnestic
#evt="hypotenSAEEvent";timevar="t_hypotenSAEEvent";
weight=NULL
aept=function(evt, timevar,wt.type=weight){
 fm=paste("Surv(", timevar, ",", evt, ")~trt", sep="");fm
 if(is.null(wt.type)){ fit=with(imp1, coxph(eval(parse(text=fm))))
 rts=sapply(1:imp1$m, function(s){db=complete(imp1, s) %>% mutate(fwt=1);  evtrate(evt, timevar,  wt=db$fwt, datin=db)})}else{
  fit=with(imp1, coxph(eval(parse(text=fm)), robust=TRUE, weights=eval(parse(text=wt.type))))
  rts=sapply(1:imp1$m, function(s){db=complete(imp1, s);  evtrate(evt, timevar,  wt=db[,wt.type], datin=db)})}
  est=summary(pool(fit))$estimate
  hr=exp(est);hr
  rates=apply(rts, 1, mean)
  num=numraw(evt=evt, datin=dboot) 
  outp=data.frame(evt, n.acei=num[1], rate.acei=rates[1], n.arb=num[2], rate.arb=rates[2], hr)
  return(outp)}

getaehr=function(){
  ae=aept(evt="AnySAEEvent", timevar="t_AnySAEEvent");
  ae1=aept(evt="hypotenSAEEvent" , timevar="t_hypotenSAEEvent") ;ae1## warnings
  ae2=aept(evt="SyncopeSAEEvent" , timevar="t_SyncopeSAEEvent") ## warnings
  ae3=aept(evt="BradySAEEvent" , timevar="t_BradySAEEvent") # warnings
  ae4=aept(evt="electroSAEEvent" , timevar="t_electroSAEEvent")
  ae5=aept(evt="InjFallSAEEvent" , timevar="t_InjFallSAEEvent")
  ae6=aept(evt="AKISAEEvent" , timevar="t_AKISAEEvent" )
  ae7=aept(evt="hypotenallEvent" , timevar="t_hypotenallEvent") 
  ae8=aept(evt="SyncopeallEvent" , timevar="t_SyncopeallEvent") 
  ae9=aept(evt="BradyallEvent" , timevar="t_BradyallEvent")
  #table(db$trt, db$BradyallEvent)
  ae10=aept(evt="electroallEvent" , timevar="t_electroallEvent")
  ae11=aept(evt="InjFallallEvent" , timevar="t_InjFallallEvent")
  ae12=aept(evt="AKIallEvent" , timevar="t_AKIallEvent")
  ae13=aept(evt="LowNaEvent" , timevar="t_LowNaEvent")
  #ae14=NA;#jjlab2=jjlab1 %>% mutate(event="High",Outcome="Serum sodium >150 mmol/liter")
  ae15=aept(evt="LowKEvent" , timevar="t_LowKEvent")
  ae16=aept(evt="HighKEvent" , timevar="t_HighKEvent")
  ae17=aept(evt="OHmeasuredEvent" , timevar="t_LowKEvent" )
  ae18=aept(evt="OHbothEvent" , timevar="t_OHbothEvent")
  return(rbind(ae,  ae1, ae2, ae3, ae4,  ae6,  ae7,ae8, ae9, ae10,  ae12, ae13, ae15, ae16, ae17, ae18) %>% mutate(rate.acei=round(rate.acei,1),
                                                                                                                        rate.arb=round(rate.arb,1)))}

Tables1.un=getaehr()
Tables1.un
weight="ATE"
Tables1.ate=getaehr()
weight="OW"
Tables1.ow=getaehr()
head(Tables1.ow)

Tables1.point=cbind(Tables1.un[, c(1,2,4,3,5,6)], Tables1.ate[, c(3,5,6)], Tables1.ow[, c(3,5,6)])

## confint for Table 3:
dcis=sb[substr(rownames(sb),1,3)=="AEs",]
dim(dcis)
est=dcis$original
se=dcis$bootSE
lcl=est-1.96*se
ucl=est+1.96*se
hr=fmtci(exp(est), exp(lcl), exp(ucl), digits=2)

hr[1:16]


names(Tables1.point)=c("Adverse event", "#ACET","#ARB","rate1", "rate2", "hrun", "rate1a", "rate2a", "hrate","rate1o", "rate2o", "hrow")
dim(Tables1.point)
Tables1.point
Tables1=mutate(Tables1.point, hrun=hr[1:16], hrate=hr[17:32], hrow=hr[33:48]) #%>% select(-hr) 
Tables1
printtab(Tables1, "Tables124.doc")


