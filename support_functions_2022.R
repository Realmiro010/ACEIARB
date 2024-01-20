
## plot histogram ps for 2 groups
hist.ps=function(data, nbins=50){ggplot(data,aes(ps))+geom_histogram(alpha=0.5,bins =nbins)+ylab(" ")+facet_wrap(~group,ncol=1)+xlab("Propensity score")}

numraw=function(evt, datin){tapply(datin[,evt], datin$trt, sum, na.rm=TRUE) }

## weighted event rate
evtrate=function(evt, timevar, datin, wt){
  dat=datin
  dat$dwt=wt
  des=svydesign(ids=~1, weights=~dwt, data=dat) 
  
  dat$Event=dat[,evt]*wt
  dat$Year=dat[, timevar]*wt/365
  dm=na.omit(select(dat, Event, Year, trt))
  
  num=tapply(dm$Event, dm$trt, sum) 
  denom=tapply(dm$Year, dm$trt, sum)
  #nraw=paste(num, "/" ,round(denom,1), sep="")
  rate=(num/denom)*100
  return(rate)}

### hr point estimate - weighted
svyhr=function(evt, timevar, datin, wt){
  
  dat=datin
  dat$dwt=wt
  des=svydesign(ids=~1, weights=~dwt, data=dat) 
  
  fm=as.formula(paste("Surv(", timevar, ",", evt, ")~trt") );  
  fit=svycoxph(fm, design=des)
  hr=exp(coef(fit)) 
  return(hr)
}

#evt="pd_mci_amnestic"; timevar="pd_mci_amnestic_days";varname="Probable Dementia or amnestic MCI (censoring death)";datin=db # wt=dm1$mwt; datin=dm1

getwt=function(evt, timevar, datin, varname, wt){
  dat=datin
  dat$dwt=wt
  des=svydesign(ids=~1, weights=~dwt, data=dat) 
  
  fm=as.formula(paste("Surv(", timevar, ",", evt, ")~trt") );fm
  fit=svycoxph(fm, design=des)
  hr=exp(coef(fit)) 
  hrci=exp(confint(fit))
  HR=fmtci(hr, hrci[,1], hrci[,2], digits=2);

  # Unweighted #events
num.raw=tapply(dat[,evt], dat$trt, sum, na.rm=TRUE) 
    # Weighted event rate
  dat$Event=dat[,evt]*wt
  dat$Year=dat[, timevar]*wt/365
  dm=na.omit(select(dat, Event, Year, trt))
  
  num=tapply(dm$Event, dm$trt, sum) 
  denom=tapply(dm$Year, dm$trt, sum)
  #nraw=paste(num, "/" ,round(denom,1), sep="")
 rate=(num/denom)*100
  out=data.frame(evt, varname, num.raw[1], rate[1], num.raw[2], rate[2], HR);out 
  #Rate=paste(num.raw, "(",round(rate,1),")", sep="")
  #out=data.frame(t(c(evt, varname,Rate, HR, fpvalue(p))))
  names(out)=c("event","Outcome",  "Events1","rateACEI", "Events2","rateARB", "HR (95% CI)")

  return(out)}






 
#evtvar="evtcomp"; dayvar="comp_days"
#dat=complete(dm1)
#tx=1
#getcif(dat=dat,tx=1, evtvar=evtvar, dayva=dayvar )
#raw.cif0=getcif(dat=dm1, tx=0, evtvar="evtcomp", dayvar="comp_days", weighted=FALSE)%>% as.data.frame %>% mutate(Year=tm1/365)

#dat=dm1;tx=0
#evtvar="evtcomp"; dayvar="comp_days"; wt.type=NULL

getcif = function(dat, tx, evtvar, dayvar, wt.type=NULL) {
  dat = dat[!is.na(dat[, dayvar]) & dat[, dayvar]>0,]
  tmpdat = dat[dat$trt==tx,]
 
if(is.null(wt.type)){wt=rep(1,nrow(tmpdat))}else{
if(wt.type=="ATE"){wt=tmpdat$ATE}
if(wt.type=="OW"){wt=tmpdat$OW}}

  tmpdat%>% select(evtcomp, comp_days) %>% arrange(comp_days) %>% filter(evtcomp>0) %>%  head
 
  tm1 = sort(unique(tmpdat[,dayvar][tmpdat[,evtvar]==1]));tm1
  ntm1 = length(tm1) #122
  tm12 = sort(unique(tmpdat[,dayvar][tmpdat[, evtvar]!=0]))
  ntm12 = length(tm12) #141
  
  atrisk1 = rep(0,ntm1)
  nj1 = rep(0,ntm1)
  haz1 = rep(0,ntm1)
  

  for (j in 1:ntm1){
    nj1[j] = sum((tmpdat[, dayvar]==tm1[j] & tmpdat[,evtvar]==1)*wt)
    atrisk1[j] = sum((tmpdat[,dayvar]>=tm1[j])*wt)
    haz1[j] = nj1[j]/atrisk1[j]
  }
  
  atrisk12 = rep(0,ntm12)
  nj12 = rep(0,ntm12)
  haz12 = rep(0,ntm12)
  s12 = rep(0,ntm12)
  
  for (j in 1:ntm12){
    nj12[j] = sum((tmpdat[,dayvar]==tm12[j] & tmpdat[,evtvar]!=0)*wt)
    atrisk12[j] = sum((tmpdat[,dayvar]>=tm12[j])*wt)
    haz12[j] = nj12[j]/atrisk12[j]
  }
  
  cbind(nj12, atrisk12, haz12) %>% head
  cbind(nj1, atrisk1, haz1) %>% head
  
  s12 = c(1,cumprod(1-haz12))
  tm12 = c(0,tm12)
  
  s1 = rep(0,ntm1)
  for (j in 1:ntm1) {
   # which(tm12<tm1[j])
  s1[j] = s12[max(which(tm12<tm1[j]))]
  }
  
cif = cumsum(s1*haz1)
out=cbind(tm1,cif)
if(max(tm1)<6){ends=cbind(tm1=6, cif=max(cif)); out=cbind(out, ends)}
  return(out)
}

#get KM estimate given weight
wtkm=function(dat, evtvar="AnyDeathEvent", dayvar="death_days", wt="fwt"){
  fm=as.formula(paste("Surv(", dayvar, ",", evtvar, ")~group") );fm
  des=svydesign(ids=~1, weights=~eval(parse(text=wt)), data=dat) 
  s1<-svykm(fm, design=des)
  (trt0=s1[[1]])
  Year=trt0$time/365
  cif=1-trt0$surv
  ev0=data.frame(Year, cif) %>% mutate(group=names(s1)[1])
  trt1=s1[[2]]
  Year=trt1$time/365
  cif=1-trt1$surv
  ev1=data.frame(Year, cif) %>% mutate(group=names(s1)[2])
  cif.death.raw=rbind(ev0,ev1) 
  return(cif.death.raw)}



#getcif(ddat, tx=1, evtvar="evtcomp", dayvar="comp_days", wt.type = "OW")
cumcif=function(evtvar, dayvar, wt.type="ATE"){
  cif1=NULL
  cif0=NULL
  for (i in 1:10){
    dat = imp.long[imp.long$.imp==i, ]
    cif1 =cbind(cif1, getcif(dat=dat,tx=1, evtvar=evtvar, dayva=dayvar, wt.type=wt.type ))
    cif0= cbind(cif0, getcif(dat=dat, tx=0, evtvar=evtvar, dayvar=dayvar, wt.type=wt.type))
  }
  mean.cif1= data.frame(Year=cif1[,1]/365, cif=apply(cif1[, 2*(1:10)],1,mean))
  mean.cif0= data.frame(Year=cif0[,1]/365, cif=apply(cif0[, 2*(1:10)],1,mean))

  cif.e=rbind(mean.cif0 %>% mutate(group="ACEI"),mean.cif1 %>% mutate(group="ARB"))
  return(cif.e)
}


so2=function(imps,wt.type=NULL){ 
  if(is.null(wt.type)) {fit=with(imps, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt))}else{
  fit=with(imps, coxph(Surv(pd_mci_amnestic_days, pd_mci_amnestic) ~ trt, robust=TRUE, weights=eval(parse(text=wt.type))))}
  hr=exp(summary(pool(fit))$estimate);
  est=summary(pool(fit))$estimate
  hr=exp(est);hr
  se=summary(pool(fit))$std.error
  HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2)
  ## pool event rate ?
  if(is.null(wt.type)){rts=sapply(1:imps$m, function(s){db=complete(imps, s);  evtrate(evt="pd_mci_amnestic", timevar="pd_mci_amnestic_days", wt=1, datin=db)})}else{
     rts=sapply(1:imps$m, function(s){db=complete(imps, s);  evtrate(evt="pd_mci_amnestic", timevar="pd_mci_amnestic_days",  wt=db[,wt.type], datin=db)})}
rates=apply(rts, 1, mean)
  num=numraw(evt="pd_mci_amnestic", datin=complete(imps)) 
  outp=data.frame(n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), HR)
  outp }

#so2(imps=filter(imp1,BLage<75), wt.type=NULL)
  #neg1=aehr(evt="CompositeEvent" , timevar="t_CompositeEvent")

nso2=function(imps,wt.type=NULL){
  #imps= filter(imp1,BLage<75)
  if(is.null(wt.type)) {fit=with(imps, coxph(Surv(t_CompositeEvent, CompositeEvent) ~ trt))}else{
    fit=with(imps, coxph(Surv(t_CompositeEvent, CompositeEvent) ~ trt, robust=TRUE, weights=eval(parse(text=wt.type))))}
  hr=exp(summary(pool(fit))$estimate);
  est=summary(pool(fit))$estimate
  hr=exp(est);hr
  se=summary(pool(fit))$std.error
  (HR=fmtci(hr, exp(est-1.96*se), exp(est+1.96*se), digits=2))
  ## pool event rate ?
  if(is.null(wt.type)){rts=sapply(1:imps$m, function(s){db=complete(imps, s);  evtrate(evt="CompositeEvent", timevar="t_CompositeEvent", wt=1, datin=db)})}else{
    rts=sapply(1:imps$m, function(s){db=complete(imps, s);  evtrate(evt="CompositeEvent", timevar="pd_mci_amnestic_days",  wt=db[,wt.type], datin=db)})}
  rates=apply(rts, 1, mean)
  num=numraw(evt="CompositeEvent", datin=complete(imps)) 
  outp=data.frame(n.acei=num[1], rate.acei=rates[1]%>%round(1), n.arb=num[2], rate.arb=rates[2]%>%round(1), HR)
  outp
  }


######## the asd function to calculate the balance##########################################################################
abs_stand_diff <- function(x_j, z, w){
  # Inputs:
  # x_j: a vecor of the covariate
  # z: a vector of treatment indicator
  # w: a vector of weights
  if (anyNA(w)) {return (NA)} else
    x_j <- as.numeric(x_j)
  
  absdiff <- abs(sum(x_j*z*w)/sum(z*w) - sum(x_j*(1-z)*w)/sum((1-z)*w))
  tstat <- absdiff/sqrt((var(x_j[which(z==1)])+var(x_j[which(z==0)]))/2)
  return (tstat)
}




  