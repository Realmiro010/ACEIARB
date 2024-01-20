
## other functions
npf=function(x){num=sum(x, na.rm=TRUE); denom=length(na.omit(x)); pct=100*num/denom; return(paste(num, " (", round(pct,1), "%)", sep="")) }

todate=function(x){as.Date(as.character(x), format="%m/%d/%Y")} # updated on 4/16/2020, however this still requires all elements to be in the same format 

todatey=function(x)
{y=as.Date(as.character(x), "%m/%d/%y")
z=as.Date(ifelse(y > Sys.Date(), format(y, "19%y-%m-%d"), format(y)))
return(z)}
todatey("11/8/17")

## string functions  ##;
findchar=function(char, allchar){allchar[grep(char, allchar, ignore.case=TRUE)]}
cutoff=function(xvec, string, pos){yy=sapply(xvec, function(s){unlist(strsplit(s,string))[pos]}); return(as.character(yy))};
#xc=c("tyu:ff", "s8:40"); 
# cutoff(xc, ":",1)
# "tyu" "s8" 

Sep=function(xvector, sep, index){
sepfun=function(x){unlist(strsplit(x, sep))[index]}
sapply(xvector, sepfun)}

#XX1=c("ghgh_45 ", "hh_34","fhgh_2003 23")
#XX=c("ghgh 45", "hh 34","fhgh 2003 23")
#Sep(XX,sep=" ",index=1)
#Sep(XX1, "_", 2)
## fill group 
fill=function(x){if(nuni(na.omit(x))>1) {stop("More than one value found for the same id!")}else if(nuni(na.omit(x))==1) {x[is.na(x)]=unique(na.omit(x))}; return(x)}
fillgroup=function(x, id){unlist(tapply(x, id, fill))}
#A=rep(1:4, each=6)
#B=c(rep(NA, 6), rep(c(NA,2,NA,NA, NA, NA), 3))
#C=c(unlist(tapply(B, A, fill)))
#D=fillgroup(B,A)
#cbind(A,B,C,D)


upall=function(x){h=unlist(strsplit(x," ")); h1=toprop(h); return(paste(h1, collapse=" "))}


### get commen elements of a few vectors (numeric or strings..), L is a list;
comout=function(L)
{len=length(L);
getc=unlist(L[1]); 
for (k in 2:len) {getc=intersect(getc, unlist(L[k]))}; 
return(getc)}
str1=1:7;str2=4:9;str3=6:8;
comout(list(str1,str2,str3))

findchars=function(strs, where){
alls=lapply(strs,function(x){findchar(x,where)});
comout(alls)}
## find strings contain all three of "kill", "jj", "wave" in strings.
strings=c("kill_wave_456_jj","wave_ty_kill_jj","ghjjkl_wavekl_wave", "kill wave", "ghjj")
findchars(c("kill","jj","wave"), strings)

trim.leading <- function (x)  sub("^\\s+", "", x)
trim.trailing <- function (x) sub("\\s+$", "", x);
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
#trim.leading("   fff")
#trim("  gh ff   ") 
# similar to sas function "propcase"
toprp=function(x){xx=trim(x);st1=toupper(substr(xx,1,1)); fol=tolower(substr(xx,2,nchar(xx))); y=paste(st1,fol,sep=""); if(is.na(x)){y=NA}; return(y)}
toprop=function(x){sapply(x, toprp)}

### Function to output one table to word
#printtab(dat, "xxx.doc")
printtab=function(dat,filename)
{rtf<-RTF(filename,width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
addParagraph(rtf, "" )
addTable(rtf,dat,font.size=9,col.justify='R',header.col.justify='R',row.names=FALSE,NA.string="-")
addNewLine(rtf)
rtf::done(rtf)}

## Function to reorder the levels of a factor;
reord=function(Fac, neword)
{X=factor(Fac);
Y=factor(X, levels(X)[neword])
return(Y)}
#d$Grade=reord(d$Grade, c(2,3,1,4))

### Function to get class of each variable in a data frame (similar to str)
vclas=function(datframe){sapply(datframe,class)}
## get the number of unique values of a vector
nuni=function(x){length(unique(x))}; ### function to count number of unique values


tablem=function(...){table(..., useNA="always")}

dupout=function(x){y=x[x %in% x[duplicated(x)]]; z=y[order(y)]; return(z)};

### take the nonmissing of two columns (if at least one is missing):
anyva=function(x,y){x[is.na(x)]=y[is.na(x)]; return(x)}## function to take the non missing of the two.

## rename ##;
rename=function(dat, old, new){names(dat)[names(dat)==old]=new; return(names(dat))}
## Count number of missings;

### find the variables in a dataframe with missingness and the number of missings.
findmiss=function(dat){M=apply(dat,2,function(x){sum(is.na(x))}); M[M>0]};
### get sample sizes for each variable in a dataframe
countn=function(dat){M=apply(dat,2,function(x){sum(!is.na(x))}); return(M)};

## format estimates and confidence interval:
#fmtci=function(est, cil, ciu, digits){paste(round(est, digits),"(",round(cil,digits),"~", round(ciu,digits),")",sep="")}
fmtci=function(est, cil, ciu, digits=NULL){
out=paste(est," (",cil,",", ciu,")",sep="");
if (!is.null(digits)){
ins=paste("%.",digits,"f",sep="");
out=paste( sprintf(ins, round(est,digits))," (",sprintf(ins, round(cil,digits)),",", sprintf(ins, round(ciu,digits)),")",sep="")} 
return(out)}




#Function to calculate relative percent change from exponentiated coefficients (RR or OR)
#it appears this is what Polina reports rather than relative risk or relative rate ratios.
frrp = function(x,fmt=TRUE){
if(x<1){
nx = -(1-x)*100
}else if(x>=1){
nx = (x-1)*100
}#end if
out=paste(round(nx,1),"%",sep="");
if (fmt==FALSE) {out=nx}
return(out)
}#end frrp

## count number of obs and number of subjects., check if there is any missing in subject id.
getns=function(dat, id) {c(nrow(dat), nuni(dat[,id]), sum(is.na(dat[,id])))}
guess=function(dat, nlevel){apply(dat,2, function(x){as.numeric(nuni(x)<=nlevel)+1})};


################################################################
library(gridExtra);
library(grid);
library(lattice)

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}




narm=function(dat){dnew=dat[apply(dat,1,function(x){sum(!is.na(x))>0}),]; return(dnew)} #remove rows missing all from dataset
retain=function(x){y=x; if (length(y)>=2){ for (i in 2:length(y)){if (is.na(y[i])) { y[i]=y[i-1]}}}; return(y)} # fill in blanks using previous value (for a vector)

naton=function(x){y=x;y[is.na(x)]=0; return(y)}; # assign "0" or FALSE to missing values.

### make column numbers for excel..
alp=c("abcdefghijklmnopqrstuvwxyz")
colnum=1:26;
for (i in 1:26){colnum[i]=substr(alp, i,i)};
b1=rep(c("", colnum[1:26]), each=26);b1
b2=rep(colnum, 26);
alps=paste(b1[1:350], b2[1:350], sep="")

#####  RTF related  #######
## to bold character 
bold=function(x){
y=paste( "{\\b ", x, "}", sep="" );
return(y)}
# to italic character
italic=function(x){
y=paste( "{\\i ", x, "}", sep="" );
return(y)}

### superscript:"0.05{\\super 2}\n")
super=function(x, s){
y=paste( x, "{\\super ", s, "}", sep="" );
return(y)}
super(0.05, 2)
super("g","5");
bold("g")

setcomp=function(set1, set2)
{if(!setequal(set1, set2)){numout=c(nuni(c(set1, set2)),length(setdiff(set1, set2)), length(setdiff(set2, set1))); 
return(numout)}else {print("Set equal"); print(nuni(set1))}
}
#a=1:99; b=1:99 ;e=8:111 ;d=c(NA, 8:111)
#setcomp(a,d); setcomp(a,e)

### From Chelsea 
plot2.cat <- function(value,catV,catH) {
  op <- par('mar','fig','mgp','cex.axis','las','srt')
  on.exit(par(op))
  if(class(catV) != "factor") break
  Hn <- dim(catH)[2]
  info <- split(cbind(catH,value),catV)
  Vn <- length(info)
  xmax <- 0
  for(i in 1:Vn) {
    info[[i]] <- append(list(values=info[[i]]),hist(info[[i]][,Hn+1],breaks=10,plot=F))
    xmax <- max(xmax,info[[i]]$counts)
  }
  h <- 7/Hn*Vn+1
  coord <- list(LR=seq(0.5/7.5,1,length=Hn+2),BT=seq(1-0.5/h,0.5/h,length=Vn+1))
  for (j in 1:Hn) {
    par(mar=c(0,0,0,0),fig=c(coord$LR[j+1],coord$LR[j+2],coord$BT[1],1),new=(j!=1))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    mtext(names(info[[1]]$values)[j],1,-1,cex=0.75)
  }
  for (i in 1:Vn){
    par(mar=c(0,0,0,0),fig=c(0,coord$LR[1],coord$BT[i+1],coord$BT[i]),las=0,new=T)
    plot(c(0, 2), c(0, 2), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    mtext(names(info)[i],2,-2,cex=0.75)
    par(mar=c(0.5,0.25,0.5,0.25),fig=c(coord$LR[1],coord$LR[2],coord$BT[i+1],coord$BT[i]),
        new=T,las=1,mgp=c(0.5,0.25,0),tck=-0.025,cex.axis=0.5)
    if(i==Vn) {
      plot(NULL, type = "n", xlim = c(0, xmax), ylim = range(info[[i]]$breaks),xlab="Count",ylab="")
      rect(0, info[[i]]$breaks[1:(length(info[[i]]$breaks) - 1)], info[[i]]$counts, info[[i]]$breaks[2:length(info[[i]]$breaks)])
      for (j in 1:Hn) {
        par(mar=c(0.5,0.25,0.5,0.25),fig=c(coord$LR[j+1],coord$LR[j+2],coord$BT[i+1],coord$BT[i]),
            new=T,las=2,mgp=c(0.5,0.25,0),tck=0,cex.axis=0.5)
        boxplot(info[[i]]$values[,Hn+1] ~ info[[i]]$values[,j],ylim=range(info[[i]]$breaks),yaxt="n")
      }
    }
    else {
      plot(NULL, type = "n", xlim = c(0, max(info[[i]]$counts)), ylim = c(range(info[[i]]$breaks)),xlab="",xaxt="n",ylab="")
      rect(0, info[[i]]$breaks[1:(length(info[[i]]$breaks) - 1)], info[[i]]$counts, info[[i]]$breaks[2:length(info[[i]]$breaks)])
      for (j in 1:Hn) {
        par(mar=c(0.5,0.25,0.5,0.25),fig=c(coord$LR[j+1],coord$LR[j+2],coord$BT[i+1],coord$BT[i]),
            new=T,las=2,mgp=c(0.5,0.25,0),tck=0,cex.axis=0.5)
        boxplot(info[[i]]$values[,Hn+1] ~ info[[i]]$values[,j],ylim=range(info[[i]]$breaks),xaxt="n",yaxt="n")
      }
    }

  }
}


#plot2.cat(value=d$DUR_EVTM, catV=factor(d$ROC), catH=d[, c("CHF", "RACE","AGE_ADMY")])




### Modifying from Chelsea's plot2.cat to make plot not by factor
plotpair1<- function(yvar, cvars, tp, dat) {
  op <- par('mar','fig','mgp','cex.axis','las','srt')
  on.exit(par(op))
 Hn <-length(cvars)
  info=dat[,c(cvars, yvar)]
value=dat[,yvar]
hv=hist(value,plot=F);
xmax=max(hv$counts)

(h <- 7/Hn+1) # 3.33

coord <- list(LR=seq(0.5/7.5,1,length=Hn+2),BT=seq(1-0.3/h,0.3/h,length=2));
#coord <- list(LR=seq(0.5/7.5,1,length=Hn+2),BT=seq(1-0.5/h,0.5/h,length=2));

  for (j in 1:Hn) {
    par(mar=c(0,0,0,0),
	fig=c(coord$LR[j+1],coord$LR[j+2],coord$BT[1],1)#;fig
	,new=(j!=1))
    plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    mtext(cvars[j],1,-1,cex=0.75)
	}
    par(mar=c(0,0,0,0),fig=c(0,coord$LR[1],coord$BT[2],coord$BT[1]),las=0,new=T)
    plot(c(0, 2), c(0, 2), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
    mtext(yvar,2,-2,cex=0.75)
    par(mar=c(0.25,0.25,0.25,0.25),fig=c(coord$LR[1],coord$LR[2],coord$BT[2],coord$BT[1]),
        new=T,las=1,mgp=c(0.5,0.25,0),tck=-0.025,cex.axis=0.65)
 hist(dat[,yvar], xlab=" ",ylab="Count", main="")
   
      for (j in 1:Hn) {
        par(mar=c(0.25,0.25,0.25,0.25),fig=c(coord$LR[j+1],coord$LR[j+2],coord$BT[2],coord$BT[1]),
            new=T,las=2,mgp=c(0.25,0.25,0.25),tck=0,cex.axis=0.65)
       
	if (tp[j]==2) {boxplot(value ~ info[,j],ylim=range(hv$breaks),yaxt="n", pch=20, cex=0.75)}else{
	plot(value ~ info[,j], pch=20, cex=0.75);abline(lm(value~info[,j]))}
      }
}

#pdf("onerow_d.pdf", width=23, height=5)
#plotpair1("DUR_EVTM", setdiff(pvars1, "DUR_EVTM"), tp=c(1,2,2,2,2,2,2,2,2,2,2,2,2,2), dat=dd)
#dev.off()
today=Sys.Date()

find.empty=function(dat, key){apply(dat[, setdiff(names(dat), key)],1, function(x){sum(!is.na(x))==0})} 

#tab1=data.frame(year=1:5, count=2:6);tab1 
#table.number=3;
#tab=tab1;
#table.title="Experiment table"
#table.description=NULL

addtable=function(tab, number=NULL, title, description=NULL, skip=3, footnote=NULL,font.size=9, ...){
addParagraph(rtf, paste(bold(paste("Table ", number,". ", title, sep="")), description));

addNewLine(rtf);
addTable(rtf,tab,col.justify='R',header.col.justify='R',font.size=font.size,row.names=FALSE,NA.string="-", ...)
if (!is.null(footnote)){addNewLine(rtf); addParagraph(rtf, footnote)};
if (skip>0) addParagraph(rtf, paste(rep("\n", skip), claapse='') ) }


detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}



breakcat=function(x, n)
{
  #  x=d$TScore.Independence;n=3
  cutp= quantile(x, probs = seq(0,1, by=1/n), na.rm=TRUE) %>% as.numeric
  cutp[which.max(cutp)]  =cutp[which.max(cutp)] +1
  cutp[which.min(cutp)]=cutp[which.min(cutp)]-1
  catx=cut(x, cutp) 
  levels(catx)=1:n
  return(catx)
}


fmtnum <- function(x) {
  if(is.na(x)){return('-')}else
  if (abs(x) > 99999) {
    return(sprintf("%.2g", x)) # round to 2 significant digits
  } else if (abs(x)>10) {
    return(round(x)) # round to 0 decimal places
  }
  else if (abs(x) >= 1) {
    return(sprintf("%.1f", x)) # round to 1 decimal place
  } else if (abs(x) >= 0.01) {
    return(sprintf("%.2f", x)) # round to 2 decimal places
  } else {
   return(sprintf("%.2g", x)) # display 2 significant digits
  }
}


mymean=function(x){
  if(sum(!is.na(x))>0){y=mean(x, na.rm=TRUE)}else {y=NA}
  return(y)
}


mysd=function(x){
  if(sum(!is.na(x))>0){y=sd(x, na.rm=TRUE)}else {y=NA}
  return(y)
}

mymax=function(x){
  if(sum(!is.na(x))>0){y=max(x, na.rm=TRUE)}
  else {y=NA}
  return(y)
}

mymin=function(x){
  if(sum(!is.na(x))>0){y=min(x, na.rm=TRUE)}else {y=NA}
  return(y)
}


mymedian=function(x){
  if(sum(!is.na(x))>0){y=median(x, na.rm=TRUE)}else {y=NA}
  return(y)}



countmiss=function(x){ paste( sum(is.na(x)), ' (', signif(mean(is.na(x))*100,2), '%)',sep='')}



findmiss=function(dat){
  M=apply(dat,2,function(x){sum(is.na(x))})
  M[M>0]
}


fisherChi = function(tabx){
  options(warn=-1)
  uap = chisq.test(tabx,correct=F)
  testtype = "Chi-square"
  #print(c("print1",testtype))
  if(sum(uap$expected<5)>0|sum(tabx)<30){
    uap=NULL
    uap = try(fisher.test(tabx),silent=TRUE)
    testtype = "Fisher"
    #print(c("print2",testtype,is.null(uap)))
    if(length(uap)==1){ #this means there was an error with fisher
      uap = chisq.test(tabx,simulate.p.value=TRUE, B = 50000)
      testtype = "Sim.Chi"
      #print(c("print3",testtype))
    }#end if
  }#end if
  options(warn=0)
  return(list(uap$p.value,testtype))
}#end fisherChi



mrange <- function(x, dg=NULL) {
  if (is.null(dg)) {
    out <- ifelse(sum(!is.na(x)) > 0, paste("(", fmtnum(mymin(x)), ", ",  fmtnum(mymax(x)), ")", sep=""), NA)
  } else {
    ins <- paste("%.", dg, "f", sep="")
    out <- ifelse(sum(!is.na(x)) > 0, paste("(", sprintf(ins, round(mymin(x), dg)), ", ", sprintf(ins, round(mymax(x), dg)), ")", sep=""), NA)
  }
  return(out)
}



miqr=function(x, dg=NULL){
if (sum(!is.na(x))>0) {nums=quantile(x, c(0.5, 0.25, 0.75), na.rm=TRUE)
  if (is.null(dg)){
return(paste(fmtnum(nums[1]), " (", fmtnum(nums[2]), ", ",  fmtnum(nums[3]), ")", sep=""))
  }else{
    ins=paste("%.",dg,"f",sep="")
    return(paste(sprintf(ins, nums[1]), " (", sprintf(ins, nums[1]) , ", " , sprintf(ins, nums[1]), ")", sep=""))
  }
} else {return(NA)}
}



msd=function(x, dg=NULL){
  if (is.null(dg)){
  out=ifelse(sum(!is.na(x))>0,paste(fmtnum(mymean(x))," (",fmtnum(mysd(x)),")",sep=""), NA)}else{
  ins=paste("%.",dg,"f",sep="")
  out=ifelse(sum(!is.na(x))>0,paste( sprintf(ins,  round(mymean(x),dg))," (",sprintf(ins, round(mysd(x),dg)),")",sep=""), NA)}
  return(out)
}


npct=function(x,dg=NULL){
  tx=table(x)
  ptx=prop.table(tx);ptx
  #convert ptx to percentages
  NP=paste(tx,   '(', signif(ptx*100,2), '%)',sep='')
  if(!is.null(dg)){
    ins=paste("%.",dg,"f",sep="")
    NP=paste(tx,   '(', sprintf(ins, ptx*100), '%)',sep='')}
  out=data.frame(levels=names(tx), NP)
  return(out)
} # disply counts as is, but percentages to 2 significant digits or specified digits

# this function computes N(%) overall and by group, and also computes row or column percents
npct.grp=function(xvar, group, pcol=TRUE, dg=NULL){
xvar=factor(xvar)
tm=table(xvar, group)
rowpct=(t(apply(tm, 1, function(x){x/sum(x)})))
colpct=apply(tm,2, function(x){x/sum(x)})
if(!is.null(dg))
{
  ins=paste("%.",dg,"f",sep="")
if (!pcol){
  tmp=matrix(paste(tm, '(', sprintf(ins, 100*rowpct), '%)', sep=''), ncol=ncol(tm))
}else{
 tmp=matrix(paste(tm, '(', sprintf(ins, 100*colpct), '%)', sep=''), ncol=ncol(tm))}

}else{
if (!pcol)
 {tmp=matrix(paste(tm, '(', signif(100*rowpct, 2), '%)', sep=''), ncol=ncol(tm))}else{
  tmp=matrix(paste(tm, '(', signif(100*colpct,2), '%)', sep=''), ncol=ncol(tm))
 }
}
tmp=as.data.frame(tmp)
names(tmp)=sort(unique(group))
taball=cbind(npct(xvar), tmp)
return(taball)
}



prfmt=function(Dtab, Var="Variable", Lev="Levels", Pval="P-value", Test="Test", Comp=TRUE){
  Dtab2=Dtab
  Dtab2[,Var][duplicated(Dtab2[,Var])]=NA
  tcode=c("c", "t", "w", "e", "f", "s","k", "j", "a", "l","q");
  tname=c("{\\super c} Chi-squared test",
          "{\\super t} T-test",
          "{\\super w} Wilcoxon rank sum test",
          "{\\super e} Exact wilcoxon rank sum test",
          "{\\super f} Fisher's exact test",
          "{\\super s} Chi-squared test by Montecarlo simulation",
          "{\\super k} Kruskal-Wallis test",
          "{\\super j} Jonckheere-Terpstra test",
          "{\\super a} ANOVA",
          "{\\super l} ANOVA test for linear trend",
          "{\\super q} Cochran-Armitage test for trend");

  if (Comp==TRUE) {
    Dtab2$SEP=bold(" - ");
    Dtab2$SEP[is.na(Dtab2[,Var]) | Dtab2[,Lev]==""]=""

 # if Lev colum is 'Missing #(%)' , make the entire row italic
it=which(Dtab2[,Lev]=='Missing #(%)')
#it=which(Dtab2[,Lev]=='Missing #(99)')
it
for (wh in it){
  for (j in 1:ncol(Dtab2)){
    Dtab2[wh,j]=ifelse(is.na(Dtab2[wh, j]), NA, italic(Dtab2[wh,j]))
  }}

Dtab2[,Var][is.na ( Dtab2[,Var] ) ]=""

    Dtab2[,Var]=paste(bold(Dtab2[,Var]), Dtab2$SEP, Dtab2[,Lev], sep=""); # collapse the variable and level into one column
    Dtab2=Dtab2[,!names(Dtab2) %in% c("SEP",Lev)] # remove the sep column
  }
  if (Pval %in% names(Dtab) & Test %in% names(Dtab)){
    P=Dtab2[,Pval];
    TT=tolower(substr(Dtab2[,Test],1,1));
    alltest=unique(na.omit(TT))
    newp=super(P, TT);
    newp[is.na(P)]=NA;
    Dtab2[, Pval]=newp
    Dtab2=Dtab2[,names(Dtab2)!=Test]
    alltest2=plyr::mapvalues(alltest, tcode, tname, warn_missing=FALSE)
    Testnote=paste(paste(alltest2, collapse=", "),".", sep="")
    return(list(Dtab2, Testnote))
  }
  else {return(Dtab2)}
}

fpvalue = function(xpval){
  if (is.na(xpval)){fpv="NA"}
  else if(xpval<0.001){
    fpv = "<0.001"}
  else if(xpval>0.055){
    fpv =sprintf("%.2f", round(xpval,2))
  }
  else{fpv =sprintf("%.3f",round(xpval,3))}#end if
  return(fpv)
}#end fpvalue

bold=function(x){
  y=paste( "{\\b ", x, "}", sep="" );
  return(y)
}


italic=function(x){
y=paste( "{\\i ", x, "}", sep="" );
return(y)}

super=function(x, s){
  y=paste( x, "{\\super ", s, "}", sep="" );
  return(y)
}

countn=function(dat){M=apply(dat,2,function(x){sum(!is.na(x))}); return(M)};
