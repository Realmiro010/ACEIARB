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

# connect-S plot
nv<-length(colnames(ASD1))
ncov<-length(rownames(ASD1))

mydata <- data.frame(Subgroups = rep(colnames(ASD1), each=nrow(ASD1)),
                     VarName = rep(rownames(ASD1), nv),
                     ASMD1 = apply(ASMD1,1,mean),
                     ASMD2 = apply(ASMD2,1,mean),
                     ASMDb = apply(ASMDb,1,mean))

mydata$VarName <- factor(mydata$VarName, levels=unique(mydata$VarName))
mydata$Subgroups=mapvalues(mydata$Subgroups,unique(mydata$Subgroups), c("Overall","Age<75", "Age>75", "Male","Female","Black", "Not Black",  "Standard", "Intensive","No MCI","MCI")) 
mydata$Subgroups <- factor(mydata$Subgroups, levels=unique(mydata$Subgroups))
head(mydata)


## only doing this to get nicer labels
unique(mydata$VarName)
tp1=tp[covars %in% setdiff(predvars, subgroup)]
tp1=tp[covars %in% predvars]
varlabs1=varlabs[covars %in% predvars]
length(varlabs1)
mydata$VarName

tt=myNOgrtab(xdat = dm1[,predvars],  typev=tp1, dnames=varlabs1, allSums = FALSE, annot = FALSE,  onelinebinary = TRUE, showbinarylevel = FALSE)
tlabs=tt[[1]][, c("Variable", "Levels")]
labs=tlabs$Levels
labs[labs==""]=tlabs$Variable[labs==""]
reflev=c("NH White", "Current", "Less than high school diploma","Doctors office/outpatient clinic")
labs2=setdiff(labs, reflev)
labs2=mapvalues(labs2, c("Former", "Never"), c("Former smoker","Never smoked"))
cbind(unique(as.character(mydata$VarName)), labs2)


toplot<- mutate(mydata, ASMD=ASMD2, lcolor = ifelse(ASMD <=0.1 |is.na(ASMD), "0",
                                                    ifelse((ASMD <=0.15), "1",
                                                           ifelse((ASMD <=0.20), "2","3"))),
                VarName=rep(labs2, nv) %>% fct_inorder())

toplot$lcolor <- factor(toplot$lcolor, levels=c("0","1","2","3"))
group.colors <- c("0" = "White", "1" = "grey85", "2" ="grey50","3"="black")
levels(toplot$VarName)
toplot$VarName1=reord(toplot$VarName, c(45-c(0:44)))
head(toplot)

g1 <- ggplot(toplot, aes(x =VarName1 , y =Subgroups , fill=lcolor)) +
  geom_dotplot(binaxis = "y",binwidth = 0.2, stackdir = "center")+labs(x="Covariate", y="Subgroup"  )+
  scale_fill_manual(values=group.colors,name="ASMD",  labels=c("<0.1", "0.1-0.15", "0.15-0.20", ">0.20"), drop=FALSE)+coord_flip()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text = element_text( size = 12 ),
        axis.title = element_text( size = 16, face = "bold" ),legend.position="bottom",legend.title = element_text(size=18),
        legend.text=element_text(size=16), legend.key.size = unit(2,"line"))+
  annotate("text", x = rep(ncov+1,nv), y = seq(1, nv,1), label = format(nsubg,zero.print = T),size = 5)+
  
  #annotate("text", rep(ncov+2,nv), y = seq(1, nv,1), label = format(VIFW,zero.print = T),size = 5)+
  #  annotate("text", x = nv+1, y = nv+1, label ="N",size = 5)+
  expand_limits(x= c(0, length(levels(toplot$VarName)) + 3),y= c(0, length(levels(toplot$Subgroups)) + 1.5) )

(splot <- g1+ theme(panel.border = element_rect(linetype = "solid", fill = NA)))
ggsave("Connect_S_OW24.pdf", plot=splot, width=14, height=14)

