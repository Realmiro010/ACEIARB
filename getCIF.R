
##&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
###           CIF plots                       #************************************************************************
##    &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
mytheme=theme(panel.background = element_blank(), plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.line = element_line(linetype = 1), 
              legend.key=element_rect(fill = "white", colour = "white"),  
              legend.key.size=unit(2,"lines"), legend.text=element_text(size=14) ,legend.box="berticle", plot.title = element_text(size=16,face="bold", margin=margin(0,0,15,0))  , 
              axis.text=element_text(size=14), axis.title=element_text(size=14))
## Death ###
dm1$fwt=1#- unweighted
subgroup=c("Age75", "BLFemale", "Black","IntensiveTrt", "BLMCI")

cif.death.raw=wtkm(dat=dm1,evtvar="AnyDeathEvent", dayvar="death_days", wt="fwt")

cif.death=NULL # weighted/pooled - ATE
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifi=wtkm(dat=dat,evtvar="AnyDeathEvent", dayvar="death_days", wt="ATE")
cif.death =cbind(cif.death, cifi$cif)}
cif.death.ate=cif.death.raw %>% mutate(cif=apply(cif.death,1,mean))

cif.death=NULL# weighted/pooled
for (i in 1:10){dat = imp.long[imp.long$.imp==i, ]
cifi=wtkm(dat=dat,evtvar="AnyDeathEvent", dayvar="death_days", wt="OW")
cif.death =cbind(cif.death, cifi$cif)}
cif.death.ow=cif.death.raw %>% mutate(cif=apply(cif.death,1,mean))
tablem(dm1$evtcomp)
tablem(dm1$pd_mci_protocolcomp)

### All primary/secondarey outcomes 
###  - pd_mci_amnestic- unweighted
raw.cif1=getcif(dat=dm1, tx=1, evtvar="evtcomp", dayvar="comp_days")%>% as.data.frame %>% mutate(Year=tm1/365)
raw.cif0=getcif(dat=dm1, tx=0, evtvar="evtcomp", dayvar="comp_days")%>% as.data.frame %>% mutate(Year=tm1/365)
cif.e1.raw=rbind(raw.cif1 %>% mutate(group="ARB"), raw.cif0 %>% mutate(group="ACEI"))%>%select(Year, cif, group)

# Combine with death then plot
ddplot=rbind(cif.e1.raw %>% mutate(type="Probable Dementia or amnestic MCI"), cif.death.raw %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)

dlab=Table3

myfunc=function(x){gsub( ","," - ", x)}
myfunc.vector=function(xv){sapply(xv, myfunc)}

dlab=dlab %>% mutate_at(vars(starts_with("HR")), myfunc.vector)

leg1=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Probable Dementia or amnestic MCI (censoring death)"])
#leg2=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Death"])
head(ddplot)

#dev.off()  

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group)) +scale_x_continuous(limits = c(0,6.3), breaks=0:6)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Probable Dementia or Amnestic MCI", hjust=0, size=5)
(Cifraw=p3+ggtitle("Before weighting")+mytheme)

##Event wighted
cif.e1.ate=cumcif("evtcomp", "comp_days")
cif.e1.ow=cumcif("evtcomp", "comp_days", "OW")

ddplot=rbind(cif.e1.ate %>% mutate(type="Probable Dementia or amnestic MCI"), cif.death.ate %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Probable Dementia or amnestic MCI (censoring death)"])
#leg2=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Probable Dementia or Amnestic MCI", hjust=0, size=5)
(Cifweighted1=p3+ggtitle("After weighting (ATE)")+mytheme)

ddplot=rbind(cif.e1.ow %>% mutate(type="Probable Dementia or amnestic MCI"), cif.death.ow %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Probable Dementia or amnestic MCI (censoring death)"])
#leg2=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Probable Dementia or Amnestic MCI", hjust=0, size=5)
(Cifweighted2=p3+ggtitle("After weighting (OW)")+mytheme)

ggarrange(Cifraw, Cifweighted1, Cifweighted2,nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
#ggsave("CIF_Primary2022_24mos.pdf", width = 19, height=7, dpi=300)
#ggsave("CIF_Primary2022_24mos.png", width = 19, height=7, dpi=300)


## Secondary outcomes 
## PD MCI only
raw.cif1=getcif(dat=dm1, tx=1, evtvar="mci_protocolcomp", dayvar="mci_protocol_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
raw.cif0=getcif(dat=dm1, tx=0, evtvar="mci_protocolcomp", dayvar="mci_protocol_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
cif.e2.raw=rbind(raw.cif1 %>% mutate(group="ARB"), raw.cif0 %>% mutate(group="ACEI"))%>%select(Year, cif, group)

# Combine with death then plot
ddplot=rbind(cif.e2.raw %>% mutate(type="Protocol-defined MCI alone"), cif.death.raw %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)

leg1=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Protocol-defined MCI alone"])
#leg2=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22,    label="Protocol-defined MCI alone", hjust=0, size=5)
(Cifraw=p3+mytheme)

cif.e2.ate=cumcif("mci_protocolcomp", "mci_protocol_death_days")
cif.e2.ow=cumcif("mci_protocolcomp", "mci_protocol_death_days" ,"OW")

ddplot=rbind(cif.e2.ate %>% mutate(type="Protocol-defined MCI alone"), cif.death.ate %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Protocol-defined MCI alone"])
#leg2=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Protocol-defined MCI alone", hjust=0, size=5)
(Cifweighted1=p3+mytheme)

ddplot=rbind(cif.e2.ow %>% mutate(type="Protocol-defined MCI alone"), cif.death.ow %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Protocol-defined MCI alone"])
#leg2=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, 
                                                                                                                                                     label="Protocol-defined MCI alone", hjust=0, size=5)
(Cifweighted2=p3+ggtitle(" After IP  Weight (OW)")+mytheme)

Cifraw.pdmci=Cifraw
Cifwt1.pdmci=Cifweighted1
Cifwt2.pdmci=Cifweighted2
ggarrange(Cifraw, Cifweighted1, Cifweighted2,nrow=1, ncol=3, common.legend=TRUE, legend="bottom")

#ggsave("CIF_MCIonly2022_24mos.pdf", width = 19, height=7, dpi=300)
#ggsave("CIF_MCIonly2022_24mos.png", width = 19, height=7, dpi=300)

### Amnestic MCI alone
raw.cif1=getcif(dat=dm1, tx=1, evtvar="mci_amnesticcomp", dayvar="mci_amnestic_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
raw.cif0=getcif(dat=dm1, tx=0, evtvar="mci_amnesticcomp", dayvar="mci_amnestic_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
cif.e3.raw=rbind(raw.cif1 %>% mutate(group="ARB"), raw.cif0 %>% mutate(group="ACEI"))%>%select(Year, cif, group)

# Combine with death then plot
ddplot=rbind(cif.e3.raw %>% mutate(type="Amnestic MCI alone"), cif.death.raw %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)

leg1=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Amnestic MCI alone"])
#leg2=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, 
                                                                                                                                                     label="Amnestic MCI alone", hjust=0, size=5)
(Cifraw=p3+ggtitle(" ")+mytheme)

cif.e3.ate=cumcif("mci_amnesticcomp", "mci_amnestic_death_days")
cif.e3.ow=cumcif("mci_amnesticcomp", "mci_amnestic_death_days" ,"OW")

ddplot=rbind(cif.e3.ate %>% mutate(type="Amnestic MCI alone"), cif.death.ate %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Amnestic MCI alone"])
#leg2=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Amnestic MCI alone", hjust=0, size=5)
(Cifweighted1=p3+ggtitle(" ")+mytheme)

ddplot=rbind(cif.e3.ow %>% mutate(type="Amnestic MCI alone"), cif.death.ow %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Amnestic MCI alone"])
#leg2=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, 
                                                                                                                                                     label="Amnestic MCI alone", hjust=0, size=5)
(Cifweighted2=p3+ggtitle(" ")+mytheme)
Cifraw.mci=Cifraw
Cifwt1.mci=Cifweighted1
Cifwt2.mci=Cifweighted2
ggarrange(Cifraw, Cifweighted1, Cifweighted2,nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
#ggsave("CIF_Amnestic_MCIonly2022_24mos.pdf", width = 19, height=7, dpi=300)
#ggsave("CIF_Amnestic_MCIonly2022_24mos.png", width = 19, height=7, dpi=300)

### PD alone
raw.cif1=getcif(dat=dm1, tx=1, evtvar="pdcomp", dayvar="pd_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
raw.cif0=getcif(dat=dm1, tx=0, evtvar="pdcomp", dayvar="pd_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
cif.e4.raw=rbind(raw.cif1 %>% mutate(group="ARB"), raw.cif0 %>% mutate(group="ACEI"))%>%select(Year, cif, group)

# Combine with death then plot
ddplot=rbind(cif.e4.raw %>% mutate(type="Probable Dementia alone"), cif.death.raw %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)

leg1=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Probable Dementia alone"])
#leg2=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Probable Dementia alone", hjust=0, size=5)
(Cifraw=p3+ggtitle(" ")+mytheme)

cif.e4.ate=cumcif("pdcomp", "pd_death_days")
cif.e4.ow=cumcif("pdcomp", "pd_death_days" ,"OW")

ddplot=rbind(cif.e4.ate %>% mutate(type="Probable Dementia alone"), cif.death.ate %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Probable Dementia alone"])
#leg2=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Probable Dementia alone", hjust=0, size=5)
(Cifweighted1=p3+ggtitle(" ")+mytheme)

ddplot=rbind(cif.e4.ow %>% mutate(type="Probable Dementia alone"), cif.death.ow %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Probable Dementia alone"])
#leg2=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, 
                                                                                                                                                     label="Probable Dementia alone", hjust=0, size=5)
(Cifweighted2=p3+ggtitle(" ")+mytheme)
Cifraw.pd=Cifraw
Cifwt1.pd=Cifweighted1
Cifwt2.pd=Cifweighted2
ggarrange(Cifraw, Cifweighted1, Cifweighted2,nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
#ggsave("CIF_PD_only2022_24mos.pdf", width = 19, height=7, dpi=300)
#ggsave("CIF_PD_only2022.png_24mos", width = 19, height=7, dpi=300)

## PD and MCI
#"pd_mci_protocolcomp", "pd_mci_protocol_death_days" "Probable Dementia or protocol-defined MCI"
raw.cif1=getcif(dat=dm1, tx=1, evtvar="pd_mci_protocolcomp", dayvar="pd_mci_protocol_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
raw.cif0=getcif(dat=dm1, tx=0, evtvar="pd_mci_protocolcomp", dayvar="pd_mci_protocol_death_days")%>% as.data.frame %>% mutate(Year=tm1/365)
cif.e5.raw=rbind(raw.cif1 %>% mutate(group="ARB"), raw.cif0 %>% mutate(group="ACEI"))%>%select(Year, cif, group)

# Combine with death then plot
ddplot=rbind(cif.e5.raw %>% mutate(type="Probable Dementia or protocol-defined MCI"), cif.death.raw %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)

leg1=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Probable Dementia or protocol-defined MCI"])
#leg2=paste("HR(95% CI):", dlab$HR.un[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, 
                                                                                                                                                     label="Probable Dementia or protocol-defined MCI", hjust=0, size=5)
(Cifraw=p3+ggtitle(" ")+mytheme)

cif.e5.ate=cumcif("pd_mci_protocolcomp", "pd_mci_protocol_death_days")
cif.e5.ow=cumcif("pd_mci_protocolcomp", "pd_mci_protocol_death_days" ,"OW")

ddplot=rbind(cif.e5.ate %>% mutate(type="Probable Dementia or protocol-defined MCI"), cif.death.ate %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Probable Dementia or protocol-defined MCI"])
#leg2=paste("HR(95% CI):", dlab$HR.ate[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, label="Probable Dementia or protocol-defined MCI", hjust=0, size=5)
(Cifweighted1=p3+ggtitle(" ")+mytheme)

ddplot=rbind(cif.e5.ow %>% mutate(type="Probable Dementia or protocol-defined MCI"), cif.death.ow %>% mutate (type="Death"))
ddplot$type=fct_inorder(ddplot$type)
tablem(ddplot$type)
leg1=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Probable Dementia or protocol-defined MCI"])
#leg2=paste("HR(95% CI):", dlab$HR.ow[dlab$Variables=="Death"])

p1=ggplot(ddplot, aes(Year, cif))+geom_step(aes(alpha=type,linetype=group))+scale_x_continuous(limits = c(0,6.3), breaks=0:7)+scale_linetype_manual(values=c(2,1), name="")
p2=p1+scale_y_continuous(limits=c(0,0.25), labels = percent_format(accuracy=1), name="Cumulative incidence") +scale_alpha_manual(values=c(1, 0.35), name=" ", guide=FALSE)
p3=p2+annotate("text",x=0, y=0.19, label=leg1, size=5, hjust=0)+annotate("text",x=5, y=0.01, label="All cause death", size=5, color="gray")+annotate("text",x=0, y=0.22, 
                                                                                                                                                     label="Probable Dementia or protocol-defined MCI", hjust=0, size=5)
(Cifweighted2=p3+ggtitle(" ")+mytheme)

Cifraw.pdpdmci=Cifraw
Cifwt1.pdpdmci=Cifweighted1
Cifwt2.pdpdmci=Cifweighted2
#ggarrange(Cifraw, Cifweighted1, Cifweighted2,nrow=1, ncol=3, common.legend=TRUE, legend="bottom")
#ggsave("CIF_PD_MCI_2022_24mos.pdf", width = 19, height=7, dpi=300)
#ggsave("CIF_PD_MCI_2022_24mos.png", width = 19, height=7, dpi=300)

ggarrange(Cifraw.pdpdmci+ggtitle("Before IP Weighting"), Cifwt1.pdpdmci+ggtitle("After IP Weighting (ATE)"), Cifraw.pd, Cifwt1.pd, Cifraw.mci, Cifwt1.mci, Cifraw.pdmci, Cifwt1.pdmci, nrow=4, ncol=2, common.legend=TRUE,legend='right')
ggsave("CIF_2nd_outcomes_12mos_efigure9.pdf", width = 13, height=18, dpi=300)

ggarrange(Cifraw.pdpdmci+ggtitle("Before IP Weighting"), Cifwt2.pdpdmci+ggtitle("After IP Weighting (OW)"), Cifraw.pd, Cifwt2.pd, Cifraw.mci, Cifwt2.mci, Cifraw.pdmci, Cifwt2.pdmci, nrow=4, ncol=2, common.legend=TRUE,legend='right')
ggsave("CIF_2nd_outcomes_12mosOW_efigure11.pdf", width = 13, height=18, dpi=300)
