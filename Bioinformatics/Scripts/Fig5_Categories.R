library(viridis)
library (vioplot)
library(ggfortify)
library(multcomp)
library(cluster)

Project="/home/huavan/CrossRegulation"
setwd(Project) 


## Comparison of the 3 genome assemblies analyses

#### Vectors
Strains<-c("iso1","OreR","CanS","") # for the barplot
sra<-c('iso1'='iso1.r0.SRR11846566','OreR'='OreR.r1.SRR25922470','CanS'="CanS.r1.SRR5687217")

# Open TE Classification file
Classif<-read.table("Data/RepBaseDmel.Classif.txt")
colnames(Classif)<-c("TE","Superfamily","Class")

######
###### Each table is read

for (n in c(1:3)){
#filename<-paste("Summary.2/Dmel.",sra[n],".250.0.TE.summary.txt.2.txt",sep='')

filename<-paste("ResultsDmel.",names(sra[n]),"/AllTEs/Dmel.",sra[n],".250.0.TE.summary.txt.2.txt",sep='')
tab<-read.table(filename, header=TRUE,sep='\t')

### Duplicate Status and Burst columns
tab$Status2 <- tab$Status
tab$Burst2 <- tab$Burst
### Change text used for labels
tab$Status[grepl("Old pi regulating",tab$Status)]<-"A"
tab$Status[grepl("expressed",tab$Status)]<-"B"
#tab$Status[grepl("Old pi not expressed",tab$Status)]<-"B1"
#tab$Status[grepl("Old pi expressed but not regulating",tab$Status)]<-"B2"
tab$Status[grepl("No old pi copy",tab$Status)]<-"C"
tab$Status[grepl("No pi copy",tab$Status)]<-"D"
tab$Status[grepl("No young Euc",tab$Status)]<-"E"

# simplify classification of Burst (Keep two categories)
tab$Burst[grepl("No",tab$Burst) | grepl("Old",tab$Burst)]<-"NoRecentBurst"
tab$Burst[!grepl("No",tab$Burst) & !grepl("Old",tab$Burst)]<-"RecentBurst"

# Consider only the highBurst (>15)
#tab$Burst[!grepl("High",tab$Burst) ]<-"NoRecentBurst"
#tab$Burst[grepl("High",tab$Burst) ]<-"RecentBurst"
write.table(tab[,-c(11,27,28,29)],paste("ResultsDmel/Table.",sra[n],".All.CpNb.txt",sep=""), sep='\t',row.names=FALSE,col.names=TRUE, quote=FALSE)

### Create 2D table
#tab1.2<-as.data.frame(table(tab1$Status,tab1$Burst))
assign(paste(names(sra)[n],".2",sep=''),as.data.frame(table(tab$Status,tab$Burst)))
assign(names(sra)[n],tab)
}
print(iso1.2)
print(OreR.2)
print(CanS.2)


#### Merge the data, calculate mean and sd
#### Keep dissociate category B as B1 B2
tabtot<- merge(CanS.2,OreR.2, by=c("Var1","Var2"))
tabtot<- merge(tabtot,iso1.2, by=c("Var1","Var2"))
colnames(tabtot)<-c("Status","Burst","CanS","OreR","iso1")
tabtot$Mean<- round(rowMeans(tabtot[,c(3,4,5)]),2)
tabtot$sd <- round(apply(tabtot[,c(3,4,5)], 1, sd, na.rm=TRUE),2)
# Order logically
tabtot <- tabtot[order(as.numeric(rownames(tabtot))),]
#tabtot <- tabtot[c(9,10,5,6,1,2,3,4,7,8),c(1,2,5,4,3,6,7)]
write.table(tabtot,"ResultsDmel/Table.All.CpNb.txt", sep='\t',row.names=TRUE,col.names=TRUE, quote=FALSE)

########### GLM stats

library(multcomp)
#tt <- read.table("Table.All.CpNb.txt", sep="\t", header=TRUE)
tt <- read.table("ResultsDmel/Table.All.CpNb.txt", sep="\t", header=TRUE)
#tt <- read.table("/Users/huavan/Dropbox/Siddharth/Simulicron-master/CrossRegulation/Summary.2/Table.All.CpNb.txt", sep="\t", header=TRUE)
tt.rb <- tt$Burst == "RecentBurst"

mytt <- rbind(
    data.frame(Status=tt$Status[tt.rb], Strain="iso1", RecentBurst=tt$iso1[tt.rb], NoRecentBurst=tt$iso1[!tt.rb]), 
    data.frame(Status=tt$Status[tt.rb], Strain="OreR", RecentBurst=tt$OreR[tt.rb], NoRecentBurst=tt$OreR[!tt.rb]),
    data.frame(Status=tt$Status[tt.rb], Strain="CanS", RecentBurst=tt$CanS[tt.rb], NoRecentBurst=tt$CanS[!tt.rb])
)

#mytt$Status <- factor(mytt$Status, levels=c("Reg.", "No reg.", "No old pi copy", "No pi copy", "No young Euc"))
mytt$Status <- factor(mytt$Status, levels=c("A", "B", "C", "D", "E"))

mm <- glm(cbind(RecentBurst, NoRecentBurst) ~ Status, data=mytt, family="binomial")

contrasts <- rbind(
    Reg.vs.noReg = c(1,-1, 0, 0, 0),
    Reg.vs.noOld = c(1, 0,-1, 0, 0),
    Reg.vs.noPi  = c(1, 0, 0,-1, 0), 
    noReg.vs.noOld=c(0, 1,-1, 0, 0),
    noReg.vs.noPi =c(0, 1, 0,-1, 0), 
    noOld.vs.noPi =c(0, 0, 1,-1, 0))
    
glmres<-multcomp::glht(mm, linfct = multcomp::mcp(Status = contrasts)) |>
        summary() |> print()
pval<-glmres$test$pvalues
sink("ResultsDmel/Stat.All.CpNb.txt")
    multcomp::glht(mm, linfct = multcomp::mcp(Status = contrasts)) |>
        summary() |> print()
sink(NULL)



# Calculate proportion of each category (averaged on the 3 datasets)
Prop.iso1<-aggregate(Freq~Var1, data=iso1.2,sum)
Prop.iso1$Prop<-round(Prop.iso1$Freq/sum(Prop.iso1$Freq)*100,2)
Prop.OreR<-aggregate(Freq~Var1, data=OreR.2,sum)
Prop.OreR $Prop<-round(Prop.OreR $Freq/sum(Prop.OreR $Freq)*100,2)
Prop.CanS<-aggregate(Freq~Var1, data=CanS.2,sum)
Prop.CanS $Prop<-round(Prop.CanS $Freq/sum(Prop.CanS $Freq)*100,2)

TabTot2<- Reduce(function(...) merge(..., all = TRUE, by = "Var1"), list(Prop.iso1, Prop.OreR, Prop.CanS))
TabTot2$Mean<-rowMeans(TabTot2[,c("Prop.x","Prop.y","Prop")])
#TabTot2$Label<-paste(TabTot2 $Var1,' (',round(TabTot2 $Mean,2),' %)',sep='')
TabTot2$Label<-TabTot2 $Var1
TabTot2

#### Merge the data, category of each TE for each assembly
tabFinal<- merge(CanS[,c(1,2,3,25,26)],OreR[,c(1,2,3,25,26)], all=TRUE,by=c("TE","Class","Superfamily"))
tabFinal <- merge(tabFinal,iso1[,c(1,2,3,25,26)], all=TRUE,by=c("TE","Class","Superfamily"))
colnames(tabFinal)<-c("TE","Class","Superfamily","CanS.Status","CanS.Burst","OreR.Status","OreR.Burst","iso1.Status","iso1.Burst")
write.table(tabFinal,"ResultsDmel/TableFinal.All.Category.txt", sep='\t',row.names=FALSE,col.names=TRUE, quote=FALSE)


##########################################
#### Figure5_Dmel
##########################################
# Dissociate Recent and NotThatRecent categories for the barplot
#tabtotHighRecent<-tabtot[tabtot$Burst=='RecentHighBurst',c(5,4,3)]
tabtotRecent<-tabtot[tabtot$Burst =='RecentBurst',c(5,4,3)] #+ tabtotHighRecent
tabtotNotRecentOnly<-tabtot[tabtot$Burst!='RecentBurst',c(5,4,3)]
tabtotNotRecent<-tabtotNotRecentOnly + tabtotRecent #+ tabtotHighRecent
rownames(tabtotRecent)<-tabtot$Status[tabtot$Burst =='RecentBurst']
rownames(tabtotNotRecent)<-tabtot$Status[tabtot$Burst !='RecentBurst']
vtotRecent<-tabtotRecent
vtotRecent$na <- rep("",5)
vtotRecent[vtotRecent==0]<-""

vtotNotRecent<-tabtotNotRecentOnly
vtotNotRecent$na <- rep("",5)
vtotNotRecent[vtotNotRecent ==0]<-""

Strains<-c("CanS","OreR","iso1","")
ncat<-length(unique(tabtot$Status))

# Define colors
colbase<-viridis(10*ncat)
colvecNR<-colbase[c(1,2,3,11,12,13,21,22,23,31,32,33,41,42,43)]
colvecR<-colbase[rev(c(7,8,9,17,18,19,27,28,29,37,38,39,47,48,49))]
colvecTR<-gsub('.{2}$',"7D",colvecR)

#pdf("/Users/huavan/Dropbox/Siddharth/Simulicron-master/CrossRegulation/Figure6.pdf", width=10,height=5)
pdf("ResultsDmel/Figure5.pdf", width=10,height=5)
par(mar=c(7,18,5,2))

plot(NULL, xlim=c(1,ncat*length(Strains)), ylim=c(0,100),bty='n',xaxt='n',ylab='Percentage of TE families',xlab='',cex.lab=1,las=1)
abline(h=seq(0,100,10), col="gray")

barplot(t(as.matrix(tabtotNotRecent))/t(as.matrix(tabtotNotRecent))*100, beside = TRUE, legend.txt=NULL, cex.names=1,las=0, line=2,  col= "gray50",horiz=F,add=TRUE,yaxt='n',xaxt='n')
barplot(t(as.matrix(tabtotRecent))/t(as.matrix(tabtotNotRecent))*100, beside = TRUE, yaxt='n',legend.txt=NULL,col= colvecR,horiz=F,add=TRUE,xaxt='n')
atx<-seq(-2,(4*ncat),4)+0.5
atx[1]<-atx[1]+1
adjx<-c(1,rep(0.5,5))
#barplot(t(as.matrix(tabtotHighRecent)), beside = TRUE,  yaxt='n', names.arg=rep("",length(unique(tabtot$Status))),col= colvecR,horiz=TRUE,add=TRUE,yaxt='n')
text(labels=rev(as.vector(t(vtotRecent))),y=5,x=c((4*ncat):1)+0.5,cex=0.7)
text(labels=rev(as.vector(t(vtotNotRecent))),y=97,x=c((4*ncat):1)+0.5,cex=0.7)
mtext(text=rep(Strains, ncat), side=3, at=c((4*ncat):1),adj=1,line=-0.3,las=1,cex=0.7)
#mtext(text=c("",TabTot2$Label), side=1, at=atx,adj=adjx,line=0,las=1,cex=1)
mtext(text=c(" ",as.vector(TabTot2$Label)), side=1, at=atx,adj=adjx,line=0,las=1,cex=1)
mtext(text=c("",paste('(',round(TabTot2$Mean,2),' %)'),sep=''), side=1, at=atx,adj=adjx,line=0.7,las=1,cex=0.7)
mtext(text=c("Young Euchromatic copies","YES", "YES","YES","YES","NO"), side=1, at=atx,adj=adjx,line=2,las=1,cex=0.7)
mtext(text=c("Sequences in piRNA clusters","YES", "YES","YES","NO",""), side=1, at=atx,adj=adjx,line=3,las=1,cex=0.7)
mtext(text=c("Old sequences in piRNA clusters","YES", "YES","NO","",""), side=1, at=atx,adj=adjx,line=4,las=1,cex=0.7)
mtext(text=c("Old sequences in piRNA clusters potentially able to regulate","YES", "NO","","",""), side=1, at=atx,adj=adjx,line=5,las=1,cex=0.7)
par(xpd=T)
legend("topleft", bty='n',inset=c(-0.37,0), legend=c("  No Recent Burst"), pch=c(15), pt.cex=4,col=c("gray50"),cex=0.8)
legend("left", bty='n',inset=c(-0.37,0), legend=c(rep("",3),"  Recent Burst",rep("",3)), pch=c(15), pt.cex=4,col=c(colvecR),cex=0.8)

segments(x0=c(2.5,2.5),x1=c(6.5,10.5),y1=c(115,125),y0=c(115,125), lwd=2)
mtext(text=c('***','***'), side=3,line=c(1.2,2.6),at=c(4.5,6.5))
dev.off()

# pdf("ResultsDmel/Figure5_Categories.pdf", width=8,height=6)
# par(mar=c(3,2,1,1))
# colbase[rev(c(7,8,9,17,18,19,27,28,29,37,38,39,47,48,49))]
# colvecTot<-c(colbase[rev(c(17,27,37,39,47))],"gray90","gray95")

# tab4plot<-t(as.matrix(tabFinal[,c("iso1.Status","CanS.Status","OreR.Status")]))
# tab4plot[is.na(tab4plot)]<-"F"
# tab4plot2<-rbind(table(tab4plot[1,]),table(tab4plot[2,]))
# tab4plot2<-rbind(tab4plot2,table(tab4plot[3,]))
# row.names(tab4plot2)<-c("iso1","CanS","OreR")
# tab4plot2<-as.matrix(t(tab4plot2))

# leg<-c("Regulation", "Expressed but not mapping", "Not expressed", "No old piRNA Cluster copy", "No piRNA Cluster copy","No young euchromatic copy","No euchromatic copy")
# plot(NULL, xlim=c(0,6), ylim=c(0,200),bty='n',xaxt='n',ylab='Number of TE families',xlab='',cex.lab=1,las=1)

# barplot(c(179,179,179), beside=F, cex.names=1,las=0,  col= c("gray95"),horiz=F,add=T,yaxt='n',xaxt='n')
# barplot(tab4plot2, beside=F, names.arg=c("iso1","CanS","OreR"),cex.names=1,las=0,  col= colvecTot,horiz=F,add=T,yaxt='n')
# #barplot(table(tab4plot[2,]), beside=F,legend.txt=NULL, cex.names=1,las=0,   col= colvecTot,horiz=F,add=TRUE,yaxt='n',xaxt='n')
# legend("right", bty='n', legend=rev(leg), pch=c(15), pt.cex=2,col= rev(colvecTot),cex=1)

# dev.off()


# # #pdf("/Users/huavan/Dropbox/Siddharth/Simulicron-master/CrossRegulation/Figure6-sup2.pdf", width=16,height=20)

# pdf("ResultsDmel/Figure5-sup2.pdf", width=16,height=20)
# par(mfrow=c(5,2),mar=c(6,12,3,2))
# colvector<-viridis(5)
# colvector2=c("grey",colvector[1],"grey",colvector[2],"grey",colvector[3],"grey20","grey",colvector[5])

# tab<-iso1
# a<-table(tab$Status2)
# anames<-paste(names(a)," (", a,")",sep='')
# ds<-c("All TEs","RecentBurst", "NoRecentBurst")
# vioplot(tab$RootDepth ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'RootDepth',sep=' '),drop=TRUE)
# vioplot(tab$mrcaDepth ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'mrcaDepth',sep=' '),drop=TRUE)
# vioplot(tab$NbEucYoung ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'NbEucYoung',sep=' '),drop=TRUE)
# vioplot(tab$NbPiYoung ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'NbPiYoung',sep=' '),drop=TRUE)
# vioplot(tab$NbEucOther ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'NbEucOther',sep=' '),drop=TRUE)
# vioplot(tab$NbPiOther ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'NbPiOther',sep=' '),drop=TRUE)
# vioplot(tab$NbEucOld ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'NbEucOld',sep=' '),drop=TRUE)
# vioplot(tab$NbPiOld ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'NbPiOld',sep=' '),drop=TRUE)
# vioplot(tab$TotalNbReads ~ tab$Burst+tab$Class,las=2,xlab='',col=colvector2, main=paste(ds[1],'TotalNbReads',sep=' '),drop=TRUE)
# vioplot(tab$FreqMinus ~ tab$Burst+tab$Class,las=2,xlab='',col= colvector2, main=paste(ds[1],'FreqMinus',sep=' '),drop=TRUE)

# dev.off()




# # 

 # #pdf("/Users/huavan/Dropbox/Siddharth/Simulicron-master/CrossRegulation/Figure6-PCA.pdf",width=12,height=8)
 # pdf("ResultsDmel/Figure5-PCA.pdf",width=12,height=8)
 # Assembly <-list(CanS,OreR,iso1)
 # for (n in 1:length(Assembly)){
	
 # assembly<-Assembly[[n]]
 # pchvector=c('NoRecentBurst'=1,'RecentBurst'=2)
 # pchvector2=c('DNA'=1,'LINE'=2,'LTR'=0,'RC'=3)
 # tabpca<-assembly[,c( "TotalNbReads", "FreqMinus","RootDepth", "mrcaDepth", "NbEucOld" ,"NbEucOther" ,"NbPiOld", "NbPiOther","NbPiYoung","NbEucYoung", "MaxYEOP")]
 # # Avoid 0 before log transformation 
 # #1
 # tabpca$TotalNbReads <- log(tabpca$TotalNbReads) #1
 # #2
 # tabpca$FreqMinus <- tabpca$FreqMinus #2
 # #3
 # tabpca$RootDepth <- tabpca$RootDepth #3
 # #4
 # tabpca$mrcaDepth <- tabpca$mrcaDepth #4
 # #5
 # tabpca$NbEucOld <- log(tabpca$NbEucOld +1) #5
 # #6
 # tabpca$NbEucOther <- log(tabpca$NbEucOther+1) #6
 # #7
 # tabpca$NbPiOld <- log(tabpca$NbPiOld +1) # 1: common status in two dataset
 # #8
 # tabpca$NbPiOther <- log(tabpca$NbPiOther +1) #8
 # #9
 # tabpca$NbPiYoung <- log(tabpca$NbPiYoung +1) #9
 # #10
 # tabpca$NbEucYoung <- log(tabpca$NbEucYoung +1) #10
 # #11
 # tabpca$MaxYEOP <- log(tabpca$MaxYEOP +1) #11
 # #12
 # tabpca$TotEuc <- tabpca$NbEucYoung+tabpca$NbEucOther+tabpca$NbEucOld # 3: same status for all three dataset
 # #13
 # tabpca$PrEucYoung <- tabpca$NbEucYoung/tabpca$TotEuc # 3: same status for all three dataset
 # #14
 # tabpca$TotPi <- tabpca$NbPiYoung+tabpca$NbPiOther+tabpca$NbPiOld # 3: same status for all three dataset
 # #15
 # tabpca$PrPiOld <- tabpca$NbPiYoung/tabpca$TotPi # 3: same status for all three dataset
 # #16
 # tabpca$PrPiEucYoung <- tabpca$NbPiYoung/(tabpca$NbPiYoung+tabpca$NbEucYoung +1)# 3: same status for all three dataset

 # tab4pca<-tabpca[,c(1:4,11,16)]
 # res<-prcomp(tab4pca, scale.=TRUE)
 # autoplot(res, data= assembly, colour="Status",loadings = TRUE, loadings.label = TRUE,loadings.label.size =5, shape="Burst",size=4, frame=TRUE,frame.colour="Status")+scale_colour_viridis(
  # alpha = 1,
  # begin = 0,
  # end = 1,
  # direction = 1,
  # discrete = TRUE,
  # option = "D",
  # aesthetics = "color"
 # ) + 
 # geom_point(shape = pchvector[assembly $Burst],size = 4,colour = "black")+

  # theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", linewidth=1), axis.text=element_text(size=16),axis.title=element_text(size=16), legend.text=element_text(size=14),legend.title=element_text(size=16),plot.margin=margin(5,5,2,2))
 # #ggsave(paste("/Users/huavan/Dropbox/Siddharth/Simulicron-master/CrossRegulation/Figure6-PCA",Strains[n],".pdf",sep=''),width=20, height=12)
 # ggsave(paste("ResultsDmel/Figure6-PCA",Strains[n],".pdf",sep=''),width=20, height=12)
 # #biplot(res)
 # #plsda(assembly,assembly$status)
 # }

 # dev.off()


# # #table(iso1$Status, iso1 $Burst)
# # #plot(iso1 $TotalNbReads, iso1 $FreqMinus,col= colvector[iso1 $Status],pch=pchvector[iso1 $Burst],bg="grey")
# # #table(iso1 $Status, iso1 $Class)
# # #fisher.test(table(iso1 $Status, iso1 $Class))

# # #table(iso1 $Burst, iso1 $Class)


