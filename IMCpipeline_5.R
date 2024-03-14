####PART 5 SIMPLIFIED NETWORK VISUALIZATION

#Won Jin Ho 
#R version 4.2.2 Patched (2022-11-10 r83330)
#Platform: x86_64-pc-linux-gnu (64-bit)
#Running under: Ubuntu 18.04.5 LTS

rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

####REQUIRED LIBRARIES####
library(readxl)
library(stringr)
library(matrixcalc)
library(matrixStats)
library(Hmisc)
library(reshape2)
library(dplyr)
library(plotrix)
library(multcomp)
library(flowCore)
library(sf)
library(clusterSim)
library(limma)
library(corrplot)
library(packcircles)
library(cowplot)
library(autoimage)
library(ggplot2)
library(ggpubr)
library(ggiraphExtra)
library(ggvoronoi)
library(ggridges)
library(gridExtra)
library(igraph)
library(qgraph)
library(circlize)
library(scales)
library(RColorBrewer)
library(ComplexHeatmap)
library(pheatmap)
library(pals)
library(plot3D)
library(akima)
library(basetheme)

####LOAD DATA####
output<-readRDS('backup_output.rds')
clusterMergeFile = paste0(work,"/Config/merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

##Levels

samplevels=c(paste0("P_",1:342))
samplestoexclude <- c(output$meta_data$sample_id[output$meta_data$Tumor=="CTRL"]) # CTRL 
typelevels=c("CTRL","LUM","TNC","HER2")
sitelevels=c("PBC","PBC_1YR","PBC_tx","PBC_recur","GI","PANC","LIVER","BRAIN","OVARY","PLEURA","SPINE","LUNG","LN","CTRL")
tumorlevels=c("IDC","ILC","IMC","CTRL")
TMAlevels=c("TMA788","TMA789","TMA801","TMA974","TMA975")
SPClevels=c("SPC01","SPC02","SPC03","SPC04","SPC06","SPC07","SPC09","SPC10",
            "SPC11","SPC12","SPC13","SPC14","SPC15","SPC16","SPC17","SPC18",
            "SPC19","SPC20","SPC21","SPC22","SPC23","SPC24","SPC25","SPC26",
            "CTRL")
clusterlevels=c("Tc", "ThN", "ThEM", "Treg",
                "NK_T", "NK", "B_I", "B_II",
                "Gran", "DC",
                "Mac_I", "Mac_II", "Mac_III", "Mac_IV",
                "Str_PS", "Str_V", "Str_VC", "Str_VS",
                "CK8", "CK8_E", "CK8_V",
                "CK14", "CK14_E", "CK14_EV",
                "CK8_14_E", "CK8_14_EV",
                "CKlo", "CKlo_EV", "CKother",
                "Neuron", "NA")
clusterlevels_macs=c("CK8", "CK8_E", "CK8_V",
                     "CK14", "CK14_E", "CK14_EV",
                     "CK8_14_E", "CK8_14_EV",
                     "CKlo", "CKlo_EV", "CKother", 
                     "Neuron", "NA",
                     "Str_PS", "Str_V", "Str_VC", "Str_VS",
                     "NK_T", "NK", "B_I", "B_II",                     
                     "Tc", "ThN", "ThEM", "Treg", 
                     "Gran", "DC",
                     "Mac_I", "Mac_II", "Mac_III", "Mac_IV")

##colors

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(cluster_merging$new_cluster)))
hex <- hue_pal()(9)
colorassignedbroad <- c(rep(hex[1],3), #ck8
                        rep(hex[2],3), #ck14
                        rep(hex[3],2), #ck8and14
                        rep(hex[4],3), #ckother
                        rep(hex[5],2), #other
                        rep(hex[6],4), #str
                        rep(hex[7],8), #lym
                        rep(hex[8],2), #myeloid1
                        rep(hex[9],4)) #myeloid2
clusternames<-clusterlevels
names(colorassigned)<-clusternames

####NETWORK VISUALIZATION####
###by a distance object (of matrix) of shortest distance between cell types

# load backup data 
expr_LUMrm<- readRDS('expr_LUMrm.rds')
expr_LUMpri_mat<- readRDS('expr_LUMpri_mat.rds')
expr_LUMmet_mat<- readRDS('expr_LUMmet_mat.rds')

expr_TNCrm<- readRDS('expr_TNCrm.rds')
expr_TNCpri_mat<- readRDS('expr_TNCpri_mat.rds')
expr_TNCmet_mat<- readRDS('expr_TNCmet_mat.rds')

pct_LUM<- readRDS('pct_LUM.rds')
pct_TNC<- readRDS('pct_TNC.rds')

#for LUM
mat_LUM=aggregate(x=expr_LUMrm, by=list(rownames(expr_LUMrm)), FUN=mean, na.rm=T)
groupnames<-mat_LUM$Group.1
mat_LUM<-as.matrix(mat_LUM[,2:ncol((mat_LUM))])
rownames(mat_LUM)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_LUM)<-str_sub(colnames(mat_LUM),6,) #simplify colnames
mat_LUMex<-mat_LUM[rownames(mat_LUM)[c(rownames(mat_LUM) %nin% "30")],colnames(mat_LUM)[c(colnames(mat_LUM) %nin% "30")]] #excluding Other subtypes
dist_LUM<-mat_LUMex

#for LUM pri
mat_LUM_pri=aggregate(x=expr_LUMpri_mat, by=list(rownames(expr_LUMpri_mat)), FUN=mean, na.rm=T)
groupnames<-mat_LUM_pri$Group.1
mat_LUM_pri<-as.matrix(mat_LUM_pri[,2:ncol((mat_LUM_pri))])
rownames(mat_LUM_pri)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_LUM_pri)<-str_sub(colnames(mat_LUM_pri),6,) #simplify colnames
mat_LUM_priex<-mat_LUM_pri[rownames(mat_LUM_pri)[c(rownames(mat_LUM_pri) %nin% "30")], colnames(mat_LUM_pri)[c(colnames(mat_LUM_pri) %nin% "30")]] #excluding Other subtypes
dist_LUM_pri<-mat_LUM_priex

#for LUM met
mat_LUM_met=aggregate(x=expr_LUMmet_mat, by=list(rownames(expr_LUMmet_mat)), FUN=mean, na.rm=T)
groupnames<-mat_LUM_met$Group.1
mat_LUM_met<-as.matrix(mat_LUM_met[,2:ncol((mat_LUM_met))])
rownames(mat_LUM_met)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_LUM_met)<-str_sub(colnames(mat_LUM_met),6,) #simplify colnames
mat_LUM_metex<-mat_LUM_met[rownames(mat_LUM_met)[c(rownames(mat_LUM_met) %nin% "30")], colnames(mat_LUM_met)[c(colnames(mat_LUM_met) %nin% "30")]] #excluding Other subtypes
dist_LUM_met<-mat_LUM_metex

#for TNC
mat_TNC=aggregate(x=expr_TNCrm, by=list(rownames(expr_TNCrm)), FUN=mean, na.rm=T)
groupnames<-mat_TNC$Group.1
mat_TNC<-as.matrix(mat_TNC[,2:ncol((mat_TNC))])
rownames(mat_TNC)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_TNC)<-str_sub(colnames(mat_TNC),6,) #simplify colnames
mat_TNCex<-mat_TNC[rownames(mat_TNC)[c(rownames(mat_TNC) %nin% "30")], colnames(mat_TNC)[c(colnames(mat_TNC) %nin% "30")]] #excluding Other subtypes
dist_TNC<-mat_TNCex

#for TNC pri
mat_TNC_pri=aggregate(x=expr_TNCpri_mat, by=list(rownames(expr_TNCpri_mat)), FUN=mean, na.rm=T)
groupnames<-mat_TNC_pri$Group.1
mat_TNC_pri<-as.matrix(mat_TNC_pri[,2:ncol((mat_TNC_pri))])
rownames(mat_TNC_pri)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_TNC_pri)<-str_sub(colnames(mat_TNC_pri),6,) #simplify colnames
mat_TNC_priex<-mat_TNC_pri[rownames(mat_TNC_pri)[c(rownames(mat_TNC_pri) %nin% "30")], colnames(mat_TNC_pri)[c(colnames(mat_TNC_pri) %nin% "30")]] #excluding Other subtypes
dist_TNC_pri<-mat_TNC_priex

#for TNC met
mat_TNC_met=aggregate(x=expr_TNCmet_mat, by=list(rownames(expr_TNCmet_mat)), FUN=mean, na.rm=T)
groupnames<-mat_TNC_met$Group.1
mat_TNC_met<-as.matrix(mat_TNC_met[,2:ncol((mat_TNC_met))])
rownames(mat_TNC_met)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_TNC_met)<-str_sub(colnames(mat_TNC_met),6,) #simplify colnames
mat_TNC_metex<-mat_TNC_met[rownames(mat_TNC_met)[c(rownames(mat_TNC_met) %nin% "30")], colnames(mat_TNC_met)[c(colnames(mat_TNC_met) %nin% "30")]] #excluding Other subtypes
dist_TNC_met<-mat_TNC_metex

#color of legends revised
allcelltypes<-clusterlevels
legendctype<-as.data.frame(cbind(paste0("ctype",1:length(allcelltypes)),allcelltypes))
legendctype$maintype<-1
legendctype$maintype[str_detect(legendctype$allcelltypes,"T")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"B")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"NK")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"DC")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Mac")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Gran")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Str")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"CK")]<-"CK"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Neuron")]<-"Other"
legendctype$maintype[str_detect(legendctype$allcelltypes,"NA")]<-"Other"


markerlist_ck = c("CK8","CK14","ECAD","VIM","KI67","HLADR")

exprtbl <- data.frame(fsApply(output$fcs,exprs)[, markerlist_ck], sample_id = output$sample_ids, cluster = output$cell_clustering1m)

exprtbl$type <- factor(output$meta_data$Phenotype[match(exprtbl$sample_id,output$meta_data$sample_id)], levels=typelevels)

exprtbl_lum <- exprtbl[exprtbl$type=="LUM",]
exprtbl_tnc <- exprtbl[exprtbl$type=="TNC",]

exprtbl <- exprtbl %>% group_by(cluster) %>% summarise_at(markerlist_ck, mean, na.rm=TRUE)

exprmtx <-as.matrix(exprtbl[,2:(length(markerlist_ck)+1)])
rownames(exprmtx) <-unlist(exprtbl[,1])

exprdf<-as.data.frame(scale(exprmtx)[legendctype$allcelltypes,]) #create scaled df

exprdf$legendctype<-legendctype$V1 #add legendctype

#create color bars
CK8colors<-num2col(exprdf$CK8, pal=coolwarm(100))
CK14colors<-num2col(exprdf$CK14, pal=coolwarm(100))
ECADcolors<-num2col(exprdf$ECAD, pal=coolwarm(100))
VIMcolors<-num2col(exprdf$VIM, pal=coolwarm(100))


exprdf$VIMcolors<-VIMcolors
exprdf$ECADcolors<-ECADcolors
exprdf$CK14colors<-CK14colors
exprdf$CK8colors<-CK8colors


legendctypeex<-legendctype[legendctype$maintype!="Other",] #excluding Other subtypes
cols <- cbbPalette
colorlist <- cols[as.numeric(as.factor(legendctypeex$maintype))]

colorlistVIM <- c(c(rep("#A4A4A4",8),rep("#A4A4A4",2),rep("#CC79A7",4),rep("#A4A4A4",4)),exprdf$VIMcolors[19:29])

colorlistCK14 <- c(c(rep("#A4A4A4",8),rep("#A4A4A4",2),rep("#CC79A7",4),rep("#A4A4A4",4)),exprdf$CK14colors[19:29])


#generate visualizations of distance relationships using network qgraph

dev.off()
pdf("Distance_Relationships_rev.pdf",width=8,height=6)

#combining primary and met but stratified by LUM and TNC

#for LUM

xx<-dist_LUM
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlistLUM<-colorlist[as.numeric(colnames(mat_LUM))] #set color
names(colorlistLUM)<-as.character(colnames(mat_LUM)) #set color reference names


g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="LUM 1/yy (from each cluster node) - ALL",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - ALL",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#LUM colored by VIM

colorlistLUM<-colorlistVIM[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="LUM 1/yy (from each cluster node) - ALL - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistLUM<-colorlistVIM[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - ALL - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)


#LUM colored by CK14

colorlistLUM<-colorlistCK14[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="LUM 1/yy (from each cluster node) - ALL - CK14",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistLUM<-colorlistCK14[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - ALL - CK14",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)


#for TNC

xx<-dist_TNC
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),
       as.character(sort(as.numeric(colnames(xx))))]

colorlistTNC<-colorlist[as.numeric(colnames(mat_TNC))] #set color
names(colorlistTNC)<-as.character(colnames(mat_TNC)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="TNC 1/yy (from each cluster node) - ALL",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="TNC 1/t(yy) (to each cluster node) - ALL",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#TNC colored by VIM

colorlistTNC<-colorlistVIM[as.numeric(colnames(mat_TNC))]
names(colorlistTNC)<-as.character(colnames(mat_TNC))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="TNC 1/yy (from each cluster node) - ALL - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistTNC<-colorlistVIM[as.numeric(colnames(mat_TNC))]
names(colorlistTNC)<-as.character(colnames(mat_TNC))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="TNC 1/t(yy) (to each cluster node) - ALL - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#within LUM but primary vs. metastatic

xx<-dist_LUM_pri
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlistLUM<-colorlist[as.numeric(colnames(mat_LUM))] #set color
names(colorlistLUM)<-as.character(colnames(mat_LUM)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="LUM 1/yy (from each cluster node) - PRI",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - PRI",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#LUM colored by VIM

colorlistLUM<-colorlistVIM[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="LUM 1/yy (from each cluster node) - PRI - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistLUM<-colorlistVIM[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - PRI - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#LUM colored by CK14

colorlistLUM<-colorlistCK14[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="LUM 1/yy (from each cluster node) - PRI - CK14",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistLUM<-colorlistCK14[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - PRI - CK14",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)


xx<-dist_LUM_met
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]

colorlistLUM<-colorlist[as.numeric(colnames(mat_LUM))] #set color
names(colorlistLUM)<-as.character(colnames(mat_LUM)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="LUM 1/yy (from each cluster node) - MET",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - MET",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#LUM colored by VIM

colorlistLUM<-colorlistVIM[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="LUM 1/yy (from each cluster node) - MET - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistLUM<-colorlistVIM[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - MET - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#LUM colored by CK14

colorlistLUM<-colorlistCK14[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="LUM 1/yy (from each cluster node) - MET - CK14",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistLUM<-colorlistCK14[as.numeric(colnames(mat_LUM))]
names(colorlistLUM)<-as.character(colnames(mat_LUM))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="LUM 1/t(yy) (to each cluster node) - MET - CK14",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)



#within TNC but primary vs. metastatic

xx<-dist_TNC_pri
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),
       as.character(sort(as.numeric(colnames(xx))))]

colorlistTNC<-colorlist[as.numeric(colnames(mat_TNC))] #set color
names(colorlistTNC)<-as.character(colnames(mat_TNC)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="TNC 1/yy (from each cluster node) - PRI",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="TNC 1/t(yy) (to each cluster node) - PRI",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#TNC colored by VIM

colorlistTNC<-colorlistVIM[as.numeric(colnames(mat_TNC))]
names(colorlistTNC)<-as.character(colnames(mat_TNC))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="TNC 1/yy (from each cluster node) - PRI - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistTNC<-colorlistVIM[as.numeric(colnames(mat_TNC))]
names(colorlistTNC)<-as.character(colnames(mat_TNC))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="TNC 1/t(yy) (to each cluster node) - PRI - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

xx<-dist_TNC_met
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),
       as.character(sort(as.numeric(colnames(xx))))]

colorlistTNC<-colorlist[as.numeric(colnames(mat_TNC))] #set color
names(colorlistTNC)<-as.character(colnames(mat_TNC)) #set color reference names

g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="TNC 1/yy (from each cluster node) - MET",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="TNC 1/t(yy) (to each cluster node) - MET",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"], 
       col = cols , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

#TNC colored by VIM

colorlistTNC<-colorlistVIM[as.numeric(colnames(mat_TNC))]
names(colorlistTNC)<-as.character(colnames(mat_TNC))
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy, 
          title="TNC 1/yy (from each cluster node) - MET - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=F, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

colorlistTNC<-colorlistVIM[as.numeric(colnames(mat_TNC))]
names(colorlistTNC)<-as.character(colnames(mat_TNC))
g<-qgraph(1/t(yy), DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/t(yy), 
          title="TNC 1/t(yy) (to each cluster node) - MET - VIM",
          layout='spring', repulsion=10, 
          maximum=max(1/dist_LUM,1/dist_TNC, na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray", edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistTNC[g$graphAttributes$Nodes$names])
legend("bottomleft", 
       legend=c("Mac","Other Types"), 
       col = c("#CC79A7","#A4A4A4") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.75,
       text.col="black" , horiz = F)

dev.off()

pdf("Distance_Relationships_VIMlegend.pdf",width=3,height=2)
legend.scale(zlim=c(min(exprdf$VIM[19:29]),max(exprdf$VIM[19:29])),
             col=coolwarm(20),
             horizontal = T,
             breaks=seq(min(exprdf$VIM[19:29]),max(exprdf$VIM[19:29]),length.out = 21),
             axis.args = list(at = c(min(exprdf$VIM[19:29]),max(exprdf$VIM[19:29]))))
dev.off()

pdf("Distance_Relationships_CK14legend.pdf",width=3,height=2)
legend.scale(zlim=c(min(exprdf$CK14[19:29]),max(exprdf$CK14[19:29])),
             col=coolwarm(20),
             horizontal = T,
             breaks=seq(min(exprdf$CK14[19:29]),max(exprdf$CK14[19:29]),length.out = 21),
             axis.args = list(at = c(min(exprdf$CK14[19:29]),max(exprdf$CK14[19:29]))))
dev.off()


###DIFF ANALYSIS OF DISTANCES: from Tc (ctype1)

dist_ttest<-function(group1=expr_Prm,
                     group2=expr_NPrm,
                     ctype1="ctype1",
                     ctype2="ctype2"){
  c1toc2_g1 = as.vector(group1[rownames(group1)==ctype1,ctype2])
  c1toc2_g2 = as.vector(group2[rownames(group2)==ctype1,ctype2])
  res<-t.test(x=c1toc2_g1,
              y=c1toc2_g2,
              paired=F)
  medians<-c(median(c1toc2_g1, na.rm=T),median(c1toc2_g2, na.rm=T))
  names(medians)<-c("P_median","NP_median")
  estimate<-res$estimate
  names(estimate)<-c("P_mean","NP_mean")
  diff_med<-medians["P_median"]-medians["NP_median"]
  diff_mean<-estimate["P_mean"]-estimate["NP_mean"]
  combine<-cbind(t(estimate), t(medians), diff_med=diff_med, diff_mean=diff_mean, pvalue=res$p.value)
  rownames(combine)<-paste(ctype1,ctype2,sep="_")
  print(combine)
}


res0<-data.frame(P_median=NA, NP_median=NA, P_mean=NA, NP_mean=NA, diff_med=NA, diff_mean=NA, pvalue=NA)
for(i in 1:length(legendctype$V1)){
  res_save<-dist_ttest(group1=expr_Prm0, group2=expr_NPrm0, ctype1="ctype1", ctype2=as.character(legendctype$V1[i])) 
  res0<-rbind(res0,res_save)
}
res_ctype1<-res0[2:nrow(res0),]  

res_ctype1$padj<-p.adjust(res_ctype1$pvalue, method = "BH")

write.csv(res_ctype1,"Results_distances_Tc.csv")

###VIOLIN PLOTS OF DISTANCES

df_P<-melt(expr_P)
df_P$progression <- "P"
df_NP<-melt(expr_NP)
df_NP$progression <- "NP"
df_combined<-rbind(df_P,df_NP)
colnames(df_combined)<-c("index","target","distance","progression")
df_combined<-df_combined[!is.na(df_combined$distance),]

df_combined_c1ind<-df_combined[df_combined$index=="ctype1",]


pdf('Distance_violinplots.pdf',width=3,height=3)

#Treg cells
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype3",], aes(x=progression, y=distance, fill=progression))+
  ylim(0,500)+ggtitle('Treg')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#B cells
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype4",], aes(x=progression, y=distance, fill=progression))+
  ylim(0,500)+ggtitle('B cells')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Gran
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype5",], aes(x=progression, y=distance, fill=progression))+
  ylim(0,150)+ggtitle('Gran')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mono
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype6",], aes(x=progression, y=distance, fill=progression))+
  ylim(0,500)+ggtitle('Mono')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype7",], aes(x=progression, y=distance, fill=progression))+
  ylim(0,250)+ggtitle('Mac')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

#Mac M2
ggplot(df_combined_c1ind[df_combined_c1ind$target=="ctype8",], aes(x=progression, y=distance, fill=progression))+
  ylim(0,300)+ggtitle('MacM2')+
  geom_violin()+theme_bw()+theme(panel.grid = element_blank(), axis.text = element_text(color="black"))

dev.off()
