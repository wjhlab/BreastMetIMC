####PART 3 TOP NEIGHBOR ANALYSIS

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

####LOAD CLUSTERED DATA####
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


####EXTRACTING MEAN EXPRESSIONS OF KEY MARKERS AND COLORLEGENDS BASED ON EXPRESSION####

##list out the cell types and create a legends data frame

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

##Focused cluster heatmaps

markerlist_mac = c("CD163","CD206","HLADR","CD86","PDL1","ARG1")
exprtbl <- data.frame(fsApply(output$fcs,exprs)[, markerlist_mac], sample_id = output$sample_ids, cluster = output$cell_clustering1m)
exprtbl$type <- factor(output$meta_data$Phenotype[match(exprtbl$sample_id,output$meta_data$sample_id)], levels=typelevels)
exprtbl_lum <- exprtbl[exprtbl$type=="LUM",]
exprtbl_tnc <- exprtbl[exprtbl$type=="TNC",]

exprtbl <- exprtbl %>% group_by(cluster) %>% summarise_at(markerlist_mac, mean, na.rm=TRUE)
exprmtx <- as.matrix(exprtbl[,2:(length(markerlist_mac)+1)])
rownames(exprmtx) <- unlist(exprtbl[,1])

exprtbl_lum <- exprtbl_lum %>% group_by(cluster) %>% summarise_at(markerlist_mac, mean, na.rm=TRUE)
exprmtx_lum <-as.matrix(exprtbl_lum[,2:(length(markerlist_mac)+1)])
rownames(exprmtx_lum) <-unlist(exprtbl_lum[,1])

exprtbl_tnc <- exprtbl_tnc %>% group_by(cluster) %>% summarise_at(markerlist_mac, mean, na.rm=TRUE)
exprmtx_tnc <-as.matrix(exprtbl_tnc[,2:(length(markerlist_mac)+1)])
rownames(exprmtx_tnc) <-unlist(exprtbl_tnc[,1])

pdf("Clusterheatmap_focused_Mac.pdf", width=5, height=5)
pheatmap(exprmtx[clusterlevels[str_detect(clusterlevels,"Mac")],], 
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
pheatmap(exprmtx_lum[clusterlevels[str_detect(clusterlevels,"Mac")],], main = "LUM",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
pheatmap(exprmtx_tnc[clusterlevels[str_detect(clusterlevels,"Mac")],], main = "TNC",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
dev.off()

markerlist_ck = c("CK8","CK14","ECAD","VIM","KI67","HLADR")
exprtbl <- data.frame(fsApply(output$fcs,exprs)[, markerlist_ck], sample_id = output$sample_ids, cluster = output$cell_clustering1m)
exprtbl$type <- factor(output$meta_data$Phenotype[match(exprtbl$sample_id,output$meta_data$sample_id)], levels=typelevels)
exprtbl_lum <- exprtbl[exprtbl$type=="LUM",]
exprtbl_tnc <- exprtbl[exprtbl$type=="TNC",]

exprtbl <- exprtbl %>% group_by(cluster) %>% summarise_at(markerlist_ck, mean, na.rm=TRUE)
exprmtx <-as.matrix(exprtbl[,2:(length(markerlist_ck)+1)])
rownames(exprmtx) <-unlist(exprtbl[,1])

exprtbl_lum <- exprtbl_lum %>% group_by(cluster) %>% summarise_at(markerlist_ck, mean, na.rm=TRUE)
exprmtx_lum <-as.matrix(exprtbl_lum[,2:(length(markerlist_ck)+1)])
rownames(exprmtx_lum) <-unlist(exprtbl_lum[,1])

exprtbl_tnc <- exprtbl_tnc %>% group_by(cluster) %>% summarise_at(markerlist_ck, mean, na.rm=TRUE)
exprmtx_tnc <-as.matrix(exprtbl_tnc[,2:(length(markerlist_ck)+1)])
rownames(exprmtx_tnc) <-unlist(exprtbl_tnc[,1])

pdf("Clusterheatmap_focused_CK.pdf", width=5, height=5)
pheatmap(exprmtx[clusterlevels[str_detect(clusterlevels,"CK")],], main = "ALL",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
pheatmap(exprmtx_lum[clusterlevels[str_detect(clusterlevels,"CK")],], main = "LUM",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
pheatmap(exprmtx_tnc[clusterlevels[str_detect(clusterlevels,"CK")],], main = "TNC",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
dev.off()

markerlist_str = c("COL","PDPN","VIM","SMA","HLADR")
exprtbl <- data.frame(fsApply(output$fcs,exprs)[, markerlist_str], sample_id = output$sample_ids, cluster = output$cell_clustering1m)
exprtbl$type <- factor(output$meta_data$Phenotype[match(exprtbl$sample_id,output$meta_data$sample_id)], levels=typelevels)
exprtbl_lum <- exprtbl[exprtbl$type=="LUM",]
exprtbl_tnc <- exprtbl[exprtbl$type=="TNC",]

exprtbl <- exprtbl %>% group_by(cluster) %>% summarise_at(markerlist_str, mean, na.rm=TRUE)
exprmtx <-as.matrix(exprtbl[,2:(length(markerlist_str)+1)])
rownames(exprmtx) <-unlist(exprtbl[,1])

exprtbl_lum <- exprtbl_lum %>% group_by(cluster) %>% summarise_at(markerlist_str, mean, na.rm=TRUE)
exprmtx_lum <-as.matrix(exprtbl_lum[,2:(length(markerlist_str)+1)])
rownames(exprmtx_lum) <-unlist(exprtbl_lum[,1])

exprtbl_tnc <- exprtbl_tnc %>% group_by(cluster) %>% summarise_at(markerlist_str, mean, na.rm=TRUE)
exprmtx_tnc <-as.matrix(exprtbl_tnc[,2:(length(markerlist_str)+1)])
rownames(exprmtx_tnc) <-unlist(exprtbl_tnc[,1])

pdf("Clusterheatmap_focused_Str.pdf", width=5, height=5)
pheatmap(exprmtx[clusterlevels[str_detect(clusterlevels,"Str")],], main = "ALL",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
pheatmap(exprmtx_lum[clusterlevels[str_detect(clusterlevels,"Str")],], main = "LUM",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
pheatmap(exprmtx_tnc[clusterlevels[str_detect(clusterlevels,"Str")],], main = "TNC",
         scale="column",
         cluster_rows = F,
         cellwidth=10,
         cellheight=10)
dev.off()





##Create scaled color bars
markerlist_ck = c("CK8","CK14","ECAD","VIM","KI67","HLADR")

exprtbl <- data.frame(fsApply(output$fcs,exprs)[, markerlist_ck], sample_id = output$sample_ids, cluster = output$cell_clustering1m)

exprtbl$type <- factor(output$meta_data$Phenotype[match(exprtbl$sample_id,output$meta_data$sample_id)], levels=typelevels)

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

exprdf<-exprdf[exprdf$legendctype %nin% c("ctype30","ctype31"),] #exclude NA and Neurons

####CORR PLOTS####

exprtbl <- data.frame(fsApply(output$fcs,exprs), sample_id = output$sample_ids, cluster = output$cell_clustering1m)
exprtbl$majortype <- exprtbl$cluster
exprtbl$majortype[str_detect(exprtbl$cluster,"Mac")] <- "Mac"
exprtbl$majortype[str_detect(exprtbl$cluster,"CK")] <- "CK"
exprtbl$majortype[str_detect(exprtbl$cluster,"Str")] <- "Str"

exprtbl_mac <- exprtbl[exprtbl$majortype == "Mac",]
exprtbl_mac$type <- factor(output$meta_data$Phenotype[match(exprtbl_mac$sample_id,output$meta_data$sample_id)], levels=typelevels)
exprtbl_mac_lum <- exprtbl_mac[exprtbl_mac$type=="LUM",]
exprtbl_mac_tnc <- exprtbl_mac[exprtbl_mac$type=="TNC",]

corr_mac <- cor(as.matrix(exprtbl_mac[,markerlist_mac]))
corr_mac_lum <- cor(as.matrix(exprtbl_mac_lum[,markerlist_mac]))
corr_mac_tnc <- cor(as.matrix(exprtbl_mac_tnc[,markerlist_mac]))

pdf("Corr_mac.pdf",width=3,height=3)
corrplot(corr_mac, title= "ALL",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
corrplot(corr_mac_lum, title= "LUM",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
corrplot(corr_mac_tnc, title= "TNC",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
dev.off()

exprtbl_ck <- exprtbl[exprtbl$majortype == "CK",]
exprtbl_ck$type <- factor(output$meta_data$Phenotype[match(exprtbl_ck$sample_id,output$meta_data$sample_id)], levels=typelevels)
exprtbl_ck_lum <- exprtbl_ck[exprtbl_ck$type=="LUM",]
exprtbl_ck_tnc <- exprtbl_ck[exprtbl_ck$type=="TNC",]

corr_ck <- cor(as.matrix(exprtbl_ck[,markerlist_ck]))
corr_ck_lum <- cor(as.matrix(exprtbl_ck_lum[,markerlist_ck]))
corr_ck_tnc <- cor(as.matrix(exprtbl_ck_tnc[,markerlist_ck]))

pdf("Corr_ck.pdf",width=3,height=3)
corrplot(corr_ck, title= "ALL",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
corrplot(corr_ck_lum, title= "LUM",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
corrplot(corr_ck_tnc, title= "TNC",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
dev.off()

exprtbl_str <- exprtbl[exprtbl$majortype == "CK",]
exprtbl_str$type <- factor(output$meta_data$Phenotype[match(exprtbl_str$sample_id,output$meta_data$sample_id)], levels=typelevels)
exprtbl_str_lum <- exprtbl_str[exprtbl_str$type=="LUM",]
exprtbl_str_tnc <- exprtbl_str[exprtbl_str$type=="TNC",]

corr_str <- cor(as.matrix(exprtbl_str[,markerlist_str]))
corr_str_lum <- cor(as.matrix(exprtbl_str_lum[,markerlist_str]))
corr_str_tnc <- cor(as.matrix(exprtbl_str_tnc[,markerlist_str]))

pdf("Corr_str.pdf",width=3,height=3)
corrplot(corr_str, title= "ALL",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
corrplot(corr_str_lum, title= "LUM",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
corrplot(corr_str_tnc, title= "TNC",
         method="square", 
         col=rev(brewer.rdbu(100)), 
         type="lower", diag=F, 
         tl.col="black", 
         #addCoef.col = 'black', number.cex = .75,
         addgrid.col= 'black')
dev.off()



####TOP 1-3 NEIGHBORS FROM CK+ CLUSTERS####

#extract expression and neighbors
exprall <- fsApply(output$fcs1, exprs)
exprall <- data.frame(exprall[,c("CK8","CK14","VIM","ECAD","CellId","Num_Neighbors","NN1","NN2","NN3")],
                      objtype=output$cell_clustering1m,
                      sample_id=output$sample_ids)

rownames(exprall)<-paste(exprall$sample_id,exprall$CellId,sep="_")

exprall$NN1 <- paste(exprall$sample_id,exprall$NN1,sep="_")
exprall$NN2 <- paste(exprall$sample_id,exprall$NN2,sep="_")
exprall$NN3 <- paste(exprall$sample_id,exprall$NN3,sep="_")

NN1ind <- match(exprall$NN1,rownames(exprall))
exprall$N1type <- exprall[NN1ind,]$objtype

NN2ind <- match(exprall$NN2,rownames(exprall))
exprall$N2type <- exprall[NN2ind,]$objtype

NN3ind <- match(exprall$NN3,rownames(exprall))
exprall$N3type <- exprall[NN3ind,]$objtype

exprall<-exprall[exprall$objtype!="NA",]

exprall$SPC <- factor(output$meta_data$Case[match(exprall$sample_id,output$meta_data$sample_id)], levels = SPClevels)
exprall$tumor <- factor(output$meta_data$Tumor[match(exprall$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
exprall$type <- factor(output$meta_data$Phenotype[match(exprall$sample_id,output$meta_data$sample_id)], levels=typelevels)
exprall$site <- factor(output$meta_data$Anno[match(exprall$sample_id,output$meta_data$sample_id)], levels=sitelevels)
exprall$TMA <- factor(output$meta_data$TMA[match(exprall$sample_id,output$meta_data$sample_id)], levels=TMAlevels)
exprall$met <- factor(output$meta_data$Met[match(exprall$sample_id,output$meta_data$sample_id)])

#choose the cell types desired for visualizing the neighboring relationships to, from CK+ cells
#double check to change the pdf names below
immunecells<-c("Mac_I","Mac_II","Mac_III","Mac_IV") #focus on macs  #immunecells<-c("Tc","ThN","Treg","NK","B_I","B_II")

##subset for luminal type first

#subsetting luminal
NN_LUM <- exprall[exprall$type=="LUM",]
NN_LUMp <- NN_LUM[NN_LUM$met==0,]
NN_LUMm <- NN_LUM[NN_LUM$met==1,]

NN_LUMm_lung <- NN_LUMm[NN_LUMm$site=="LUNG",]
NN_LUMm_liver <- NN_LUMm[NN_LUMm$site=="LIVER",]
NN_LUMm_brain <- NN_LUMm[NN_LUMm$site=="BRAIN",]

#summarize patients with CKs that have nearest macrophage neighbors

NN_LUMp_CKtoMac <- NN_LUMp[NN_LUMp$objtype %in% c("CK8","CK8_E","CK8_V","CK8_14_E","CK8_14_EV","CK8_14_EV"),] #take only CK cells
NN_LUMp_CKtoMac_N1 <- NN_LUMp_CKtoMac[NN_LUMp_CKtoMac$N1type %in% immunecells,] #get only mac neighbors
NN_LUMp_CKtoMac_N1 <- NN_LUMp_CKtoMac_N1 %>% group_by(SPC, N1type) %>% summarise(tot = table(N1type))
NN_LUMp_CKtoMac_N2 <- NN_LUMp_CKtoMac[NN_LUMp_CKtoMac$N2type %in% immunecells,]
NN_LUMp_CKtoMac_N2 <- NN_LUMp_CKtoMac_N2 %>% group_by(SPC, N2type) %>% summarise(tot = table(N2type))
NN_LUMp_CKtoMac_N3 <- NN_LUMp_CKtoMac[NN_LUMp_CKtoMac$N3type %in% immunecells,]
NN_LUMp_CKtoMac_N3 <- NN_LUMp_CKtoMac_N3 %>% group_by(SPC, N3type) %>% summarise(tot = table(N3type))
colnames(NN_LUMp_CKtoMac_N1) <- c("SPC", "Ntype", "tot")
colnames(NN_LUMp_CKtoMac_N2) <- c("SPC", "Ntype", "tot")
colnames(NN_LUMp_CKtoMac_N3) <- c("SPC", "Ntype", "tot")
NN_LUMp_CKtoMac_Ntop3 <- rbind(NN_LUMp_CKtoMac_N1,NN_LUMp_CKtoMac_N2,NN_LUMp_CKtoMac_N3) #combine all CK to macs
NN_LUMp_CKtoMac_Ntot <- NN_LUMp_CKtoMac_Ntop3 %>% group_by(SPC) %>% summarise(tot = sum(tot)) #total up all neighbors by mac in general
NN_LUMp_CKtoMac_Ntop3 <- NN_LUMp_CKtoMac_Ntop3 %>% group_by(SPC, Ntype) %>% summarise(tot = sum(tot)) #total up all neighbors by mac subtype

NN_LUMp_CKtoMac_tot<-table(NN_LUMp_CKtoMac$SPC) #total CK cells in SPC
NN_LUMp_CKtoMac_tot<-NN_LUMp_CKtoMac_tot[NN_LUMp_CKtoMac_tot!=0] 

colnames(NN_LUMp_CKtoMac_Ntop3) <- c("SPC","Ntype","Ntot")
NN_LUMp_CKtoMac_Ntop3$tot <- NN_LUMp_CKtoMac_tot[as.character(NN_LUMp_CKtoMac_Ntop3$SPC)]
NN_LUMp_CKtoMac_Ntop3$normtot <- NN_LUMp_CKtoMac_Ntop3$Ntot/NN_LUMp_CKtoMac_Ntop3$tot*100 #normalize mac neighbor totals by CK totals
colnames(NN_LUMp_CKtoMac_Ntot) <- c("SPC","Ntot")
NN_LUMp_CKtoMac_Ntot$tot <- NN_LUMp_CKtoMac_tot[as.character(NN_LUMp_CKtoMac_Ntot$SPC)] #for all macs
NN_LUMp_CKtoMac_Ntot$normtot <- NN_LUMp_CKtoMac_Ntot$Ntot/NN_LUMp_CKtoMac_Ntot$tot*100 #normalize

NN_LUMp_CKtoMac_median <- NN_LUMp_CKtoMac_Ntop3 %>% group_by(Ntype) %>% summarise(median = median(normtot)) #find medians
NN_LUMp_CKtoMac_Ntop3$median <- NN_LUMp_CKtoMac_median$median[match(NN_LUMp_CKtoMac_Ntop3$Ntype, NN_LUMp_CKtoMac_median$Ntype)]

bins <- NN_LUMp_CKtoMac_Ntop3$normtot > NN_LUMp_CKtoMac_Ntop3$median #binarize
NN_LUMp_CKtoMac_Ntop3$bins <- bins*1 
write.csv(NN_LUMp_CKtoMac_Ntop3,'Results_NNLUMpCKtoMacbySPC.csv') 

NN_LUMp_CKtoMac_Ntot$bins <- (NN_LUMp_CKtoMac_Ntot$normtot > median(NN_LUMp_CKtoMac_Ntot$normtot))*1 #same for all macs
write.csv(NN_LUMp_CKtoMac_Ntot,'Results_NNLUMpCKtoMactotbySPC.csv') 

#compile data
LUMp1 <- table(NN_LUMp[NN_LUMp$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm1 <- table(NN_LUMm[NN_LUMm$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung1 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8",]) * 10000
LUMm_liver1 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8",]) * 10000
LUMm_brain1 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8",]) * 10000

LUMp2 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_E",]$N1type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm2 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_E",]$N1type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung2 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_E",]$N1type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_E",]) * 10000
LUMm_liver2 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_E",]$N1type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_E",]) * 10000
LUMm_brain2 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_E",]$N1type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_E",]) * 10000

LUMp3 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_V",]$N1type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_V",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm3 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_V",]$N1type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_V",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung3 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_V",]$N1type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_V",]) * 10000
LUMm_liver3 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_V",]$N1type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_V",]) * 10000
LUMm_brain3 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_V",]$N1type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_V",]) * 10000

LUMp4 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_14_E",]$N1type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm4 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_14_E",]$N1type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung4 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_E",]$N1type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_E",]) * 10000
LUMm_liver4 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_E",]$N1type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_E",]) * 10000
LUMm_brain4 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_E",]$N1type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_E",]) * 10000

LUMp5 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_14_EV",]$N1type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm5 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_14_EV",]$N1type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung5 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_EV",]$N1type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_EV",]) * 10000
LUMm_liver5 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_EV",]$N1type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_EV",]) * 10000
LUMm_brain5 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_EV",]$N1type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_EV",]) * 10000

NNdf1 <- data.frame(CK8_pri=as.vector(LUMp1),
                    CK8_E_pri=as.vector(LUMp2),
                    CK8_V_pri=as.vector(LUMp3),
                    CK8_met=as.vector(LUMm1),
                    CK8_met_LUNG=as.vector(LUMm_lung1),
                    CK8_met_LIVER=as.vector(LUMm_liver1),
                    CK8_met_BRAIN=as.vector(LUMm_brain1),
                    CK8_E_met=as.vector(LUMm2),
                    CK8_E_met_LUNG=as.vector(LUMm_lung2),
                    CK8_E_met_LIVER=as.vector(LUMm_liver2),
                    CK8_E_met_BRAIN=as.vector(LUMm_brain2),
                    CK8_V_met=as.vector(LUMm3),
                    CK8_V_met_LUNG=as.vector(LUMm_lung3),
                    CK8_V_met_LIVER=as.vector(LUMm_liver3),
                    CK8_V_met_BRAIN=as.vector(LUMm_brain3),
                    CK8_14_E_pri=as.vector(LUMp4),
                    CK8_14_E_met=as.vector(LUMm4),
                    CK8_14_E_met_LUNG=as.vector(LUMm_lung4),
                    CK8_14_E_met_LIVER=as.vector(LUMm_liver4),
                    CK8_14_E_met_BRAIN=as.vector(LUMm_brain4),
                    CK8_14_EV_pri=as.vector(LUMp5),
                    CK8_14_EV_met=as.vector(LUMm5),
                    CK8_14_EV_met_LUNG=as.vector(LUMm_lung5),
                    CK8_14_EV_met_LIVER=as.vector(LUMm_liver5),
                    CK8_14_EV_met_BRAIN=as.vector(LUMm_brain5))


rownames(NNdf1) <- clusterlevels
NNdf1<-NNdf1[immunecells,]
rownames(NNdf1) <- paste("N1",rownames(NNdf1),sep="_")
NNdf1[is.na(NNdf1)]<-0

LUMp1 <- table(NN_LUMp[NN_LUMp$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm1 <- table(NN_LUMm[NN_LUMm$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung1 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8",]) * 10000
LUMm_liver1 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8",]) * 10000
LUMm_brain1 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8",]) * 10000

LUMp2 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_E",]$N2type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm2 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_E",]$N2type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung2 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_E",]$N2type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_E",]) * 10000
LUMm_liver2 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_E",]$N2type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_E",]) * 10000
LUMm_brain2 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_E",]$N2type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_E",]) * 10000

LUMp3 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_V",]$N2type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_V",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm3 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_V",]$N2type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_V",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung3 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_V",]$N2type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_V",]) * 10000
LUMm_liver3 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_V",]$N2type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_V",]) * 10000
LUMm_brain3 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_V",]$N2type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_V",]) * 10000

LUMp4 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_14_E",]$N2type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm4 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_14_E",]$N2type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung4 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_E",]$N2type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_E",]) * 10000
LUMm_liver4 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_E",]$N2type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_E",]) * 10000
LUMm_brain4 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_E",]$N2type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_E",]) * 10000

LUMp5 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_14_EV",]$N2type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm5 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_14_EV",]$N2type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung5 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_EV",]$N2type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_EV",]) * 10000
LUMm_liver5 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_EV",]$N2type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_EV",]) * 10000
LUMm_brain5 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_EV",]$N2type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_EV",]) * 10000

NNdf2 <- data.frame(CK8_pri=as.vector(LUMp1),
                    CK8_E_pri=as.vector(LUMp2),
                    CK8_V_pri=as.vector(LUMp3),
                    CK8_met=as.vector(LUMm1),
                    CK8_met_LUNG=as.vector(LUMm_lung1),
                    CK8_met_LIVER=as.vector(LUMm_liver1),
                    CK8_met_BRAIN=as.vector(LUMm_brain1),
                    CK8_E_met=as.vector(LUMm2),
                    CK8_E_met_LUNG=as.vector(LUMm_lung2),
                    CK8_E_met_LIVER=as.vector(LUMm_liver2),
                    CK8_E_met_BRAIN=as.vector(LUMm_brain2),
                    CK8_V_met=as.vector(LUMm3),
                    CK8_V_met_LUNG=as.vector(LUMm_lung3),
                    CK8_V_met_LIVER=as.vector(LUMm_liver3),
                    CK8_V_met_BRAIN=as.vector(LUMm_brain3),
                    CK8_14_E_pri=as.vector(LUMp4),
                    CK8_14_E_met=as.vector(LUMm4),
                    CK8_14_E_met_LUNG=as.vector(LUMm_lung4),
                    CK8_14_E_met_LIVER=as.vector(LUMm_liver4),
                    CK8_14_E_met_BRAIN=as.vector(LUMm_brain4),
                    CK8_14_EV_pri=as.vector(LUMp5),
                    CK8_14_EV_met=as.vector(LUMm5),
                    CK8_14_EV_met_LUNG=as.vector(LUMm_lung5),
                    CK8_14_EV_met_LIVER=as.vector(LUMm_liver5),
                    CK8_14_EV_met_BRAIN=as.vector(LUMm_brain5))


rownames(NNdf2) <- clusterlevels
NNdf2<-NNdf2[immunecells,]
rownames(NNdf2) <- paste("N2",rownames(NNdf2),sep="_")
NNdf2[is.na(NNdf2)]<-0

LUMp1 <- table(NN_LUMp[NN_LUMp$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm1 <- table(NN_LUMm[NN_LUMm$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung1 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8",]) * 10000
LUMm_liver1 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8",]) * 10000
LUMm_brain1 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8",]) * 10000

LUMp2 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_E",]$N3type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm2 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_E",]$N3type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung2 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_E",]$N3type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_E",]) * 10000
LUMm_liver2 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_E",]$N3type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_E",]) * 10000
LUMm_brain2 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_E",]$N3type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_E",]) * 10000

LUMp3 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_V",]$N3type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_V",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm3 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_V",]$N3type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_V",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung3 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_V",]$N3type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_V",]) * 10000
LUMm_liver3 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_V",]$N3type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_V",]) * 10000
LUMm_brain3 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_V",]$N3type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_V",]) * 10000

LUMp4 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_14_E",]$N3type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm4 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_14_E",]$N3type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung4 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_E",]$N3type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_E",]) * 10000
LUMm_liver4 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_E",]$N3type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_E",]) * 10000
LUMm_brain4 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_E",]$N3type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_E",]) * 10000

LUMp5 <- table(NN_LUMp[NN_LUMp$objtype=="CK8_14_EV",]$N3type)[clusterlevels] / nrow(NN_LUMp[NN_LUMp$objtype=="CK8_14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm5 <- table(NN_LUMm[NN_LUMm$objtype=="CK8_14_EV",]$N3type)[clusterlevels] / nrow(NN_LUMm[NN_LUMm$objtype=="CK8_14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
LUMm_lung5 <- table(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_EV",]$N3type)[clusterlevels] / nrow(NN_LUMm_lung[NN_LUMm_lung$objtype=="CK8_14_EV",]) * 10000
LUMm_liver5 <- table(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_EV",]$N3type)[clusterlevels] / nrow(NN_LUMm_liver[NN_LUMm_liver$objtype=="CK8_14_EV",]) * 10000
LUMm_brain5 <- table(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_EV",]$N3type)[clusterlevels] / nrow(NN_LUMm_brain[NN_LUMm_brain$objtype=="CK8_14_EV",]) * 10000

NNdf3 <- data.frame(CK8_pri=as.vector(LUMp1),
                    CK8_E_pri=as.vector(LUMp2),
                    CK8_V_pri=as.vector(LUMp3),
                    CK8_met=as.vector(LUMm1),
                    CK8_met_LUNG=as.vector(LUMm_lung1),
                    CK8_met_LIVER=as.vector(LUMm_liver1),
                    CK8_met_BRAIN=as.vector(LUMm_brain1),
                    CK8_E_met=as.vector(LUMm2),
                    CK8_E_met_LUNG=as.vector(LUMm_lung2),
                    CK8_E_met_LIVER=as.vector(LUMm_liver2),
                    CK8_E_met_BRAIN=as.vector(LUMm_brain2),
                    CK8_V_met=as.vector(LUMm3),
                    CK8_V_met_LUNG=as.vector(LUMm_lung3),
                    CK8_V_met_LIVER=as.vector(LUMm_liver3),
                    CK8_V_met_BRAIN=as.vector(LUMm_brain3),
                    CK8_14_E_pri=as.vector(LUMp4),
                    CK8_14_E_met=as.vector(LUMm4),
                    CK8_14_E_met_LUNG=as.vector(LUMm_lung4),
                    CK8_14_E_met_LIVER=as.vector(LUMm_liver4),
                    CK8_14_E_met_BRAIN=as.vector(LUMm_brain4),
                    CK8_14_EV_pri=as.vector(LUMp5),
                    CK8_14_EV_met=as.vector(LUMm5),
                    CK8_14_EV_met_LUNG=as.vector(LUMm_lung5),
                    CK8_14_EV_met_LIVER=as.vector(LUMm_liver5),
                    CK8_14_EV_met_BRAIN=as.vector(LUMm_brain5))

rownames(NNdf3) <- clusterlevels
NNdf3<-NNdf3[immunecells,]
rownames(NNdf3) <- paste("N3",rownames(NNdf3),sep="_")
NNdf3[is.na(NNdf3)]<-0

#collate dataframes for all top 3 neighbors
NNdf<-rbind(NNdf1,NNdf2,NNdf3)

#clean up low count data (anything with less than 1000 count in the subset can be removed)- primary has plenty

mfilterselect<-function(objtype=objtype, 
                        minnum=1000, 
                        tissuetypes=tissuetypes){
  tofilter<-names(which(table(NN_LUMm[NN_LUMm$objtype==objtype,]$site)<minnum))[names(which(table(NN_LUMm[NN_LUMm$objtype==objtype,]$site)<1000)) %in% tissuetypes]
  ifelse(tofilter!=0,return(paste(objtype,"met",tofilter,sep="_")),return("nothing_to_filter"))
  #return(paste(objtype,"met",tofilter,sep="_"))
}

totest<-list("CK8","CK8_E","CK8_V","CK8_14_E","CK8_14_EV","CK8_14_EV")
lowcountsexclude<-c()
for(i in 1:length(totest)){
  lowcountsexclude<-c(lowcountsexclude,mfilterselect(totest[i],1000,c("LUNG","LIVER","BRAIN")))}

colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]

NNdf<-NNdf[,colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]]

#annotations

markerstoannotate <- c("VIM","ECAD","CK8","CK14")

anno_expr_p<-NN_LUMp %>% group_by(objtype,met) %>% summarise_at(vars(markerstoannotate), mean)
anno_expr_p<-anno_expr_p[anno_expr_p$objtype %in% c("CK8","CK8_E","CK8_V","CK8_14_E","CK8_14_EV"),]
anno_expr_p$site<- "PBC"
anno_expr_p$pri<- "PRI"
anno_expr_p$name<- paste(anno_expr_p$objtype,"pri",sep="_")
anno_expr_m<-NN_LUMm %>% group_by(objtype,met,site) %>% summarise_at(vars(markerstoannotate), mean)
anno_expr_m<-anno_expr_m[anno_expr_m$objtype %in% c("CK8","CK8_E","CK8_V","CK8_14_E","CK8_14_EV"),]
anno_expr_m<-anno_expr_m[anno_expr_m$site %in% c("LUNG","LIVER","BRAIN"),]
anno_expr_m$pri<- "MET"
anno_expr_m$name <- paste(anno_expr_m$objtype,"met",anno_expr_m$site,sep="_")
anno_expr_m2<-anno_expr_m %>% group_by(objtype,met) %>% summarise_at(vars(markerstoannotate), mean)
anno_expr_m2$site <- "MET"
anno_expr_m2$pri<- "MET"
anno_expr_m2$name <- paste(anno_expr_m2$objtype,"met",sep="_")
anno_expr<-rbind(anno_expr_p,
                 anno_expr_m,
                 anno_expr_m2)

ck8<-anno_expr$CK8
ck14<-anno_expr$CK14
vim<-anno_expr$VIM
ecad<-anno_expr$ECAD
pri<-anno_expr$pri
site<-anno_expr$site
names<-anno_expr$name

colannos<-data.frame(vim=vim,
                     ecad=ecad,
                     ck8=ck8,
                     ck14=ck14,
                     met=pri,
                     site=site,
                     row.names=names)
colannos_lum <- colannos

colannos_colors_p <- c("black","lightgray")
names(colannos_colors_p) <- c("PRI","MET")
colannos_colors_s <- brewer.set1(5)
names(colannos_colors_s) <- c("PBC","MET","LUNG","LIVER","BRAIN")
colannos_colors_e <- num2col(colannos$ecad, colorRampPalette(c("black","yellow"))(100))
names(colannos_colors_e) <- colannos$ecad
colannos_colors_v <- num2col(colannos$vim, colorRampPalette(c("black","yellow"))(100))
names(colannos_colors_v) <- colannos$vim
colannos_colors_8 <- num2col(colannos$ck8, colorRampPalette(c("black","green"))(100), ref=c(0,6.5)) #note range set
names(colannos_colors_8) <- colannos$ck8
colannos_colors_14 <- num2col(colannos$ck14, colorRampPalette(c("black","cyan"))(100), ref=c(0,5.2)) #note range set
names(colannos_colors_14) <- colannos$ck14
colannos_colors <- list(met=colannos_colors_p, 
                        site=colannos_colors_s,
                        ecad=sort(colannos_colors_e), 
                        vim=sort(colannos_colors_v),
                        ck8=sort(colannos_colors_8),
                        ck14=sort(colannos_colors_14))

#draw heatmap
NNdf_lum<-NNdf[rowSums(NNdf)!=0,] #get rid of rows with all zeros
colannos_lum<-colannos_lum[rownames(colannos_lum) %in% colnames(NNdf_lum),]

pdf("NN_LUM_CK8_N1-3_Mac.pdf", width=7, height=12)
pheatmap(NNdf_lum, 
         main="Top Mac Neighbors of Luminal Cancer CK8+ Clusters",
         cluster_rows = T, cluster_columns = F, scale = "row", 
         annotation_col = colannos_lum, 
         annotation_colors = colannos_colors,
         cutree_rows = 4,
         color = rev(brewer.rdbu(100)),
         border_color = NA,
         treeheight_row = 5,
         treeheight_col = 5,
         cellheight = 10, cellwidth = 10,
         show_colnames = F)
dev.off()


##repeat for TNC

#subset TNC
NN_TNC <- exprall[exprall$type=="TNC",]
NN_TNCp <- NN_TNC[NN_TNC$met==0,]
NN_TNCm <- NN_TNC[NN_TNC$met==1,]

NN_TNCm_lung <- NN_TNCm[NN_TNCm$site=="LUNG",]
NN_TNCm_liver <- NN_TNCm[NN_TNCm$site=="LIVER",]
NN_TNCm_brain <- NN_TNCm[NN_TNCm$site=="BRAIN",]

#summarize patients with CKs that have nearest macrophage neighbors

NN_TNCp_CKtoMac <- NN_TNCp[NN_TNCp$objtype %in% c("CK8","CK8_E","CK8_V","CK8_14_E","CK8_14_EV","CK8_14_EV"),] #take only CK cells
NN_TNCp_CKtoMac_N1 <- NN_TNCp_CKtoMac[NN_TNCp_CKtoMac$N1type %in% immunecells,] #get only mac neighbors
NN_TNCp_CKtoMac_N1 <- NN_TNCp_CKtoMac_N1 %>% group_by(SPC, N1type) %>% summarise(tot = table(N1type))
NN_TNCp_CKtoMac_N2 <- NN_TNCp_CKtoMac[NN_TNCp_CKtoMac$N2type %in% immunecells,]
NN_TNCp_CKtoMac_N2 <- NN_TNCp_CKtoMac_N2 %>% group_by(SPC, N2type) %>% summarise(tot = table(N2type))
NN_TNCp_CKtoMac_N3 <- NN_TNCp_CKtoMac[NN_TNCp_CKtoMac$N3type %in% immunecells,]
NN_TNCp_CKtoMac_N3 <- NN_TNCp_CKtoMac_N3 %>% group_by(SPC, N3type) %>% summarise(tot = table(N3type))
colnames(NN_TNCp_CKtoMac_N1) <- c("SPC", "Ntype", "tot")
colnames(NN_TNCp_CKtoMac_N2) <- c("SPC", "Ntype", "tot")
colnames(NN_TNCp_CKtoMac_N3) <- c("SPC", "Ntype", "tot")
NN_TNCp_CKtoMac_Ntop3 <- rbind(NN_TNCp_CKtoMac_N1,NN_TNCp_CKtoMac_N2,NN_TNCp_CKtoMac_N3) #combine all CK to macs
NN_TNCp_CKtoMac_Ntot <- NN_TNCp_CKtoMac_Ntop3 %>% group_by(SPC) %>% summarise(tot = sum(tot)) #total up all neighbors by mac in general
NN_TNCp_CKtoMac_Ntop3 <- NN_TNCp_CKtoMac_Ntop3 %>% group_by(SPC, Ntype) %>% summarise(tot = sum(tot)) #total up all neighbors by mac subtype

NN_TNCp_CKtoMac_tot<-table(NN_TNCp_CKtoMac$SPC) #total CK cells in SPC
NN_TNCp_CKtoMac_tot<-NN_TNCp_CKtoMac_tot[NN_TNCp_CKtoMac_tot!=0] 

colnames(NN_TNCp_CKtoMac_Ntop3) <- c("SPC","Ntype","Ntot")
NN_TNCp_CKtoMac_Ntop3$tot <- NN_TNCp_CKtoMac_tot[as.character(NN_TNCp_CKtoMac_Ntop3$SPC)]
NN_TNCp_CKtoMac_Ntop3$normtot <- NN_TNCp_CKtoMac_Ntop3$Ntot/NN_TNCp_CKtoMac_Ntop3$tot*100 #normalize mac neighbor totals by CK totals
colnames(NN_TNCp_CKtoMac_Ntot) <- c("SPC","Ntot")
NN_TNCp_CKtoMac_Ntot$tot <- NN_TNCp_CKtoMac_tot[as.character(NN_TNCp_CKtoMac_Ntot$SPC)] #for all macs
NN_TNCp_CKtoMac_Ntot$normtot <- NN_TNCp_CKtoMac_Ntot$Ntot/NN_TNCp_CKtoMac_Ntot$tot*100 #normalize

NN_TNCp_CKtoMac_median <- NN_TNCp_CKtoMac_Ntop3 %>% group_by(Ntype) %>% summarise(median = median(normtot)) #find medians
NN_TNCp_CKtoMac_Ntop3$median <- NN_TNCp_CKtoMac_median$median[match(NN_TNCp_CKtoMac_Ntop3$Ntype, NN_TNCp_CKtoMac_median$Ntype)]

bins <- NN_TNCp_CKtoMac_Ntop3$normtot > NN_TNCp_CKtoMac_Ntop3$median #binarize
NN_TNCp_CKtoMac_Ntop3$bins <- bins*1 
write.csv(NN_TNCp_CKtoMac_Ntop3,'Results_NNTNCpCKtoMacbySPC.csv') 

NN_TNCp_CKtoMac_Ntot$bins <- (NN_TNCp_CKtoMac_Ntot$normtot > median(NN_TNCp_CKtoMac_Ntot$normtot))*1 #same for all macs
write.csv(NN_TNCp_CKtoMac_Ntot,'Results_NNTNCpCKtoMactotbySPC.csv') 

#compile data
TNCp1 <- table(NN_TNCp[NN_TNCp$objtype=="CK14",]$N1type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm1 <- table(NN_TNCm[NN_TNCm$objtype=="CK14",]$N1type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung1 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14",]$N1type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14",]) * 10000
TNCm_brain1 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14",]$N1type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14",]) * 10000

TNCp2 <- table(NN_TNCp[NN_TNCp$objtype=="CK14_E",]$N1type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm2 <- table(NN_TNCm[NN_TNCm$objtype=="CK14_E",]$N1type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung2 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_E",]$N1type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_E",]) * 10000
TNCm_brain2 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_E",]$N1type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_E",]) * 10000

TNCp3 <- table(NN_TNCp[NN_TNCp$objtype=="CK14_EV",]$N1type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm3 <- table(NN_TNCm[NN_TNCm$objtype=="CK14_EV",]$N1type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung3 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_EV",]$N1type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_EV",]) * 10000
TNCm_brain3 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_EV",]$N1type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_EV",]) * 10000

TNCp4 <- table(NN_TNCp[NN_TNCp$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm4 <- table(NN_TNCm[NN_TNCm$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung4 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK8",]) * 10000
TNCm_brain4 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK8",]$N1type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK8",]) * 10000

NNdf1 <- data.frame(CK14_pri=as.vector(TNCp1),
                    CK14_E_pri=as.vector(TNCp2),
                    CK14_EV_pri=as.vector(TNCp3),
                    CK14_met=as.vector(TNCm1),
                    CK14_met_LUNG=as.vector(TNCm_lung1),
                    CK14_met_BRAIN=as.vector(TNCm_brain1),
                    CK14_E_met=as.vector(TNCm2),
                    CK14_E_met_LUNG=as.vector(TNCm_lung2),
                    CK14_E_met_BRAIN=as.vector(TNCm_brain2),
                    CK14_EV_met=as.vector(TNCm3),
                    CK14_EV_met_LUNG=as.vector(TNCm_lung3),
                    CK14_EV_met_BRAIN=as.vector(TNCm_brain3),
                    CK8_pri=as.vector(TNCp4),
                    CK8_met=as.vector(TNCm4),
                    CK8_met_LUNG=as.vector(TNCm_lung4),
                    CK8_met_BRAIN=as.vector(TNCm_brain4))

rownames(NNdf1) <- clusterlevels
NNdf1<-NNdf1[immunecells,]
rownames(NNdf1) <- paste("N1",rownames(NNdf1),sep="_")
NNdf1[is.na(NNdf1)]<-0

TNCp1 <- table(NN_TNCp[NN_TNCp$objtype=="CK14",]$N2type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm1 <- table(NN_TNCm[NN_TNCm$objtype=="CK14",]$N2type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung1 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14",]$N2type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14",]) * 10000
TNCm_brain1 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14",]$N2type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14",]) * 10000

TNCp2 <- table(NN_TNCp[NN_TNCp$objtype=="CK14_E",]$N2type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm2 <- table(NN_TNCm[NN_TNCm$objtype=="CK14_E",]$N2type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung2 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_E",]$N2type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_E",]) * 10000
TNCm_brain2 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_E",]$N2type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_E",]) * 10000

TNCp3 <- table(NN_TNCp[NN_TNCp$objtype=="CK14_EV",]$N2type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm3 <- table(NN_TNCm[NN_TNCm$objtype=="CK14_EV",]$N2type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung3 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_EV",]$N2type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_EV",]) * 10000
TNCm_brain3 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_EV",]$N2type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_EV",]) * 10000

TNCp4 <- table(NN_TNCp[NN_TNCp$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm4 <- table(NN_TNCm[NN_TNCm$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung4 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK8",]) * 10000
TNCm_brain4 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK8",]$N2type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK8",]) * 10000

NNdf2 <- data.frame(CK14_pri=as.vector(TNCp1),
                    CK14_E_pri=as.vector(TNCp2),
                    CK14_EV_pri=as.vector(TNCp3),
                    CK14_met=as.vector(TNCm1),
                    CK14_met_LUNG=as.vector(TNCm_lung1),
                    CK14_met_BRAIN=as.vector(TNCm_brain1),
                    CK14_E_met=as.vector(TNCm2),
                    CK14_E_met_LUNG=as.vector(TNCm_lung2),
                    CK14_E_met_BRAIN=as.vector(TNCm_brain2),
                    CK14_EV_met=as.vector(TNCm3),
                    CK14_EV_met_LUNG=as.vector(TNCm_lung3),
                    CK14_EV_met_BRAIN=as.vector(TNCm_brain3),
                    CK8_pri=as.vector(TNCp4),
                    CK8_met=as.vector(TNCm4),
                    CK8_met_LUNG=as.vector(TNCm_lung4),
                    CK8_met_BRAIN=as.vector(TNCm_brain4))

rownames(NNdf2) <- clusterlevels
NNdf2<-NNdf2[immunecells,]
rownames(NNdf2) <- paste("N2",rownames(NNdf2),sep="_")
NNdf2[is.na(NNdf2)]<-0

TNCp1 <- table(NN_TNCp[NN_TNCp$objtype=="CK14",]$N3type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm1 <- table(NN_TNCm[NN_TNCm$objtype=="CK14",]$N3type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung1 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14",]$N3type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14",]) * 10000
TNCm_brain1 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14",]$N3type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14",]) * 10000

TNCp2 <- table(NN_TNCp[NN_TNCp$objtype=="CK14_E",]$N3type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm2 <- table(NN_TNCm[NN_TNCm$objtype=="CK14_E",]$N3type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14_E",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung2 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_E",]$N3type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_E",]) * 10000
TNCm_brain2 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_E",]$N3type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_E",]) * 10000

TNCp3 <- table(NN_TNCp[NN_TNCp$objtype=="CK14_EV",]$N3type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm3 <- table(NN_TNCm[NN_TNCm$objtype=="CK14_EV",]$N3type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK14_EV",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung3 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_EV",]$N3type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK14_EV",]) * 10000
TNCm_brain3 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_EV",]$N3type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK14_EV",]) * 10000

TNCp4 <- table(NN_TNCp[NN_TNCp$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_TNCp[NN_TNCp$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm4 <- table(NN_TNCm[NN_TNCm$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_TNCm[NN_TNCm$objtype=="CK8",]) * 10000 # number of first nearest neighbor for every 10K cells (normalized by total first)
TNCm_lung4 <- table(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_TNCm_lung[NN_TNCm_lung$objtype=="CK8",]) * 10000
TNCm_brain4 <- table(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK8",]$N3type)[clusterlevels] / nrow(NN_TNCm_brain[NN_TNCm_brain$objtype=="CK8",]) * 10000

NNdf3 <- data.frame(CK14_pri=as.vector(TNCp1),
                    CK14_E_pri=as.vector(TNCp2),
                    CK14_EV_pri=as.vector(TNCp3),
                    CK14_met=as.vector(TNCm1),
                    CK14_met_LUNG=as.vector(TNCm_lung1),
                    CK14_met_BRAIN=as.vector(TNCm_brain1),
                    CK14_E_met=as.vector(TNCm2),
                    CK14_E_met_LUNG=as.vector(TNCm_lung2),
                    CK14_E_met_BRAIN=as.vector(TNCm_brain2),
                    CK14_EV_met=as.vector(TNCm3),
                    CK14_EV_met_LUNG=as.vector(TNCm_lung3),
                    CK14_EV_met_BRAIN=as.vector(TNCm_brain3),
                    CK8_pri=as.vector(TNCp4),
                    CK8_met=as.vector(TNCm4),
                    CK8_met_LUNG=as.vector(TNCm_lung4),
                    CK8_met_BRAIN=as.vector(TNCm_brain4))

rownames(NNdf3) <- clusterlevels
NNdf3<-NNdf3[immunecells,]
rownames(NNdf3) <- paste("N3",rownames(NNdf3),sep="_")
NNdf3[is.na(NNdf3)]<-0

#collate dataframes for all top 3 neighbors
NNdf<-rbind(NNdf1,NNdf2,NNdf3)

#clean up low count data (anything with less than 1000 count in the subset can be removed)- primary has plenty

mfilterselect<-function(objtype=objtype, 
                        minnum=1000, 
                        tissuetypes=tissuetypes){
  tofilter<-names(which(table(NN_LUMm[NN_LUMm$objtype==objtype,]$site)<minnum))[names(which(table(NN_LUMm[NN_LUMm$objtype==objtype,]$site)<1000)) %in% tissuetypes]
  ifelse(tofilter!=0,return(paste(objtype,"met",tofilter,sep="_")),return("nothing_to_filter"))
  #return(paste(objtype,"met",tofilter,sep="_"))
}

totest<-list("CK14","CK14_E","CK14_EV","CK8")
lowcountsexclude<-c()
for(i in 1:length(totest)){
  lowcountsexclude<-c(lowcountsexclude,mfilterselect(totest[i],1000,c("LUNG","BRAIN")))}

colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]

NNdf<-NNdf[,colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]]


#annotations

markerstoannotate <- c("VIM","ECAD","CK8","CK14")

anno_expr_p<-NN_TNCp %>% group_by(objtype,met) %>% summarise_at(vars(all_of(markerstoannotate)), mean)
anno_expr_p<-anno_expr_p[anno_expr_p$objtype %in% c("CK14","CK14_E","CK14_EV","CK8"),]
anno_expr_p$site<- "PBC"
anno_expr_p$pri<- "PRI"
anno_expr_p$name<- paste(anno_expr_p$objtype,"pri",sep="_")
anno_expr_m<-NN_TNCm %>% group_by(objtype,met,site) %>% summarise_at(vars(markerstoannotate), mean)
anno_expr_m<-anno_expr_m[anno_expr_m$objtype %in% c("CK14","CK14_E","CK14_EV","CK8"),]
anno_expr_m<-anno_expr_m[anno_expr_m$site %in% c("LUNG","BRAIN"),]
anno_expr_m$pri<- "MET"
anno_expr_m$name <- paste(anno_expr_m$objtype,"met",anno_expr_m$site,sep="_")
anno_expr_m2<-anno_expr_m %>% group_by(objtype,met) %>% summarise_at(vars(markerstoannotate), mean)
anno_expr_m2$site <- "MET"
anno_expr_m2$pri<- "MET"
anno_expr_m2$name <- paste(anno_expr_m2$objtype,"met",sep="_")
anno_expr<-rbind(anno_expr_p,
                 anno_expr_m,
                 anno_expr_m2)

ck8<-anno_expr$CK8
ck14<-anno_expr$CK14
vim<-anno_expr$VIM
ecad<-anno_expr$ECAD
pri<-anno_expr$pri
site<-anno_expr$site
names<-anno_expr$name

colannos<-data.frame(vim=vim,
                     ecad=ecad,
                     ck8=ck8,
                     ck14=ck14,
                     met=pri,
                     site=site,
                     row.names=names)
colannos_tnc<-colannos

colannos_colors_p <- c("black","lightgray")
names(colannos_colors_p) <- c("PRI","MET")
colannos_colors_s <- brewer.set1(5)[c(1:3,5)]
names(colannos_colors_s) <- c("PBC","MET","LUNG","BRAIN")
colannos_colors_e <- num2col(colannos$ecad, colorRampPalette(c("black","yellow"))(100))
names(colannos_colors_e) <- colannos$ecad
colannos_colors_v <- num2col(colannos$vim, colorRampPalette(c("black","yellow"))(100))
names(colannos_colors_v) <- colannos$vim
colannos_colors_8 <- num2col(colannos$ck8, colorRampPalette(c("black","green"))(100), ref=c(0,6.5)) #note range set
names(colannos_colors_8) <- colannos$ck8
colannos_colors_14 <- num2col(colannos$ck14, colorRampPalette(c("black","cyan"))(100), ref=c(0,5.2)) #note range set
names(colannos_colors_14) <- colannos$ck14
colannos_colors <- list(met=colannos_colors_p, 
                        site=colannos_colors_s,
                        ecad=sort(colannos_colors_e), 
                        vim=sort(colannos_colors_v),
                        ck8=sort(colannos_colors_8),
                        ck14=sort(colannos_colors_14))

#draw heatmap
NNdf_tnc<-NNdf[rowSums(NNdf)!=0,] #get rid of rows with all zeros
colannos_tnc<-colannos_tnc[rownames(colannos_tnc) %in% colnames(NNdf_tnc),]

pdf("NN_TNC_CK14_N1-3_Mac.pdf", width=7, height=14)
pheatmap(NNdf_tnc, 
         main="Top Mac Neighbors of TNC Cancer CK14+ Clusters",
         cluster_rows = T, cluster_columns = F, scale = "row", 
         annotation_col = colannos_tnc, 
         annotation_colors = colannos_colors,
         color = rev(brewer.rdbu(100)),
         cutree_rows = 4,
         border_color = NA,
         treeheight_row = 5,
         treeheight_col = 5,
         cellheight = 10, cellwidth = 10,
         show_colnames = F)
dev.off()


###JOINT HEATMAP

colnames(NNdf_lum)<-paste0(colnames(NNdf_lum),"_lum")
colnames(NNdf_tnc)<-paste0(colnames(NNdf_tnc),"_tnc")
NNdf_all <- cbind(NNdf_lum,NNdf_tnc)

colannos_lum$type<-factor("LUM")
rownames(colannos_lum)<-paste0(rownames(colannos_lum),"_lum")
colannos_tnc$type<-factor("TNC")
rownames(colannos_tnc)<-paste0(rownames(colannos_tnc),"_tnc")
colannos_all<-rbind(colannos_lum,colannos_tnc)
colannos_all<-colannos_all[,c("vim","ecad","ck8","ck14","met","type","site")] #rearrange the order in which the annos show

colannos_colors_p <- c("black","lightgray")
names(colannos_colors_p) <- c("PRI","MET")
colannos_colors_s <- brewer.set1(5)
names(colannos_colors_s) <- c("PBC","MET","LUNG","LIVER","BRAIN")
colannos_colors_e <- num2col(colannos_all$ecad, colorRampPalette(c("black","yellow"))(100))
names(colannos_colors_e) <- colannos_all$ecad
colannos_colors_v <- num2col(colannos_all$vim, colorRampPalette(c("black","yellow"))(100))
names(colannos_colors_v) <- colannos_all$vim
colannos_colors_8 <- num2col(colannos_all$ck8, colorRampPalette(c("black","green"))(100)) 
names(colannos_colors_8) <- colannos_all$ck8
colannos_colors_14 <- num2col(colannos_all$ck14, colorRampPalette(c("black","cyan"))(100)) 
names(colannos_colors_14) <- colannos_all$ck14
colannos_colors_type <- c("black","lightgray")
names(colannos_colors_type)<-c("LUM","TNC")

colannos_colors_joint <- list(met=colannos_colors_p, 
                              site=colannos_colors_s,
                              type=colannos_colors_type,
                              ecad=sort(colannos_colors_e), 
                              vim=sort(colannos_colors_v),
                              ck8=sort(colannos_colors_8),
                              ck14=sort(colannos_colors_14))

pdf("NN_all_N1-3_Mac.pdf", width=7, height=14)
pheatmap(NNdf_all, 
         main="Top Mac Neighbors of CK+ Clusters
         
         ",
         cluster_rows = T, cluster_columns = F, scale = "row", 
         annotation_col = colannos_all, 
         annotation_colors = colannos_colors_joint,
         color = rev(brewer.rdbu(100)),
         cutree_rows = 4,
         border_color = NA,
         treeheight_row = 5,
         treeheight_col = 5,
         cellheight = 15, cellwidth = 10,
         show_colnames = F)
dev.off()

