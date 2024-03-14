####PART 4 DISTANCE RELATIONSHIPS

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
library(spatstat)

####LOAD DATA####
output<-readRDS('backup_output.rds')
clusterMergeFile = paste0(work,"/Config/merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)
counts_table <- table(output$cell_clustering1m, output$sample_ids)
counts <- as.data.frame.matrix(counts_table)

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

####DISTANCE RELATIONSHIPS####

#how many of each cell types are there in the dataset?

ggdft <- melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
              id.vars = "cluster", value.name = "counts", 
              variable.name = "sample_id")
ggdft$sample_id <- factor(ggdft$sample_id, levels=samplevels)
ggdft$cluster <- factor(ggdft$cluster, levels=clusterlevels)
ggdft$SPC <- factor(output$meta_data$Case[match(ggdft$sample_id,output$meta_data$sample_id)])
ggdft$tumor <- factor(output$meta_data$Tumor[match(ggdft$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
ggdft$type <- factor(output$meta_data$Phenotype[match(ggdft$sample_id,output$meta_data$sample_id)], levels=typelevels)
ggdft$site <- factor(output$meta_data$Anno[match(ggdft$sample_id,output$meta_data$sample_id)], levels=sitelevels)
ggdft$met <- factor(output$meta_data$Met[match(ggdft$sample_id,output$meta_data$sample_id)])
ggdft$paired <- factor(output$meta_data$Paired[match(ggdft$sample_id,output$meta_data$sample_id)])
ggdft$TMA <- output$meta_data$TMA[match(ggdft$sample_id,output$meta_data$sample_id)]

totalcounts<-ggdft %>% group_by(cluster, SPC, tumor, type, site, met, paired) %>% summarize_at(vars(counts),funs(sum))

write.csv(totalcounts,"Totalcounts.csv")

totalcounts_celltype <- ggdft %>% group_by(cluster, met) %>% summarize_at(vars(counts),funs(sum))

write.csv(totalcounts_celltype, "Totalcounts_celltype_met.csv")

totalcounts_tumortype <- ggdft %>% group_by(cluster, type) %>% summarize_at(vars(counts),funs(sum))

write.csv(totalcounts_tumortype, "Totalcounts_tumortype_met.csv")

#percentage of each respective total

totalcounts_met <- totalcounts_celltype[totalcounts_celltype$met=="1",]
totalcounts_primary <- totalcounts_celltype[totalcounts_celltype$met=="0",]
totalmet<-sum(totalcounts_met$counts)
totalprimary<-sum(totalcounts_primary$counts)
pct_met <- totalcounts_met$counts/totalmet*100
names(pct_met)<- totalcounts_met$cluster
pct_primary <- totalcounts_primary$counts/totalprimary*100
names(pct_primary)<- totalcounts_primary$cluster

ggdf_pctm <- melt(pct_met);ggdf_pctm$met<-"met"
ggdf_pctp <- melt(pct_primary);ggdf_pctp$met<-"primary"
ggdf_pct<-rbind(ggdf_pctp,ggdf_pctm)
ggdf_pct$cluster<-c(rownames(ggdf_pctp),rownames(ggdf_pctm))
rownames(ggdf_pct)<-1:nrow(ggdf_pct)
ggdf_pct$met<-factor(ggdf_pct$met, levels=c("primary","met"))

bp <- ggplot(ggdf_pct, aes(x = met, y = value, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))
pdf('Abundance_stackedbar_met.pdf', width=3, height=4); bp; dev.off()


totalcounts_LUM <- totalcounts_tumortype[totalcounts_tumortype$type=="LUM",]
totalcounts_TNC <- totalcounts_tumortype[totalcounts_tumortype$type=="TNC",]
totalLUM<-sum(totalcounts_LUM$counts)
totalTNC<-sum(totalcounts_TNC$counts)
pct_LUM <- totalcounts_LUM$counts/totalLUM*100
names(pct_LUM)<- totalcounts_LUM$cluster
pct_TNC <- totalcounts_TNC$counts/totalTNC*100
names(pct_TNC)<- totalcounts_TNC$cluster

saveRDS(pct_LUM, 'pct_LUM.rds')
saveRDS(pct_TNC, 'pct_TNC.rds')

ggdf_pctLUM <- melt(pct_LUM);ggdf_pctLUM$type<-"LUM"
ggdf_pctTNC <- melt(pct_TNC);ggdf_pctTNC$type<-"TNC"
ggdf_pcttype<-rbind(ggdf_pctLUM,ggdf_pctTNC)
ggdf_pcttype$cluster<-c(rownames(ggdf_pctLUM),rownames(ggdf_pctTNC))
rownames(ggdf_pcttype)<-1:nrow(ggdf_pcttype)
ggdf_pcttype$type<-factor(ggdf_pcttype$type, levels=c("LUM","TNC"))

bp <- ggplot(ggdf_pcttype, aes(x = type, y = value, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=2))
pdf('Abundance_stackedbar_type.pdf', width=3, height=4); bp; dev.off()

#identify which cell types in either of the data subsets are very rare - have less than 0.01%
#also include NA cluster = ctype14
exclude_LUM<-which(pct_LUM<0.01)
exclude_LUM<-legendctype$V1[match(names(exclude_LUM), legendctype$allcelltypes)]
exclude_LUM<-union(exclude_LUM, "ctype31")
exclude_TNC<-which(pct_TNC<0.01)
exclude_TNC<-legendctype$V1[match(names(exclude_TNC), legendctype$allcelltypes)]
exclude_TNC<-union(exclude_TNC, "ctype31")


#get out expression levels, X, and Y coords

expr <- fsApply(output$fcs1, exprs) #create expression matrix

expr0<-data.frame(expr[,c(union(output$subtype_markers,output$functional_markers),"CellId","X_coord","Y_coord")],
                  cluster=output$cell_clustering1m,
                  sample_id=output$sample_ids)
expr0$SPC <- factor(output$meta_data$Case[match(expr0$sample_id,output$meta_data$sample_id)])
expr0$tumor <- factor(output$meta_data$Tumor[match(expr0$sample_id,output$meta_data$sample_id)], levels=tumorlevels)
expr0$type <- factor(output$meta_data$Phenotype[match(expr0$sample_id,output$meta_data$sample_id)], levels=typelevels)
expr0$site <- factor(output$meta_data$Anno[match(expr0$sample_id,output$meta_data$sample_id)], levels=sitelevels)
expr0$met <- factor(output$meta_data$Met[match(expr0$sample_id,output$meta_data$sample_id)])
expr0$paired <- factor(output$meta_data$Paired[match(expr0$sample_id,output$meta_data$sample_id)])
expr0$TMA <- output$meta_data$TMA[match(expr0$sample_id,output$meta_data$sample_id)]

expr0<-expr0[expr0$type %in% c("LUM","TNC"),]

expr0<-expr0[expr0$cluster!="NA",]

expr0_LUM<-expr0[expr0$type=="LUM",]
expr0_TNC<-expr0[expr0$type=="TNC",]

###CREATE DISTANCE MATRICES FOR PROGRESSION CRITERIA IN LUM or TNC ONLY

LUM<-unique(expr0[expr0$type=="LUM",]$sample_id)
TNC<-unique(expr0[expr0$type=="TNC",]$sample_id)


##LUM

expr_LUM<-c()

for(k in 1:length(LUM)){
  expr_k<-expr0[expr0$sample_id==LUM[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy <- matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy) <- legendctype$allcelltypes
  dummy <- as.data.frame(dummy)
  
  #since cell type and expression marker are named the same for CK8 and CK14, we need to rename them before running distance calculations on merged dfs 
  colnames(expr_k)[which(colnames(expr_k)=="CK8")]<-"CK8ex" #replace CK8 expression colname before merging to avoid auto renaming the column for distances to CK8
  colnames(expr_k)[which(colnames(expr_k)=="CK14")]<-"CK14ex" #replace CK14 expression colname before merging to avoid auto renaming the column for distances to CK14
  
  expr_k <- data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  
  expr_LUM <- rbind(expr_LUM, expr_k)
}
  
colnames(expr_LUM)<-c(colnames(expr0_LUM),legendctype$V1)
expr_LUM_m<-as.matrix(expr_LUM[,colnames(expr_LUM)[str_detect(colnames(expr_LUM),"ctype")]])
expr_LUM$ctype_no<-paste0("ctype",match(expr_LUM$cluster,legendctype$allcelltypes))
rownames(expr_LUM_m)<-expr_LUM$ctype_no
expr_LUM_m[is.infinite(expr_LUM_m)]<-NA

expr_LUM<-expr_LUM_m
expr_LUM[expr_LUM == 0] <- NA

saveRDS(expr_LUM,'backup_dist_LUM_rev.rds')

#load previously saved distance matrix

expr_LUM<- readRDS('backup_dist_LUM_rev.rds')

#create expression + distance data frames

expr_LUMcombined <- cbind(expr0_LUM,expr_LUM)

expr_LUM_CK<-expr_LUMcombined[expr_LUMcombined$cluster %in% c("CK8","CK8_E","CK8_V","CK14","CK14_E","CK14_EV","CK8_14_E","CK8_14_EV","CKlo","CKlo_EV","CKother"),]

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_LUMrm0<-expr_LUM[,colnames(expr_LUM) %nin% colnames(expr_LUM)[colSums(expr_LUM, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types

expr_LUMrm<-expr_LUMrm0[rownames(expr_LUMrm0) %nin% exclude_LUM, colnames(expr_LUMrm0) %nin% exclude_LUM]

#create new matrices with primary and mets only
expr_LUMcombined$ctype<-rownames(expr_LUM) #add a column for "ctypes" to add as rownames later

expr_LUMpri<-expr_LUMcombined[expr_LUMcombined$met==0,] #filter out mets
expr_LUMpri_mat<-as.matrix(expr_LUMpri[,paste0("ctype",1:31)]) #select the ctype columns only
rownames(expr_LUMpri_mat)<-expr_LUMpri$ctype #add back the ctypes as rownames
expr_LUMpri_mat0<-expr_LUMpri_mat[,colnames(expr_LUMpri_mat) %nin% colnames(expr_LUMpri_mat)[colSums(expr_LUMpri_mat, na.rm = T)==0]]
expr_LUMpri_mat<-expr_LUMpri_mat0[rownames(expr_LUMpri_mat0) %nin% exclude_LUM, colnames(expr_LUMpri_mat0) %nin% exclude_LUM]

expr_LUMmet<-expr_LUMcombined[expr_LUMcombined$met==1,] #filter out nonmets
expr_LUMmet_mat<-as.matrix(expr_LUMmet[,paste0("ctype",1:31)]) #select the ctype columns only
rownames(expr_LUMmet_mat)<-expr_LUMmet$ctype #add back the ctypes as rownames
expr_LUMmet_mat0<-expr_LUMmet_mat[,colnames(expr_LUMmet_mat) %nin% colnames(expr_LUMmet_mat)[colSums(expr_LUMmet_mat, na.rm = T)==0]]
expr_LUMmet_mat<-expr_LUMmet_mat0[rownames(expr_LUMmet_mat0) %nin% exclude_LUM, colnames(expr_LUMmet_mat0) %nin% exclude_LUM]

saveRDS(expr_LUMrm, 'expr_LUMrm.rds')
saveRDS(expr_LUMpri_mat, 'expr_LUMpri_mat.rds')
saveRDS(expr_LUMmet_mat, 'expr_LUMmet_mat.rds')


##TNC

expr_TNC<-c()

for(k in 1:length(TNC)){
  expr_k<-expr0[expr0$sample_id==TNC[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy <- matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy) <- legendctype$allcelltypes
  dummy <- as.data.frame(dummy)
  
  #since cell type and expression marker are named the same for CK8 and CK14, we need to rename them before running distance calculations on merged dfs 
  colnames(expr_k)[which(colnames(expr_k)=="CK8")]<-"CK8ex" #replace CK8 expression colname before merging to avoid auto renaming the column for distances to CK8
  colnames(expr_k)[which(colnames(expr_k)=="CK14")]<-"CK14ex" #replace CK14 expression colname before merging to avoid auto renaming the column for distances to CK14
  
  expr_k <- data.frame(expr_k,dummy)
  
  X <- ppp(x=expr_k$X_coord,y=expr_k$Y_coord,marks = as.factor(expr_k$cluster), window = owin(xrange = range(expr_k$X_coord), yrange = range(expr_k$Y_coord)),checkdup = F)
  distByCelltype <- nndist(X, by=marks(X))
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for (i in 1:(ncol(distByCelltype))){
    d <- which(colnames(expr_k) == (colnames(distByCelltype)[i]))
    expr_k[d] <- distByCelltype[,i]
  }
  expr_TNC <- rbind(expr_TNC, expr_k)
}

colnames(expr_TNC)<-c(colnames(expr0_TNC),legendctype$V1)
expr_TNC_m<-as.matrix(expr_TNC[,colnames(expr_TNC)[str_detect(colnames(expr_TNC),"ctype")]])
expr_TNC$ctype_no<-paste0("ctype",match(expr_TNC$cluster,legendctype$allcelltypes))
rownames(expr_TNC_m)<-expr_TNC$ctype_no
expr_TNC_m[is.infinite(expr_TNC_m)]<-NA

expr_TNC<-expr_TNC_m
expr_TNC[expr_TNC == 0] <- NA

saveRDS(expr_TNC,'backup_dist_TNC_rev.rds')

#load previously saved distance matrix

expr_TNC<-readRDS('backup_dist_TNC_rev.rds')

#create expression + distance data frames

expr_TNCcombined <- cbind(expr0_TNC,expr_TNC)

expr_TNC_CK<-expr_TNCcombined[expr_TNCcombined$cluster %in% c("CK8","CK8_E","CK8_V","CK14","CK14_E","CK14_EV","CK8_14_E","CK8_14_EV","CKlo","CKlo_EV","CKother"),]

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_TNCrm0<-expr_TNC[,colnames(expr_TNC) %nin% colnames(expr_TNC)[colSums(expr_TNC, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.5%) cell types

expr_TNCrm<-expr_TNCrm0[rownames(expr_TNCrm0) %nin% exclude_TNC, colnames(expr_TNCrm0) %nin% exclude_TNC]

#create new matrices with primary and mets only
expr_TNCcombined$ctype<-rownames(expr_TNC) #add a column for "ctypes" to add as rownames later

expr_TNCpri<-expr_TNCcombined[expr_TNCcombined$met==0,] #filter out mets
expr_TNCpri_mat<-as.matrix(expr_TNCpri[,paste0("ctype",1:31)]) #select the ctype columns only
rownames(expr_TNCpri_mat)<-expr_TNCpri$ctype #add back the ctypes as rownames
expr_TNCpri_mat0<-expr_TNCpri_mat[,colnames(expr_TNCpri_mat) %nin% colnames(expr_TNCpri_mat)[colSums(expr_TNCpri_mat, na.rm = T)==0]]
expr_TNCpri_mat<-expr_TNCpri_mat0[rownames(expr_TNCpri_mat0) %nin% exclude_TNC, colnames(expr_TNCpri_mat0) %nin% exclude_TNC]

expr_TNCmet<-expr_TNCcombined[expr_TNCcombined$met==1,] #filter out nonmets
expr_TNCmet_mat<-as.matrix(expr_TNCmet[,paste0("ctype",1:31)]) #select the ctype columns only
rownames(expr_TNCmet_mat)<-expr_TNCmet$ctype #add back the ctypes as rownames
expr_TNCmet_mat0<-expr_TNCmet_mat[,colnames(expr_TNCmet_mat) %nin% colnames(expr_TNCmet_mat)[colSums(expr_TNCmet_mat, na.rm = T)==0]]
expr_TNCmet_mat<-expr_TNCmet_mat0[rownames(expr_TNCmet_mat0) %nin% exclude_TNC, colnames(expr_TNCmet_mat0) %nin% exclude_TNC]

saveRDS(expr_TNCrm, 'expr_TNCrm.rds')
saveRDS(expr_TNCpri_mat, 'expr_TNCpri_mat.rds')
saveRDS(expr_TNCmet_mat, 'expr_TNCmet_mat.rds')

###3D plots for macrophage proximity based on ECAD vs. VIM and CK14 vs. CK8

##ECAD vs. VIM
akimainterp <- function(exprdf,ctype){
  exprdf_f<-exprdf[!is.na(exprdf[[ctype]]),]
  akima::interp(x=exprdf_f$ECAD,exprdf_f$VIM,1/exprdf_f[[ctype]],duplicate="mean")
}

##CK14 vs. CK8
akimainterp2 <- function(exprdf,ctype){
  exprdf_f<-exprdf[!is.na(exprdf[[ctype]]),]
  akima::interp(x=exprdf_f$CK8,exprdf_f$CK14,1/exprdf_f[[ctype]],duplicate="mean")
}

##CK14 vs. VIM
akimainterp3 <- function(exprdf,ctype){
  exprdf_f<-exprdf[!is.na(exprdf[[ctype]]),]
  akima::interp(x=exprdf_f$VIM,exprdf_f$CK14,1/exprdf_f[[ctype]],duplicate="mean")
}

##for perspective 3D plots of ECAD vs. VIM

#TNC

int_ctype11<-akimainterp(exprdf=expr_TNC_CK,ctype="ctype11")
int_ctype12<-akimainterp(exprdf=expr_TNC_CK,ctype="ctype12")
int_ctype13<-akimainterp(exprdf=expr_TNC_CK,ctype="ctype13")
int_ctype14<-akimainterp(exprdf=expr_TNC_CK,ctype="ctype14")

pdf("Distance_3DVE_TNC.pdf",width=5, height=5)

persp3D(int_ctype11$x,int_ctype11$y, int_ctype11$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="TNC: CK+ to Mac_I")
persp3D(int_ctype12$x,int_ctype12$y, int_ctype12$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="TNC: CK+ to Mac_II")
persp3D(int_ctype13$x,int_ctype13$y, int_ctype13$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="TNC: CK+ to Mac_III")
persp3D(int_ctype14$x,int_ctype14$y, int_ctype14$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="TNC: CK+ to Mac_IV")

scatter3D(x=expr_TNC_CK$ECAD,
          y=expr_TNC_CK$VIM, 
          z=1/expr_TNC_CK$ctype11,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="TNC: CK+ to Mac_I")
scatter3D(x=expr_TNC_CK$ECAD, 
          y=expr_TNC_CK$VIM, 
          z=1/expr_TNC_CK$ctype12,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="TNC: CK+ to Mac_II")
scatter3D(x=expr_TNC_CK$ECAD, 
          y=expr_TNC_CK$VIM, 
          z=1/expr_TNC_CK$ctype13,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="TNC: CK+ to Mac_III")
scatter3D(x=expr_TNC_CK$ECAD, 
          y=expr_TNC_CK$VIM, 
          z=1/expr_TNC_CK$ctype14,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="TNC: CK+ to Mac_IV")
dev.off()

#LUM

int_ctype11<-akimainterp(exprdf=expr_LUM_CK,ctype="ctype11")
int_ctype12<-akimainterp(exprdf=expr_LUM_CK,ctype="ctype12")
int_ctype13<-akimainterp(exprdf=expr_LUM_CK,ctype="ctype13")
int_ctype14<-akimainterp(exprdf=expr_LUM_CK,ctype="ctype14")

pdf("Distance_3DVE_LUM.pdf",width=5, height=5)

persp3D(int_ctype11$x,int_ctype11$y, int_ctype11$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_I")
persp3D(int_ctype12$x,int_ctype12$y, int_ctype12$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_II")
persp3D(int_ctype13$x,int_ctype13$y, int_ctype13$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_III")
persp3D(int_ctype14$x,int_ctype14$y, int_ctype14$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="ECAD",
        #ylim=c(0,5),
        ylab="VIM",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_IV")

scatter3D(x=expr_LUM_CK$ECAD,
          y=expr_LUM_CK$VIM, 
          z=1/expr_LUM_CK$ctype11,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_I")
scatter3D(x=expr_LUM_CK$ECAD, 
          y=expr_LUM_CK$VIM, 
          z=1/expr_LUM_CK$ctype12,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_II")
scatter3D(x=expr_LUM_CK$ECAD, 
          y=expr_LUM_CK$VIM, 
          z=1/expr_LUM_CK$ctype13,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_III")
scatter3D(x=expr_LUM_CK$ECAD, 
          y=expr_LUM_CK$VIM, 
          z=1/expr_LUM_CK$ctype14,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="ECAD",
          #ylim=c(0,5),
          ylab="VIM",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_IV")
dev.off()

##for perspective 3D plots of CK8 vs. CK14 for LUMINAL only

int_ctype11b<-akimainterp2(exprdf=expr_LUM_CK,ctype="ctype11")
int_ctype12b<-akimainterp2(exprdf=expr_LUM_CK,ctype="ctype12")
int_ctype13b<-akimainterp2(exprdf=expr_LUM_CK,ctype="ctype13")
int_ctype14b<-akimainterp2(exprdf=expr_LUM_CK,ctype="ctype14")

pdf("Distance_3DVE_LUM_CK8_CK14.pdf",width=5, height=5)

persp3D(int_ctype11b$x,int_ctype11b$y, int_ctype11b$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="CK8",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_I")
persp3D(int_ctype12b$x,int_ctype12b$y, int_ctype12b$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="CK8",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_II")
persp3D(int_ctype13b$x,int_ctype13b$y, int_ctype13b$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="CK8",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_III")
persp3D(int_ctype14b$x,int_ctype14b$y, int_ctype14b$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="CK8",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_IV")

scatter3D(x=expr_LUM_CK$CK8,
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype11,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="CK8",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_I")
scatter3D(x=expr_LUM_CK$CK8, 
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype12,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="CK8",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_II")
scatter3D(x=expr_LUM_CK$CK8, 
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype13,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="CK8",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_III")
scatter3D(x=expr_LUM_CK$CK8, 
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype14,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="CK8",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_IV")
dev.off()


##for perspective 3D plots of VIM vs. CK14 for LUMINAL only

int_ctype11c<-akimainterp3(exprdf=expr_LUM_CK,ctype="ctype11")
int_ctype12c<-akimainterp3(exprdf=expr_LUM_CK,ctype="ctype12")
int_ctype13c<-akimainterp3(exprdf=expr_LUM_CK,ctype="ctype13")
int_ctype14c<-akimainterp3(exprdf=expr_LUM_CK,ctype="ctype14")

pdf("Distance_3DVE_LUM_VIM_CK14.pdf",width=5, height=5)

persp3D(int_ctype11c$x,int_ctype11c$y, int_ctype11c$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="VIM",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_I")
persp3D(int_ctype12c$x,int_ctype12c$y, int_ctype12c$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="VIM",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_II")
persp3D(int_ctype13c$x,int_ctype13c$y, int_ctype13c$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="VIM",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_III")
persp3D(int_ctype14c$x,int_ctype14c$y, int_ctype14c$z,
        phi=40, theta=60,
        #contour = list(col="black", side=c(0,.05)),
        bty="b",
        #xlim=c(0,2.5),
        xlab="VIM",
        #ylim=c(0,5),
        ylab="CK14",
        zlim=c(0,0.16),
        zlab="proximity", #(1/avg shortest dist)
        col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
        clim=c(0,0.16),
        clab="Proximity",
        main="LUM: CK+ to Mac_IV")

scatter3D(x=expr_LUM_CK$VIM,
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype11,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="VIM",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_I")
scatter3D(x=expr_LUM_CK$VIM, 
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype12,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="VIM",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_II")
scatter3D(x=expr_LUM_CK$VIM, 
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype13,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="VIM",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_III")
scatter3D(x=expr_LUM_CK$VIM, 
          y=expr_LUM_CK$CK14, 
          z=1/expr_LUM_CK$ctype14,
          phi=40, theta=60,
          #contour = list(col="black", side=c(0,.05)),
          bty="b",
          #xlim=c(0,2.5),
          xlab="VIM",
          #ylim=c(0,5),
          ylab="CK14",
          zlim=c(0,0.16),
          zlab="proximity", #(1/avg shortest dist)
          col=coolwarm(100)[c(seq(1,57,by=3),60:100)],
          clim=c(0,0.28),
          clab="Proximity",
          main="LUM: CK+ to Mac_IV")
dev.off()


####PAIRED COMPARISON####

filterunpaired<-function(data){
  data<-data[data$site2 %in% c("PBC","BRAIN","LUNG","LN","LIVER"),] #only keep sites that will be plotted
  filter<-table(data$SPC,data$site2) #enumerate how many values there are for each SPC and each site2
  filter<-filter[rowSums(filter)!=0,colSums(filter)!=0] #get rid of rows and cols with only 0
  filter[rownames(filter)[filter[,"PBC"]!=0],] #get rid of rows without "PBC"
  filter<-filter[rowSums(filter)>=2,] #get rid of those that have less than 2 values
  keep<-rownames(filter)
  return(data[data$SPC %in% keep,])
}

expr_TNC_CK$site2 <- expr_TNC_CK$site
expr_TNC_CK$site2[str_detect(expr_TNC_CK$site2,"PBC")] <- "PBC" 

subdf11<-expr_TNC_CK[!is.na(expr_TNC_CK[["ctype11"]]),]
subdf11_close<-subdf11[subdf11$ctype11<50,]
subdf11_far<-subdf11[!subdf11$ctype11<50,]

subdf12<-expr_TNC_CK[!is.na(expr_TNC_CK[["ctype12"]]),]
subdf12_close<-subdf12[subdf12$ctype12<50,]
subdf12_far<-subdf12[!subdf12$ctype12<50,]

subdf13<-expr_TNC_CK[!is.na(expr_TNC_CK[["ctype13"]]),]
subdf13_close<-subdf13[subdf13$ctype13<50,]
subdf13_far<-subdf13[!subdf13$ctype13<50,]

subdf14<-expr_TNC_CK[!is.na(expr_TNC_CK[["ctype14"]]),]
subdf14_close<-subdf14[subdf14$ctype14<50,]
subdf14_far<-subdf14[!subdf14$ctype14<50,]

df11c <- subdf11_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df11c$mac <- "Mac_I" ; df11c$prox <- "close" ; df11c<-filterunpaired(df11c)
df11f <- subdf11_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df11f$mac <- "Mac_I" ; df11f$prox <- "far"; df11f<-filterunpaired(df11f)

df12c <- subdf12_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df12c$mac <- "Mac_II" ; df12c$prox <- "close" ; df12c<-filterunpaired(df12c)
df12f <- subdf12_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df12f$mac <- "Mac_II" ; df12f$prox <- "far" ; df12f<-filterunpaired(df12f)

df13c <- subdf13_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df13c$mac <- "Mac_III" ; df13c$prox <- "close" ; df13c<-filterunpaired(df13c)
df13f <- subdf13_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df13f$mac <- "Mac_III" ; df13f$prox <- "far" ; df13f<-filterunpaired(df13f)

df14c <- subdf14_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df14c$mac <- "Mac_IV" ; df14c$prox <- "close" ; df14c<-filterunpaired(df14c)
df14f <- subdf14_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df14f$mac <- "Mac_IV" ; df14f$prox <- "far" ; df14f<-filterunpaired(df14f)

dfallmac <- rbind(df11c,df11f,df12c,df12f,df13c,df13f,df14c,df14f)
dfallmac1<-dfallmac
dfallmac1$type <- "TNC"

expr_LUM_CK$site2 <- expr_LUM_CK$site
expr_LUM_CK$site2[str_detect(expr_LUM_CK$site2,"PBC")] <- "PBC" 

subdf11<-expr_LUM_CK[!is.na(expr_LUM_CK[["ctype11"]]),]
subdf11_close<-subdf11[subdf11$ctype11<50,]
subdf11_far<-subdf11[!subdf11$ctype11<50,]

subdf12<-expr_LUM_CK[!is.na(expr_LUM_CK[["ctype12"]]),]
subdf12_close<-subdf12[subdf12$ctype12<50,]
subdf12_far<-subdf12[!subdf12$ctype12<50,]

subdf13<-expr_LUM_CK[!is.na(expr_LUM_CK[["ctype13"]]),]
subdf13_close<-subdf13[subdf13$ctype13<50,]
subdf13_far<-subdf13[!subdf13$ctype13<50,]

subdf14<-expr_LUM_CK[!is.na(expr_LUM_CK[["ctype14"]]),]
subdf14_close<-subdf14[subdf14$ctype14<50,]
subdf14_far<-subdf14[!subdf14$ctype14<50,]

df11c <- subdf11_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df11c$mac <- "Mac_I" ; df11c$prox <- "close" ; df11c<-filterunpaired(df11c)
df11f <- subdf11_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df11f$mac <- "Mac_I" ; df11f$prox <- "far"; df11f<-filterunpaired(df11f)

df12c <- subdf12_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df12c$mac <- "Mac_II" ; df12c$prox <- "close" ; df12c<-filterunpaired(df12c)
df12f <- subdf12_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df12f$mac <- "Mac_II" ; df12f$prox <- "far" ; df12f<-filterunpaired(df12f)

df13c <- subdf13_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df13c$mac <- "Mac_III" ; df13c$prox <- "close" ; df13c<-filterunpaired(df13c)
df13f <- subdf13_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df13f$mac <- "Mac_III" ; df13f$prox <- "far" ; df13f<-filterunpaired(df13f)

df14c <- subdf14_close %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df14c$mac <- "Mac_IV" ; df14c$prox <- "close" ; df14c<-filterunpaired(df14c)
df14f <- subdf14_far %>% group_by(SPC,site2) %>% summarise_at(., vars(VIM, KI67, HLADR), mean) ; df14f$mac <- "Mac_IV" ; df14f$prox <- "far" ; df14f<-filterunpaired(df14f)

dfallmac <- rbind(df11c,df11f,df12c,df12f,df13c,df13f,df14c,df14f)
dfallmac2<-dfallmac
dfallmac2$type<-"LUM"

dfallmac <- rbind(dfallmac1,dfallmac2)




ggplotmac1k <- ggplot(dfallmac[dfallmac$site2 %in% c("PBC","BRAIN","LUNG","LN","LIVER"),], aes(x=prox, y=KI67))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(site2 ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac1v <- ggplot(dfallmac[dfallmac$site2 %in% c("PBC","BRAIN","LUNG","LN","LIVER"),], aes(x=prox, y=VIM))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(site2 ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac1h <- ggplot(dfallmac[dfallmac$site2 %in% c("PBC","BRAIN","LUNG","LN","LIVER"),], aes(x=prox, y=HLADR))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(site2 ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

pdf("Expression_byProximity+.pdf", width=5, height=6.5)
ggplotmac1k
ggplotmac1v
ggplotmac1h
dev.off()

write.csv(dfallmac[dfallmac$site2 %in% c("PBC","BRAIN","LUNG","LN","LIVER"),], "Results_ExpressionbyProximity.csv")


#filter out unpaired samples for within each organ comparison
pvm<-dfallmac[dfallmac$site2 %in% c("PBC","BRAIN"),]
pvm$site2<-factor(pvm$site2)
pvm$SPC<-factor(pvm$SPC)
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,1]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"BRAIN"]]
pvm1<-pvm[pvm$mac=="Mac_I",]
pvm1<-pvm1[pvm1$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,2]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"BRAIN"]]
pvm2<-pvm[pvm$mac=="Mac_II",]
pvm2<-pvm2[pvm2$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,3]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"BRAIN"]]
pvm3<-pvm[pvm$mac=="Mac_III",]
pvm3<-pvm3[pvm3$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,4]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"BRAIN"]]
pvm4<-pvm[pvm$mac=="Mac_IV",]
pvm4<-pvm4[pvm4$SPC %in% keep,]
pvm<-rbind(pvm1,pvm2,pvm3,pvm4)

ggplotmac2k <- ggplot(pvm, aes(x=site2, y=KI67))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac2v <- ggplot(pvm, aes(x=site2, y=VIM))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac2h <- ggplot(pvm, aes(x=site2, y=HLADR))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

pdf("Expression_byProximity_PBCvBRAIN+.pdf", width=5, height=3.3)
ggplotmac2k
ggplotmac2v
ggplotmac2h
dev.off()


pvm<-dfallmac[dfallmac$site2 %in% c("PBC","LUNG"),]
pvm$site2<-factor(pvm$site2)
pvm$SPC<-factor(pvm$SPC)
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,1]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LUNG"]]
pvm1<-pvm[pvm$mac=="Mac_I",]
pvm1<-pvm1[pvm1$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,2]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LUNG"]]
pvm2<-pvm[pvm$mac=="Mac_II",]
pvm2<-pvm2[pvm2$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,3]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LUNG"]]
pvm3<-pvm[pvm$mac=="Mac_III",]
pvm3<-pvm3[pvm3$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,4]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LUNG"]]
pvm4<-pvm[pvm$mac=="Mac_IV",]
pvm4<-pvm4[pvm4$SPC %in% keep,]
pvm<-rbind(pvm1,pvm2,pvm3,pvm4)

ggplotmac3k <- ggplot(pvm, aes(x=site2, y=KI67))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac3v <- ggplot(pvm, aes(x=site2, y=VIM))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac3h <- ggplot(pvm, aes(x=site2, y=HLADR))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

pdf("Expression_byProximity_PBCvLUNG.pdf", width=5, height=3.3)
ggplotmac3k
ggplotmac3v
ggplotmac3h
dev.off()


pvm<-dfallmac[dfallmac$site2 %in% c("PBC","LN"),]
pvm$site2<-factor(pvm$site2)
pvm$SPC<-factor(pvm$SPC)
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,1]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LN"]]
pvm1<-pvm[pvm$mac=="Mac_I",]
pvm1<-pvm1[pvm1$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,2]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LN"]]
pvm2<-pvm[pvm$mac=="Mac_II",]
pvm2<-pvm2[pvm2$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,3]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LN"]]
pvm3<-pvm[pvm$mac=="Mac_III",]
pvm3<-pvm3[pvm3$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,4]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LN"]]
pvm4<-pvm[pvm$mac=="Mac_IV",]
pvm4<-pvm4[pvm4$SPC %in% keep,]
pvm<-rbind(pvm1,pvm2,pvm3,pvm4)

ggplotmac4k <- ggplot(pvm, aes(x=site2, y=KI67))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac4v <- ggplot(pvm, aes(x=site2, y=VIM))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac4h <- ggplot(pvm, aes(x=site2, y=HLADR))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

pdf("Expression_byProximity_PBCvLN.pdf", width=5, height=3.3)
ggplotmac4k
ggplotmac4v
ggplotmac4h
dev.off()


pvm<-dfallmac[dfallmac$site2 %in% c("PBC","LIVER"),]
pvm$site2<-factor(pvm$site2)
pvm$SPC<-factor(pvm$SPC)
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,1]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LIVER"]]
pvm1<-pvm[pvm$mac=="Mac_I",]
pvm1<-pvm1[pvm1$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,2]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LIVER"]]
pvm2<-pvm[pvm$mac=="Mac_II",]
pvm2<-pvm2[pvm2$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,3]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LIVER"]]
pvm3<-pvm[pvm$mac=="Mac_III",]
pvm3<-pvm3[pvm3$SPC %in% keep,]
pvm_table<-table(pvm$SPC,pvm$site2,pvm$mac)
pvm_table<-pvm_table[,,4]
pvm_table<-pvm_table[rowSums(pvm_table)!=0,colSums(pvm_table)!=0]
keep<-rownames(pvm_table)[pvm_table[,"PBC"]==pvm_table[,"LIVER"]]
pvm4<-pvm[pvm$mac=="Mac_IV",]
pvm4<-pvm4[pvm4$SPC %in% keep,]
pvm<-rbind(pvm1,pvm2,pvm3,pvm4)

ggplotmac5k <- ggplot(pvm, aes(x=site2, y=KI67))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac5v <- ggplot(pvm, aes(x=site2, y=VIM))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

ggplotmac5h <- ggplot(pvm, aes(x=site2, y=HLADR))+
  geom_line(aes(group=SPC, color=type))+
  geom_point(aes(group=SPC, color=type))+
  facet_grid(prox ~ mac, scales = "free")+
  scale_y_continuous(expand = expansion(mult = c(0.15,0.33)))+
  stat_compare_means(method="t.test",
                     paired=T, 
                     label = "p.signif", 
                     label.x = 1.4, 
                     label.y.npc="top",
                     hide.ns = T,
                     size=5)+
  theme(
    panel.border = element_rect(linewidth=0.25, color="black", fill=NA),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_line(color="black",linewidth=0.25),
    axis.text = element_text(color="black")
  )

pdf("Expression_byProximity_PBCvLIVER.pdf", width=5, height=3.3)
ggplotmac5k
ggplotmac5v
ggplotmac5h
dev.off()


#distance data needed
expr_LUMpri_CK <- expr_LUMpri[expr_LUMpri$cluster %in% c("CK8","CK8_E","CK8_V","CK14","CK14_E","CK14_EV","CK8_14_E","CK8_14_EV","CKlo","CKlo_EV","CKother"),]
expr_LUMmet_CK <- expr_LUMmet[expr_LUMmet$cluster %in% c("CK8","CK8_E","CK8_V","CK14","CK14_E","CK14_EV","CK8_14_E","CK8_14_EV","CKlo","CKlo_EV","CKother"),]
expr_TNCpri_CK <- expr_TNCpri[expr_TNCpri$cluster %in% c("CK8","CK8_E","CK8_V","CK14","CK14_E","CK14_EV","CK8_14_E","CK8_14_EV","CKlo","CKlo_EV","CKother"),]
expr_TNCmet_CK <- expr_TNCmet[expr_TNCmet$cluster %in% c("CK8","CK8_E","CK8_V","CK14","CK14_E","CK14_EV","CK8_14_E","CK8_14_EV","CKlo","CKlo_EV","CKother"),]


#input the distance df you want
subsetdf_distCK <- expr_TNCmet_CK
#this will name the heatmap and the file accordingly
name <- "TNC MET"



##Run the code to obtain a distance df and generate heatmaps

meandistfromCK <- subsetdf_distCK %>% group_by(cluster) %>% summarise(res=mean(ctype1,na.rm=T))
for (i in 2:18){
  ctypei <- legendctype$V1[i]
  mnpercluster <- subsetdf_distCK %>% group_by(cluster) %>% summarise(res=mean(get(ctypei),na.rm=T))
  meandistfromCK <- cbind(meandistfromCK, res=mnpercluster$res)
}



colnames(meandistfromCK)[2:ncol(meandistfromCK)] <- legendctype$allcelltypes[1:18]
rownames(meandistfromCK) <- meandistfromCK$cluster
meandistfromCK <- meandistfromCK[,2:ncol(meandistfromCK)]
meandistfromCK <- as.matrix(meandistfromCK)

col_fun = colorRamp2(c(0, 400, 800), c("red", "white", "blue"))

pdf(paste0("Distance_MeanCKtoTME_",name,".pdf"), width=8, height=8)
Heatmap(meandistfromCK,
        name = name,
        col = col_fun,
        border = T,
        cluster_rows = F,
        cluster_columns = F,
        column_names_side = "top",
        width = ncol(meandistfromCK)*unit(5, "mm"), 
        height = nrow(meandistfromCK)*unit(5, "mm")
)
dev.off()


