setwd("G:/Qiaolab/liudandan/human/ALL_data")
TPM <- read.table("qc_tpm(1127).txt",as.is = TRUE)
#TPM.log <- log2(TPM/10+1)
Stage.class <- sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 1 )
S2.tpm <- TPM[ ,which(Stage.class %in% c("S2"))]
S2.tpm.log <- log2(S2.tpm/10+1)
Embryo.class <- sapply( strsplit(as.character(colnames(S2.tpm)), "_"), "[[", 3 )

library(Seurat)
library(dplyr)
nbt <- CreateSeuratObject(S2.tpm, project = "S2", min.cells = 3, min.features = 200)
nbt
nbt[["percent.mt"]] <- PercentageFeatureSet(nbt, pattern = "^MT-")
VlnPlot(nbt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
nbt <- NormalizeData(nbt, normalization.method = "LogNormalize", scale.factor = 10000)
nbt <- NormalizeData(nbt)
nbt <- FindVariableFeatures(nbt, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
VariableFeaturePlot(nbt)
S2_vargene <- data.frame(name=VariableFeatures(object = nbt))
write.table(S2_vargene,"human_S2_vargene.txt",quote = F,col.names = F,row.names = F)
nbt <- ScaleData(nbt)
nbt <- RunPCA(nbt, features = VariableFeatures(object = nbt))
DimPlot(nbt, reduction = "pca")
DimHeatmap(nbt, dims = 1:15, cells = 500, balanced = TRUE)
nbt <- JackStraw(nbt, num.replicate = 100)
nbt <- ScoreJackStraw(nbt, dims = 1:20)
JackStrawPlot(nbt, dims = 1:20)
ElbowPlot(nbt)
nbt <- FindNeighbors(nbt, dims = 1:7)
nbt <- FindClusters(nbt, resolution = 0.2)
nbt <- RunUMAP(nbt, dims = 1:7)
DimPlot(nbt, reduction = "umap",label = T)


nbt.out=data.frame(nbt@reductions$umap@cell.embeddings)
nbt.out$SAMPLE=rownames(nbt.out)
nbt.out$Stage=sapply( strsplit(as.character(rownames(nbt.out)), "_"), "[[", 1 )
nbt.out$Embryo=sapply( strsplit(as.character(rownames(nbt.out)), "_"), "[[", 3 )

nbt.out$Location=sapply( strsplit(as.character(rownames(nbt.out)), "_"), "[[", 4 )

nbt.out$cell_type <-NA
table(nbt.out$Embryo)
#mycolors=c("blue","brown","purple","darkgreen","RoyalBlue1","grey","DarkGoldenrod3","springgreen","deeppink1","firebrick1","cyan","deepskyblue4","pink","lightgoldenrod3","hotpink","steelblue1","PaleGreen","indianred3","lightblue","mediumpurple1")
mycolors <- c("pink","gray","aquamarine2","green","red","navy","purple","cyan","darkorange2","cornflowerblue","deeppink","orchid1","darkturquoise","deepskyblue4","dimgray","palegreen","khaki1",
              "darkolivegreen1","slateblue","yellow","darkred","darksalmon","blue","gold4","forestgreen")
nbt.out$Embryo <- factor(nbt.out$Embryo,levels = c("E1","E2","E3","E4","E5","E7","E8","E10","E11","E12","E13","E24"))
library(ggplot2)
p<-ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= Embryo,shape=Location))
p<-p+geom_point(size=1.5)+theme_bw()+
  theme(panel.grid=element_blank())+scale_color_manual(values = mycolors,limits=c("E1","E2","E3","E4","E5","E6","E7","E8","E9","E10","E11","E12","E13","E14","E15","E16","E17","E18","E19","E20","E21","E22","E23","E24","E25"))+scale_shape_manual(values=c(17,15),limits=c("M","P"))#+scale_x_continuous(limits = c(-40,40))+scale_y_continuous(limits = c(-40,40))
p
# S2.EPI <- read.table("../ALL_data/S2.EPI_name.txt",header = T)
# S2.PE <- read.table("../ALL_data/S2.PE_name.txt",header = T)
# S2.MTE1 <- read.table("../ALL_data/S2.MTE1_name.txt",header = T)
# S2.MTE2 <- read.table("../ALL_data/S2.MTE2_name.txt",header = T)
# S2.MTE3 <- read.table("../ALL_data/S2.MTE3_name.txt",header = T)
# S2.PTE <- read.table("../ALL_data/S2.PTE_name.txt",header = T)
# nbt.out[which(rownames(nbt.out) %in% S2.EPI$name), ]$cell_type=c("S2.EPI")
# nbt.out[which(rownames(nbt.out) %in% S2.PE$name), ]$cell_type=c("S2.PE")
# nbt.out[which(rownames(nbt.out) %in% S2.MTE1$name), ]$cell_type=c("S2.MTE1")
# nbt.out[which(rownames(nbt.out) %in% S2.MTE2$name), ]$cell_type=c("S2.MTE2")
# nbt.out[which(rownames(nbt.out) %in% S2.MTE3$name), ]$cell_type=c("S2.MTE3")
# nbt.out[which(rownames(nbt.out) %in% S2.PTE$name), ]$cell_type=c("S2.PTE")
#############S2.TE1 
C4_name <- data.frame(name=CellsByIdentities(nbt,4))
C4_tpm <- S2.tpm[ ,which(colnames(S2.tpm) %in% C4_name$X4)]
colnames(C4_tpm) <- paste("S2.PE",colnames(C4_tpm),sep = "_")
write.table(C4_tpm,"S2.PE_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C4_tpm.log <- log2(C4_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C4_name$X4), ]$cell_type <- c("S2.PE")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank())
################S2.TE2
C3_name <- data.frame(name=CellsByIdentities(nbt,2))
C3_tpm <- S2.tpm[ ,which(colnames(S2.tpm) %in% C3_name$X2)]
colnames(C3_tpm) <- paste("S2.MTE3",colnames(C3_tpm),sep = "_")
write.table(C3_tpm,"S2.MTE3_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C3_tpm.log <- log2(C3_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C3_name$X2), ]$cell_type <- c("S2.MTE3")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank())
#################S2.S3.tTE
C4_name <- data.frame(name=CellsByIdentities(nbt,3))
C4_tpm <- S2.tpm[ ,which(colnames(S2.tpm) %in% C4_name$X3)]
colnames(C4_tpm) <- paste("S2.MTE2",colnames(C4_tpm),sep = "_")
write.table(C4_tpm,"S2.MTE2_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C4_tpm.log <- log2(C4_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C4_name$X3), ]$cell_type <- c("S2.MTE2")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=1.5)+theme_bw()+theme(panel.grid=element_blank())+scale_shape_manual(values=c(17,15),limits=c("M","P"))
#################S2.TE3
C1_name <- data.frame(name=CellsByIdentities(nbt,1))
C1_tpm <- S2.tpm[ ,which(colnames(S2.tpm) %in% C1_name$X1)]
#colnames(C1_tpm) <- paste("S2.TE4",colnames(C1_tpm),sep = "_")
colnames(C1_tpm) <- paste("S2.PTE",colnames(C1_tpm),sep = "_")
write.table(C1_tpm,"S2.PTE_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C1_tpm.log <- log2(C1_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C1_name$X1), ]$cell_type <- c("S2.PTE")
#################S2.S3.tTE
C4_name <- data.frame(name=CellsByIdentities(nbt,0))
C4_tpm <- S2.tpm[ ,which(colnames(S2.tpm) %in% C4_name$X0)]
colnames(C4_tpm) <- paste("S2.MTE1",colnames(C4_tpm),sep = "_")
write.table(C4_tpm,"S2.MTE1_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C4_tpm.log <- log2(C4_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C4_name$X0), ]$cell_type <- c("S2.MTE1")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=1.5)+theme_bw()+theme(panel.grid=element_blank())+scale_shape_manual(values=c(17,15),limits=c("M","P"))
#################S2.TE3
C1_name <- data.frame(name=CellsByIdentities(nbt,5))
C1_tpm <- S2.tpm[ ,which(colnames(S2.tpm) %in% C1_name$X5)]
#colnames(C1_tpm) <- paste("S2.TE4",colnames(C1_tpm),sep = "_")
colnames(C1_tpm) <- paste("S2.EPI",colnames(C1_tpm),sep = "_")
write.table(C1_tpm,"S2.EPI_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C1_tpm.log <- log2(C1_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C1_name$X5), ]$cell_type <- c("S2.EPI")






pdf("../Figure/Figure4b.5.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3.5)+theme_bw()+theme(panel.grid=element_blank())+scale_shape_manual(values=c(17,15),limits=c("M","P"))+scale_color_manual(values = c("S2.MTE1"="#25BFF7","S2.MTE2"="orchid1","S2.MTE3"="#FFF64A","S2.PTE"="limegreen","S2.EPI"="#f4a790","S2.PE"="burlywood4"),limits=c("S2.MTE1","S2.MTE2","S2.MTE3","S2.PTE","S2.EPI","S2.PE"))
dev.off()
write.table(nbt.out,"../Table/Figure2a_S2_coordinate.xls",quote = F,col.names = T,row.names = F,sep = "\t")

pdf("../Figure/revise_Human_S2_umap.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3.5)+theme_bw()+theme(panel.grid=element_blank())+scale_shape_manual(values=c(17,15),limits=c("M","P"))+scale_color_manual(values = c("S2.MTE1"="#25BFF7","S2.MTE2"="orchid1","S2.MTE3"="#FFF64A","S2.PTE"="limegreen","S2.EPI"="#f4a790","S2.PE"="burlywood4"),limits=c("S2.MTE1","S2.MTE2","S2.MTE3","S2.PTE","S2.EPI","S2.PE"))
dev.off()

marker.plot <- function(x){
  nbt.out[ ,x]=t(S2.tpm.log[x,])
  gene=nbt.out[ ,x]
  ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=0.5)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                               axis.ticks.length = unit(0.15, "cm"),
                               panel.background = element_rect(fill = "white", colour = "white"),
                               axis.title.y= element_text(size=9),
                               title= element_text(size=9),
                               legend.text= element_text(size=9),
                               axis.text.x = element_text(size=9))+
    scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
    labs(title=x)+theme(panel.grid=element_blank())+theme_bw()+
    theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
    guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
}

CDX2 <- marker.plot("CDX2")
KRT18 <- marker.plot("KRT18")
GATA2 <- marker.plot("GATA2")
GATA3 <- marker.plot("GATA3")
CLDN10 <-marker.plot("CLDN10")
CLDN4 <-marker.plot("CLDN4")
FUT4 <- marker.plot("FUT4")
library(easyGgplot2)
pdf("../Figure/Figure4b.1.pdf",width = 13,height = 6)
ggplot2.multiplot(CDX2,KRT18,GATA2,GATA3,CLDN10,FUT4,cols=3)
dev.off()

CLDN3 <- marker.plot("CLDN3")
FGFR4 <- marker.plot("FGFR4")
marker.plot("MUC16")
marker.plot("CCNE1")
marker.plot("XIST")
marker.plot("TSIX")
marker.plot("IDH1")
marker.plot("NDUFAB1")
SOX2 <- marker.plot("SOX2")
POU5F1 <- marker.plot("POU5F1")
NANOG <- marker.plot("NANOG")


pdf("../Figure/Figure4b.2.pdf",width = 13,height = 3)
ggplot2.multiplot(SOX2,NANOG,POU5F1,cols=3)
dev.off()

GATA4 <- marker.plot("GATA4")
PDGFRA <- marker.plot("PDGFRA")
FOXA2 <- marker.plot("FOXA2")
pdf("../Figure/Figure4b.3.pdf",width = 13,height = 3)
ggplot2.multiplot(GATA4,PDGFRA,FOXA2,cols=3)
dev.off()


pdf("../Figure/Figure4b.4.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= Embryo,shape=Location))+geom_point(size=1)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("pink","gray","aquamarine2","green","red","purple","cyan","cornflowerblue","deeppink","orchid1","darkturquoise","gold4"),limits=c("E1","E2","E3","E4","E5","E7","E8","E10","E11","E12","E13","E24"))+scale_shape_manual(values=c(17,15),limits=c("M","P"))
dev.off()


#############FIG S4a
TE_chart <- data.frame(CDX2=t(S2.tpm.log["CDX2",]))
TE_chart[ ,"GATA2"]=t(S2.tpm.log["GATA2",])
TE_chart[ ,"GATA3"]=t(S2.tpm.log["GATA3",])
nbt.out[ ,"TE"]=rowMeans(TE_chart)
TE <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=TE))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.7,name= "log2(TPM/10+1)")+
  labs(title="TE")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

PE_chart <- data.frame(GATA4=t(S2.tpm.log["GATA4",]))
PE_chart[ ,"PDGFRA"]=t(S2.tpm.log["PDGFRA",])
PE_chart[ ,"FOXA2"]=t(S2.tpm.log["FOXA2",])
nbt.out[ ,"PE"]=rowMeans(PE_chart)
PE <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=PE))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2,name= "log2(TPM/10+1)")+
  labs(title="PE")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

EPI_chart <- data.frame(SOX2=t(S2.tpm.log["SOX2",]))
EPI_chart[ ,"NANOG"]=t(S2.tpm.log["NANOG",])
EPI_chart[ ,"POU5F1"]=t(S2.tpm.log["POU5F1",])
nbt.out[ ,"EPI"]=rowMeans(EPI_chart)
EPI <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=EPI))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="EPI")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
library(easyGgplot2)
pdf("../Figure/Human_S2_markers_exp.pdf",width = 6,height = 2)
ggplot2.multiplot(TE,EPI,PE,cols=3)
dev.off()



nbt.out[ ,"CCR7"]=t(S2.tpm.log["CCR7",])
CCR7 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=CCR7))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="CCR7")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"CYP19A1"]=t(S2.tpm.log["CYP19A1",])
CYP19A1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=CYP19A1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="CYP19A1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"DLX5"]=t(S2.tpm.log["DLX5",])
DLX5 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=DLX5))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="DLX5")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"OVOL1"]=t(S2.tpm.log["OVOL1",])
OVOL1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=OVOL1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="OVOL1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"GCM1"]=t(S2.tpm.log["GCM1",])
GCM1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=GCM1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2,name= "log2(TPM/10+1)")+
  labs(title="GCM1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"MUC15"]=t(S2.tpm.log["MUC15",])
MUC15 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=MUC15))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="MUC15")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
library(easyGgplot2)
pdf("../Figure/Revise_Human_FigS4_PTEmarkers_exp.pdf",width = 6,height = 4)
ggplot2.multiplot(CCR7,CYP19A1,DLX5,OVOL1,GCM1,MUC15,cols=3)
dev.off()

