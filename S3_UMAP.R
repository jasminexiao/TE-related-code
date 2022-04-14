setwd("E:/Qiaolab/liudandan/human/ALL_data")
TPM <- read.table("qc_tpm(1127).txt",as.is = TRUE)
#TPM.log <- log2(TPM/10+1)
Stage.class <- sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 1 )
S3.tpm <- TPM[ ,which(Stage.class %in% c("S3"))]
S3.tpm.log <- log2(S3.tpm/10+1)
#Embryo.class <- sapply( strsplit(as.character(colnames(S3.tpm)), "_"), "[[", 3 )
#S3.tpm.no.E25 <- S3.tpm[ ,-which(Embryo.class %in% c("E25"))]
#S3.tpm.no.E25.log <- log2(S3.tpm.no.E25/10+1)
library(Seurat)
library(dplyr)
library(Seurat)
library(dplyr)
nbt <- CreateSeuratObject(S3.tpm, project = "S3", min.cells = 3, min.features = 200)
nbt
nbt[["percent.mt"]] <- PercentageFeatureSet(nbt, pattern = "^MT-")
VlnPlot(nbt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
nbt <- NormalizeData(nbt, normalization.method = "LogNormalize", scale.factor = 10000)
nbt <- NormalizeData(nbt)
nbt <- FindVariableFeatures(nbt, selection.method = "vst", nfeatures = 1800)

# plot variable features with and without labels
VariableFeaturePlot(nbt)
nbt <- ScaleData(nbt)
nbt <- RunPCA(nbt, features = VariableFeatures(object = nbt))
DimPlot(nbt, reduction = "pca")
#DimHeatmap(nbt, dims = 1:15, cells = 500, balanced = TRUE)
nbt <- JackStraw(nbt, num.replicate = 100)
nbt <- ScoreJackStraw(nbt, dims = 1:20)
JackStrawPlot(nbt, dims = 1:20)
ElbowPlot(nbt)
nbt <- FindNeighbors(nbt, dims = 1:4)
nbt <- FindClusters(nbt, resolution = 0.3)
nbt <- RunUMAP(nbt, dims = 1:4)
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
nbt.out$Embryo <- factor(nbt.out$Embryo,levels = c("E15","E16","E17","E18","E23"))
library(ggplot2)
p<-ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= Embryo,shape=Location))
p<-p+geom_point(size=1)+theme_bw()+
  theme(panel.grid=element_blank())+scale_color_manual(values = c("dimgray","palegreen","khaki1","darkolivegreen1","blue"),limits=c("E15","E16","E17","E18","E23"))#+scale_x_continuous(limits = c(-40,40))+scale_y_continuous(limits = c(-40,40))
p

#############S3.TE1 
C0_name <- data.frame(name=CellsByIdentities(nbt,1))
C0_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C0_name$X1)]
colnames(C0_tpm) <- paste("S3.PTE",colnames(C0_tpm),sep = "_")
write.table(C0_tpm,"S3.PTE_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C0_tpm.log <- log2(C0_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C0_name$X1), ]$cell_type <- c("S3.PTE")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank())
################S3.TE2
C1_name <- data.frame(name=CellsByIdentities(nbt,2))
C1_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C1_name$X2)]
colnames(C1_tpm) <- paste("S3.coISK",colnames(C1_tpm),sep = "_")
write.table(C1_tpm,"S3.coISK_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C1_tpm.log <- log2(C1_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C1_name$X2), ]$cell_type <- c("S3.coISK")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3)+theme_bw()+theme(panel.grid=element_blank())
#################S3.S3.tTE
C2_name <- data.frame(name=CellsByIdentities(nbt,0))
C2_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C2_name$X0)]
colnames(C2_tpm) <- paste("S3.MTE",colnames(C2_tpm),sep = "_")
write.table(C2_tpm,"S3.MTE_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C2_tpm.log <- log2(C2_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C2_name$X0), ]$cell_type <- c("S3.MTE")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=1.5)+theme_bw()+theme(panel.grid=element_blank())#+scale_shape_manual(values=c(17,15),limits=c("M","P"))
#################S3.TE3
C3 <- data.frame(CellsByIdentities(nbt,3))
C3_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C3$X3)]
nbt.C3 <- CreateSeuratObject(C3_tpm, project = "S3", min.cells = 3, min.features = 200)
nbt.C3
nbt.C3[["percent.mt"]] <- PercentageFeatureSet(nbt.C3, pattern = "^MT-")
VlnPlot(nbt.C3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
nbt.C3 <- NormalizeData(nbt.C3, normalization.method = "LogNormalize", scale.factor = 10000)
nbt.C3 <- FindVariableFeatures(nbt.C3, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
VariableFeaturePlot(nbt.C3)
nbt.C3 <- ScaleData(nbt.C3)
nbt.C3 <- RunPCA(nbt.C3, npcs=49,features = VariableFeatures(object = nbt.C3))
DimPlot(nbt.C3, reduction = "pca")
#DimHeatmap(nbt.C3, dims = 1:15, cells = 500, balanced = TRUE)
nbt.C3 <- JackStraw(nbt.C3, num.replicate = 100)
nbt.C3 <- ScoreJackStraw(nbt.C3, dims = 1:20)
JackStrawPlot(nbt.C3, dims = 1:20)
ElbowPlot(nbt.C3)
nbt.C3 <- FindNeighbors(nbt.C3, dims = 1:4)
nbt.C3 <- FindClusters(nbt.C3, resolution = 1)
nbt.C3 <- RunUMAP(nbt.C3, dims = 1:4)
DimPlot(nbt.C3, reduction = "umap",label = T)

C2.0_name <- data.frame(name=CellsByIdentities(nbt.C3,0))
nbt.out[which(nbt.out$SAMPLE %in% C2.0_name$X0), ]$cell_type <- c("S3.EPI")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=1.5)+theme_bw()+theme(panel.grid=element_blank())

C2.1_name <- data.frame(name=CellsByIdentities(nbt.C3,1))
nbt.out[which(nbt.out$SAMPLE %in% C2.1_name$X1), ]$cell_type <- c("S3.PE")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=1.5)+theme_bw()+theme(panel.grid=element_blank())

# pdf("../Figure/Figure4b.5.pdf",width = 6.5,height = 5)
# ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3.5)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("S3.MTE"="#BB47F5","S3.PTE"="#ADE7EC","S3.EPI"="#bf9a8e","S3.PE"="darkgoldenrod1","S3.coISK"="darkgreen"),limits=c("S3.MTE","S3.PTE","S3.EPI","S3.PE","S3.coISK"))
# dev.off()
write.table(nbt.out,"../Table/Figure2a_S3_coordinate.xls",quote = F,col.names = T,row.names = F,sep = "\t")

pdf("../Figure/revise_Human_S3_umap.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3.2)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("S3.MTE"="#BB47F5","S3.PTE"="#ADE7EC","S3.EPI"="#bf9a8e","S3.PE"="darkgoldenrod1","S3.coISK"="darkgreen"),limits=c("S3.MTE","S3.PTE","S3.EPI","S3.PE","S3.coISK"))
dev.off()


pdf("../Figure/H33_UMAP.pdf",width = 6.5,height = 5)
ggplot(nbt.H33.out,aes(x=UMAP_1,y=UMAP_2,color= Stage,shape=Location))+geom_point(size=1.5)+theme_bw()
dev.off()
marker.plot <- function(x){
  nbt.out[ ,x]=t(S3.tpm.log[x,])
  gene=nbt.out[ ,x]
  ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=1)+theme(axis.line = element_line(size = 0.5, colour = "black"),
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
pdf("../Figure/PSG5.pdf",width = 4,height = 2.5)
marker.plot("PSG5")
dev.off()
pdf("../Figure/PSG9.pdf",width = 4,height = 2.5)
marker.plot("PSG9")
dev.off()
library(easyGgplot2)
pdf("../Figure/Figure4c.1.pdf",width = 13,height = 6)
ggplot2.multiplot(CDX2,KRT18,GATA2,GATA3,CLDN10,FUT4,cols=3)
dev.off()

CLDN3 <- marker.plot("CLDN3")
FGFR4 <- marker.plot("FGFR4")
marker.plot("CD53")
marker.plot("CCNE1")
marker.plot("XIST")
marker.plot("TSIX")
marker.plot("IDH1")
marker.plot("NDUFAB1")
SOX2 <- marker.plot("SOX2")
POU5F1 <- marker.plot("POU5F1")
NANOG <- marker.plot("NANOG")


pdf("../Figure/Figure4c.2.pdf",width = 13,height = 3)
ggplot2.multiplot(SOX2,NANOG,POU5F1,cols=3)
dev.off()

GATA4 <- marker.plot("GATA4")
PDGFRA <- marker.plot("PDGFRA")
FOXA2 <- marker.plot("FOXA2")
pdf("../Figure/Figure4c.3.pdf",width = 13,height = 3)
ggplot2.multiplot(GATA4,PDGFRA,FOXA2,cols=3)
dev.off()


pdf("../Figure/Figure4c.4.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= Embryo,shape=Location))+geom_point(size=1)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("dimgray","palegreen","khaki1","darkolivegreen1","blue"),limits=c("E15","E16","E17","E18","E23"))#+scale_shape_manual(values=c(17,15),limits=c("M","P"))
dev.off()

#############FIG S4a
S3.tpm.log1 <- S3.tpm.log[ ,as.vector(rownames(nbt.out))]
TE_chart <- data.frame(CDX2=t(S3.tpm.log1["CDX2",]))
TE_chart[ ,"GATA2"]=t(S3.tpm.log1["GATA2",])
TE_chart[ ,"GATA3"]=t(S3.tpm.log1["GATA3",])
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

PE_chart <- data.frame(GATA4=t(S3.tpm.log1["GATA4",]))
PE_chart[ ,"PDGFRA"]=t(S3.tpm.log1["PDGFRA",])
PE_chart[ ,"FOXA2"]=t(S3.tpm.log1["FOXA2",])
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
  scale_colour_gradient2(low="blue",high="red",midpoint = 1.7,name= "log2(TPM/10+1)")+
  labs(title="PE")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

EPI_chart <- data.frame(SOX2=t(S3.tpm.log1["SOX2",]))
EPI_chart[ ,"NANOG"]=t(S3.tpm.log1["NANOG",])
EPI_chart[ ,"POU5F1"]=t(S3.tpm.log1["POU5F1",])
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
  scale_colour_gradient2(low="blue",high="red",midpoint = 2,name= "log2(TPM/10+1)")+
  labs(title="EPI")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
library(easyGgplot2)
pdf("../Figure/Human_S3_markers_exp.pdf",width = 6,height = 2)
ggplot2.multiplot(TE,EPI,PE,cols=3)
dev.off()
pdf("../Figure/revise_Human_S3_markers_exp.pdf",width = 6,height = 2)
ggplot2.multiplot(TE,EPI,PE,cols=3)
dev.off()

S3.tpm.log1 <- S3.tpm.log[ ,as.vector(rownames(nbt.out))]
nbt.out[ ,"HAND1"]=t(S3.tpm.log1["HAND1",])
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=HAND1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 1.7,name= "log2(TPM/10+1)")+
  labs(title="HAND1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"GCM1"]=t(S3.tpm.log1["GCM1",])
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=GCM1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.6,name= "log2(TPM/10+1)")+
  labs(title="GCM1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"PGF"]=t(S3.tpm.log1["PGF",])
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=PGF))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 4,name= "log2(TPM/10+1)")+
  labs(title="PGF")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
####################autocluster
nbt.cluster <- FindClusters(nbt, reduction.type = "pca", dims.use = 1:7, resolution = 0.44, print.output = F)
UMAPPlot(nbt.cluster)

##########ICM subtype
C3_name <- data.frame(name=WhichCells(nbt.cluster,3))
C3_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C3_name$name)]
colnames(C3_tpm) <- paste("S3.EPI",colnames(C3_tpm),sep = "_")
write.table(C3_tpm,"S3.EPI_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C3_tpm.log <- log2(C3_tpm/10+1)
C3_UMAP <- nbt.out[which(nbt.out$SAMPLE %in% C3_name$name), ]

nbt.out[which(nbt.out$SAMPLE %in% C3_name$name), ]$cell_type <- c("S3.EPI")

ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=0.5)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("#37ACCB","#FB9100","#B0C915"))
#############S3.TE1 
C0_name <- data.frame(name=WhichCells(nbt.cluster,0))
C0_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C0_name$name)]
C0_tpm.log <- log2(C0_tpm/10+1)
C0_UMAP <- nbt.out[which(nbt.out$SAMPLE %in% C0_name$name), ]

nbt.out[which(nbt.out$SAMPLE %in% (C0_UMAP[which(C0_UMAP$UMAP_1 >-10), ]$SAMPLE)), ]$cell_type <- c("S3.MTE")
C0_tpm.S3.te1 <- C0_tpm[ ,which(colnames(C0_tpm) %in% nbt.out[(nbt.out$cell_type %in% c("S3.MTE")), ]$SAMPLE)]
colnames(C0_tpm.S3.te1) <- paste("S3.MTE",colnames(C0_tpm.S3.te1),sep = "_")
write.table(C0_tpm.S3.te1,"S3.TE1_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
nbt.out[which(nbt.out$SAMPLE %in% (C0_UMAP[which(C0_UMAP$UMAP_1 < -10), ]$SAMPLE)), ]$cell_type <- c("S3.PE")
C0_tpm.PE <- C0_tpm[ ,which(colnames(C0_tpm) %in% nbt.out[(nbt.out$cell_type %in% c("S3.PE")), ]$SAMPLE)]
colnames(C0_tpm.PE) <- paste("S3.PE",colnames(C0_tpm.PE),sep = "_")
#write.table(C0_tpm.PE,"S3.PE_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=0.5)+theme_bw()+theme(panel.grid=element_blank())
################ISH
C2_name <- data.frame(name=WhichCells(nbt.cluster,2))
C2_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C2_name$name)]
colnames(C2_tpm) <- paste("S3.ISH",colnames(C2_tpm),sep = "_")
write.table(C2_tpm,"S3.ISH_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C2_tpm.log <- log2(C2_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C2_name$name), ]$cell_type <- c("S3.ISH")
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=0.5)+theme_bw()+theme(panel.grid=element_blank())

#################S3.TE2
C1_name <- data.frame(name=WhichCells(nbt.cluster,1))
C1_tpm <- S3.tpm[ ,which(colnames(S3.tpm) %in% C1_name$name)]
#colnames(C1_tpm) <- paste("S3.TE2",colnames(C1_tpm),sep = "_")
colnames(C1_tpm) <- paste("S3.PTE",colnames(C1_tpm),sep = "_")
write.table(C1_tpm,"S3.PTE_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
C1_tpm.log <- log2(C1_tpm/10+1)
nbt.out[which(nbt.out$SAMPLE %in% C1_name$name), ]$cell_type <- c("S3.PTE")
pdf("../Figure/Figure4c.5.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3.5)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("S3.MTE"="#BB47F5","S3.PTE"="#ADE7EC","S3.EPI"="#bf9a8e","S3.PE"="darkgoldenrod1","S3.ISH"="darkgreen"),limits=c("S3.MTE","S3.PTE","S3.EPI","S3.PE","S3.ISH"))
dev.off()
write.table(nbt.out,"../Table/Figure2a_S3_coordinate.xls",quote = F,col.names = T,row.names = F,sep = "\t")
S3_tpm.new <- cbind(C0_tpm.S3.te1,C1_tpm,C3_tpm,C0_tpm.PE,C2_tpm)
colnames(S3_tpm.new)
S3_MTE_PTE_ISH_tpm <- cbind(C0_tpm.S3.te1,C1_tpm,C2_tpm)
colnames(S3_MTE_PTE_ISH_tpm)
##################
#########findallmarkers
TPM.seq.class <- sapply( strsplit(as.character(colnames(S3_tpm.new)), "_"), "[[", 1 )
TPM.seq.class <- factor(TPM.seq.class,levels = c("S3.EPI","S3.PE","S3.MTE","S3.PTE","S3.ISH"))
TPM.seq <- S3_tpm.new[ ,order(TPM.seq.class)]
library(Seurat)

nbt.stage2 <- CreateSeuratObject(raw.data =TPM.seq, min.cells = 3, min.genes = 2000)
nbt.stage2 <- NormalizeData(object = nbt.stage2, normalization.method = "LogNormalize",scale.factor = 10000)
nbt.stage2 <- FindVariableGenes(object = nbt.stage2, mean.function = ExpMean, dispersion.function = LogVMR,
                                x.low.cutoff = 0.4, y.cutoff = 0.5)
length(nbt.stage2@var.genes)
nbt.stage2.all.markers=FindAllMarkers(nbt.stage2,test.use = "roc")
nbt.stage2.markers.use=subset(nbt.stage2.all.markers,avg_diff>0&power>0.4)
nbt.stage2.markers.use$cluster <- factor(nbt.stage2.markers.use$cluster,levels = c("S3.EPI","S3.PE","S3.MTE","S3.PTE","S3.ISH"))
nbt.stage2.markers.use <- nbt.stage2.markers.use[order(nbt.stage2.markers.use$cluster), ]
write.table(nbt.stage2.markers.use,"../Table/Stage3_all_type_markers.xls",quote = F,col.names = T,row.names = F,sep = "\t")

nbt.stage.TE.markers.use.gene <- as.vector(nbt.stage2.markers.use$gene)

nbt.stage.TE.markers.use.tpm <- TPM.seq[nbt.stage.TE.markers.use.gene,]
mean <- rowMeans(nbt.stage.TE.markers.use.tpm)
sd_value <- apply(nbt.stage.TE.markers.use.tpm,1,sd)

z_score <- (nbt.stage.TE.markers.use.tpm-mean)/sd_value
z_score[z_score>=1.5] <- 1.5
z_score[z_score<= -1.5] <- -1.5

ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )
aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
#colnames(aaka3)<-c("cluster")
colnames(aaka22)<-c("Class")

aka4=list(Class=c("S3.TE1"="#BB47F5","S3.TE2"="#ADE7EC","S3.EPI"="#bf9a8e","S3.PE"="darkgoldenrod1","ISH"="darkgreen"))
library('RColorBrewer')
library("pheatmap")
pheatmap( z_score,annotation_col =aaka22,annotation_colors = aka4,
          color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(50), 
          #breaks = c(seq(-2,-1.1,length=15),seq(-1,1,length=70),seq(1.1,2,length=15)),
          cluster_row = F,cluster_col =F,show_rownames = F, 
          show_colnames =F, scale = "none",legend =T,border_color =T,
          cellwidth = 0.5,cellheight = 0.1,fontsize_row=8,
          fontsize=8,fontsize_col = 5,file="../Figure/Figure4c.6.pdf"
)

##############TOP10 new markers
nbt.stage2.markers.use.select <- nbt.stage2.markers.use[which(nbt.stage2.markers.use$pct.1>0.6 & nbt.stage2.markers.use$pct.2<0.4), ]
table(nbt.stage2.markers.use.select$cluster)
MTE.markers <- nbt.stage2.markers.use.select[which(nbt.stage2.markers.use.select$cluster %in% c("S3.MTE")), ]
PTE.markers <- nbt.stage2.markers.use.select[which(nbt.stage2.markers.use.select$cluster %in% c("S3.PTE")), ]
ISH.markers <- nbt.stage2.markers.use.select[which(nbt.stage2.markers.use.select$cluster %in% c("S3.ISH")), ]
ISH.markers.SEL <- ISH.markers[1:10, ]
TOP10_markers <- rbind(MTE.markers,PTE.markers,ISH.markers.SEL)
TOP10_markers.gene <- as.vector(TOP10_markers$gene)

TOP10_markers.tpm <- S3_MTE_PTE_ISH_tpm[TOP10_markers.gene,]
mean <- rowMeans(TOP10_markers.tpm)
sd_value <- apply(TOP10_markers.tpm,1,sd)

z_score <- (TOP10_markers.tpm-mean)/sd_value
z_score[z_score>=0.5] <- 0.5
z_score[z_score<= -0.5] <- -0.5

ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )
aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
#colnames(aaka3)<-c("cluster")
colnames(aaka22)<-c("Class")

aka4=list(Class=c("S3.MTE"="#BB47F5","S3.PTE"="#ADE7EC","S3.ISH"="darkgreen"))
library('RColorBrewer')
library("pheatmap")
pheatmap( z_score,annotation_col =aaka22,annotation_colors = aka4,
          #color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(50), 
          #breaks = c(seq(-2,-1.1,length=15),seq(-1,1,length=70),seq(1.1,2,length=15)),
          color = colorRampPalette(c("blue", "black", "yellow"))(100),
          cluster_row = F,cluster_col =F,show_rownames = T, 
          show_colnames =F, scale = "none",legend =T,border_color =T,
          cellwidth = 0.6,cellheight = 9,fontsize_row=8,
          fontsize=8,fontsize_col = 5,file="../Figure/TOP10_MTE_PTE_ISH_marker_heatmap.pdf"
)
###############NAME
S3.MTE_name=data.frame(name=nbt.out[which(nbt.out$cell_type %in% c("S3.MTE")), ]$SAMPLE)
S3.EPI_name=data.frame(name=nbt.out[which(nbt.out$cell_type %in% c("S3.EPI")), ]$SAMPLE)
S3.PE_name=data.frame(name=nbt.out[which(nbt.out$cell_type %in% c("S3.PE")), ]$SAMPLE)
S3.ISH_name=data.frame(name=nbt.out[which(nbt.out$cell_type %in% c("S3.ISH")), ]$SAMPLE)
S3.PTE_name=data.frame(name=nbt.out[which(nbt.out$cell_type %in% c("S3.PTE")), ]$SAMPLE)


write.table(S3.MTE_name,"S3.MTE_name.txt",quote = F,col.names = T,row.names = F)
write.table(S3.EPI_name,"S3.EPI_name.txt",quote = F,col.names = T,row.names = F)
write.table(S3.PE_name,"S3.PE_name.txt",quote = F,col.names = T,row.names = F)
write.table(S3.ISH_name,"S3.ISH_name.txt",quote = F,col.names = T,row.names = F)
write.table(S3.PTE_name,"S3.PTE_name.txt",quote = F,col.names = T,row.names = F)

#################CELLPHONE
ISK_tpm1 <- read.table("ISK_1.UMI_TPM_no_ERCC.xls",header = T,sep = "\t")
rownames(ISK_tpm1) <- ISK_tpm1$Gene
ISK_tpm1$Gene <- NULL
ISK_tpm2 <- read.table("ISK_2.UMI_TPM_no_ERCC.xls",header = T,sep = "\t")
rownames(ISK_tpm2) <- ISK_tpm2$Gene
ISK_tpm2$Gene <- NULL
ISK_TPM <- cbind(ISK_tpm1,ISK_tpm2)
colnames(ISK_TPM)
S3_TPM <- cbind(C3_tpm,C0_tpm.PE,C0_tpm.S3.te1,C1_tpm,C2_tpm,ISK_TPM)
colnames(S3_TPM)
ensambl_id <- read.table("G:/tang lab/huyuqiong/01.RetinaProject/cyd_analysis/ensembl_symbol_geneid.txt",header = T,sep = "\t")
S3_TPM.sel <- S3_TPM[which(rownames(S3_TPM) %in% ensambl_id$gene_name), ]
seq <- as.vector(ensambl_id$gene_name)
S3_TPM.sel <- S3_TPM.sel[seq, ]
rownames(S3_TPM.sel) <- ensambl_id$ensembl
S3_TPM.sel.log <- log2(S3_TPM.sel/10+1)
S3_TPM.sel.log[1:10,1:3]
write.table(S3_TPM.sel.log,"S3_TPM.sel.Ensembl.txt",quote = F,col.names = T,row.names = T,sep="\t")

S3_chart <- data.frame(Cell=colnames(S3_TPM),cell_type=sapply( strsplit(as.character(colnames(S3_TPM)), "_"), "[[", 1))
write.table(S3_chart,"S3_chart.txt",quote = F,col.names = T,row.names = F,sep="\t")

#######################ADD H33 SAMPLE#####################
H33_tpm <- read.table("H33_tpm.select.txt",header = T)
colnames(H33_tpm) <- paste("H33",colnames(H33_tpm),sep = "_")
colnames(TPM.seq)
S3_tpm.new_H33 <- cbind(S3_tpm.new,H33_tpm)
write.table(S3_tpm.new_H33,"human_S3_H33_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
S3_tpm.new_H33.log <- log2(S3_tpm.new_H33/10+1)
library(Seurat)
library(dplyr)
nbt.H33 <- CreateSeuratObject(raw.data =S3_tpm.new_H33, min.cells = 3, min.genes = 2000)
nbt.H33 <- NormalizeData(object = nbt.H33, normalization.method = "LogNormalize",scale.factor = 10000)
nbt.H33 <- FindVariableGenes(object = nbt.H33, mean.function = ExpMean, dispersion.function = LogVMR,
                             x.low.cutoff = 0.4, y.cutoff = 0.4)####0.4 & 0.4
length(nbt.H33@var.genes)
nbt.H33@var.genes <- nbt@var.genes
nbt.H33 <- ScaleData(object = nbt.H33)##, vars.to.regress = c("nUMI")
nbt.H33 <- RunPCA(object = nbt.H33, pc.genes = nbt.H33@var.genes, do.print = TRUE, pcs.print = 1:5, 
                  genes.print = 5)
#pdf("../new_plot/Stage_PCA.pdf",width = 5,height = 4)
#DimPlot(nbt.H33,1,3, reduction.use = "pca",cols.use=c("#37ACCB","#FB9100","#B0C915"))
#dev.off()
#pdf("../new_plot/Stage_PCA_heatmap.pdf",width = 12,height = 14)
nbt.H33 <- ProjectPCA(object = nbt.H33, do.print = FALSE)
#PCHeatmap(object = nbt.H33, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
#          label.columns = FALSE, use.full = FALSE)
#dev.off()
#VizPCA(nbt.H33,1:3)
#PCA.gene <- nbt.H33@dr$pca
#write.table(PCA.gene,"PCA_gene.xls",quote = F,sep = "\t",row.names = T,col.names = T)
#DimPlot(nbt.H33,4,9, reduction.use = "pca",cols.use=c("#37ACCB","#FB9100","#B0C915"))
nbt.H33=JackStraw(nbt.H33,num.replicate = 100,display.progress = TRUE)
JackStrawPlot(nbt.H33,PCs = 1:20)
PCElbowPlot(object = nbt.H33)
nbt.H33 <- FindClusters(object = nbt.H33, dims.use = 1:7, reduction.type = "pca",resolution = 0.5,#force.recalc=TRUE,dim.use=10
                        print.output = 0, save.SNN = TRUE)#  
#nbt.H33 <- RunUMAP(nbt.H33, dims.use = 1:12,max_iter=5000)
PrintFindClustersParams(object = nbt.H33)
nbt.H33 <- RunUMAP(object = nbt.H33, dims.use = 1:7, do.fast = TRUE)
pdf("../Figure/S3_H33_auto_UMAP.pdf",width = 3.7,height = 3)
UMAPPlot(nbt.H33)
dev.off()

nbt.H33.out=data.frame(nbt.H33@dr$UMAP@cell.embeddings)
nbt.H33.out$SAMPLE=rownames(nbt.H33.out)
nbt.H33.out$Stage=sapply( strsplit(as.character(rownames(nbt.H33.out)), "_"), "[[", 1 )
nbt.H33.out$Stage <- gsub("H33","Unnormal",nbt.H33.out$Stage)
nbt.H33.out$Embryo=sapply( strsplit(as.character(rownames(nbt.H33.out)), "_"), "[[", 4 )

nbt.H33.out$Location=sapply( strsplit(as.character(rownames(nbt.H33.out)), "_"), "[[", 5 )

nbt.H33.out$cell_type <-NA
table(nbt.H33.out$Embryo)
#mycolors=c("blue","brown","purple","darkgreen","RoyalBlue1","grey","DarkGoldenrod3","springgreen","deeppink1","firebrick1","cyan","deepskyblue4","pink","lightgoldenrod3","hotpink","steelblue1","PaleGreen","indianred3","lightblue","mediumpurple1")
mycolors <- c("pink","gray","aquamarine2","green","red","navy","purple","cyan","darkorange2","cornflowerblue","deeppink","orchid1","darkturquoise","deepskyblue4","dimgray","palegreen","khaki1",
              "darkolivegreen1","slateblue","yellow","darkred","darksalmon","blue","gold4","forestgreen")

pdf("../Figure/S3_H33_UMAP.pdf",width = 4.3,height = 3)
p<-ggplot(nbt.H33.out,aes(x=UMAP_1,y=UMAP_2,color= Stage,shape=Location))
p<-p+geom_point(size=1.1)+theme_bw()+
  theme(panel.grid=element_blank())+scale_color_manual(values = c("S3.MTE"="#BB47F5","S3.PTE"="#ADE7EC","S3.EPI"="#bf9a8e","S3.PE"="darkgoldenrod1","S3.ISH"="darkgreen","Unnormal"="deeppink"),limits=c("S3.MTE","S3.PTE","S3.EPI","S3.PE","S3.ISH","Unnormal"))
p
dev.off()

pdf("../Figure/S3_H33_embryo_UMAP.pdf",width = 4,height = 3)
ggplot(nbt.H33.out,aes(x=UMAP_1,y=UMAP_2,color= Embryo,shape=Location))+geom_point(size=1.1)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("dimgray","palegreen","khaki1","darkolivegreen1","blue","black"),limits=c("E15","E16","E17","E18","E23","E33"))
dev.off()




new.marker.plot <- function(x){
  nbt.H33.out[ ,x]=t(S3_tpm.new_H33.log[x,])
  gene=nbt.H33.out[ ,x]
  ggplot(nbt.H33.out,aes(x=UMAP_1,y=UMAP_2,color=gene))+
    geom_point(size=1)+theme(axis.line = element_line(size = 0.5, colour = "black"),
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

CDX2 <- new.marker.plot("CDX2")
KRT18 <- new.marker.plot("KRT18")
GATA2 <- new.marker.plot("GATA2")
GATA3 <- new.marker.plot("GATA3")
CLDN10 <- new.marker.plot("CLDN10")
CLDN4 <- new.marker.plot("CLDN4")
FUT4 <- new.marker.plot("FUT4")
t1 <- new.marker.plot("FOSL1")
t2 <- new.marker.plot("TPRN")
t3 <- new.marker.plot("WNT7B")
t4 <- new.marker.plot("C11orf95")
t5 <- new.marker.plot("IL6")
t6 <- new.marker.plot("ZBTB7A")
t7 <- new.marker.plot("RAB17")
t8 <- new.marker.plot("SORD")
t9 <- new.marker.plot("CEMP1")

t10 <- new.marker.plot("QDPR")
t11 <- new.marker.plot("MRGPRX1")
t12 <- new.marker.plot("BCS1L")
t13 <- new.marker.plot("COX7B2")
t14 <- new.marker.plot("CD53")
t15 <- new.marker.plot("CFD")
t16 <- new.marker.plot("PRPSAP1")
t17 <- new.marker.plot("HDDC3")
t18 <- new.marker.plot("BCKDK")
t19 <- new.marker.plot("GSTO1")
t20 <- new.marker.plot("STAG3")
ggplot2.multiplot(t1,t2,t3,t4,t5,t6,t7,t8,t9,cols=4)
ggplot2.multiplot(t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,cols=4)
library(easyGgplot2)
pdf("../Figure/S3_H33_TE_markers.pdf",width = 13,height = 6)
ggplot2.multiplot(CDX2,KRT18,GATA2,GATA3,CLDN10,FUT4,cols=3)
dev.off()

CLDN3 <- new.marker.plot("CLDN3")
FGFR4 <- new.marker.plot("FGFR4")
new.marker.plot("HAND1")
new.marker.plot("PGF")
new.marker.plot("GCM1")
new.marker.plot("SDC1")
new.marker.plot("IDH1")
new.marker.plot("NDUFAB1")
SOX2 <- new.marker.plot("SOX2")
POU5F1 <- new.marker.plot("POU5F1")
NANOG <- new.marker.plot("NANOG")


pdf("../Figure/S3_H33_ICM_markers.pdf",width = 13,height = 3)
ggplot2.multiplot(SOX2,NANOG,POU5F1,cols=3)
dev.off()

GATA4 <- new.marker.plot("GATA4")
PDGFRA <- new.marker.plot("PDGFRA")
FOXA2 <- new.marker.plot("FOXA2")
pdf("../Figure/S3_H33_PE_markers.pdf",width = 13,height = 3)
ggplot2.multiplot(GATA4,PDGFRA,FOXA2,cols=3)
dev.off()



HAPLN1 <- new.marker.plot("HAPLN1")
ENPEP <- new.marker.plot("ENPEP")
HAND1 <- new.marker.plot("HAND1")
CD53 <- new.marker.plot("CD53")
CDX1 <- new.marker.plot("CDX1")

pdf("../Figure/S3_H33_MTE_markers.pdf",width = 13,height = 6)
ggplot2.multiplot(HAPLN1,ENPEP,HAND1,CD53,CDX1,cols=3)
dev.off()

CLMP <- new.marker.plot("CLMP")
ADGRL3 <- new.marker.plot("ADGRL3")
PID1 <- new.marker.plot("PID1")
IL1R1 <- new.marker.plot("IL1R1")
FST <- new.marker.plot("FST")
pdf("../Figure/S3_H33_PTE_markers.pdf",width = 13,height = 6)
ggplot2.multiplot(CLMP,ADGRL3,PID1,IL1R1,FST,cols=3)
dev.off()

ENG <- new.marker.plot("ENG")
MAGEA3 <- new.marker.plot("MAGEA3")
FUT3 <- new.marker.plot("FUT3")
METTL7A <- new.marker.plot("METTL7A")
CDKN2A <- new.marker.plot("CDKN2A")
pdf("../Figure/S3_H33_ISH_markers.pdf",width = 13,height = 6)
ggplot2.multiplot(ENG,MAGEA3,FUT3,METTL7A,CDKN2A,cols=3)
dev.off()

########################FIgure 6c
nbt.H33.out[ ,"HAND1"]=t(S3_tpm.new_H33.log["HAND1",])
HAND1 <- ggplot(nbt.H33.out,aes(x=UMAP_1,y=UMAP_2,color=HAND1))+
  geom_point(size=1)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                           axis.ticks.length = unit(0.15, "cm"),
                           panel.background = element_rect(fill = "white", colour = "white"),
                           axis.title.y= element_text(size=9),
                           title= element_text(size=9),
                           legend.text= element_text(size=9),
                           axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.0,name= "log2(TPM/10+1)")+
  labs(title="HAND1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))

nbt.H33.out.H33 <- nbt.H33.out[which(nbt.H33.out$Embryo %in% c("E33") & (nbt.H33.out$UMAP_1 >5) & (nbt.H33.out$UMAP_2 < -6)), ]
library(ggrepel)
p<-ggplot(nbt.H33.out.H33)
p<-p+geom_point(aes(UMAP_1,UMAP_2,color=Embryo,shape=Location))+ geom_text_repel(aes(UMAP_1,UMAP_2,label=rownames(nbt.H33.out.H33)))+theme+xlim(-30,30)#+scale_color_manual(values=mycolors)
p

nbt.H33.out[ ,"GCM1"]=t(S3_tpm.new_H33.log["GCM1",])
GCM1 <- ggplot(nbt.H33.out,aes(x=UMAP_1,y=UMAP_2,color=GCM1))+
  geom_point(size=1)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                           axis.ticks.length = unit(0.15, "cm"),
                           panel.background = element_rect(fill = "white", colour = "white"),
                           axis.title.y= element_text(size=9),
                           title= element_text(size=9),
                           legend.text= element_text(size=9),
                           axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.9,name= "log2(TPM/10+1)")+
  labs(title="GCM1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
library(easyGgplot2)
pdf("../Figure/Figure4_S3_H33_TE_markers.pdf",width = 6,height = 3)
ggplot2.multiplot(HAND1,GCM1,cols=2)
dev.off()

ggplot(nbt.H33.out,aes(x=UMAP_1,y=UMAP_2,color=HAND1))+
  geom_point(size=1)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                           axis.ticks.length = unit(0.15, "cm"),
                           panel.background = element_rect(fill = "white", colour = "white"),
                           axis.title.y= element_text(size=9),
                           title= element_text(size=9),
                           legend.text= element_text(size=9),
                           axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.0,name= "log2(TPM/10+1)")+
  labs(title="HAND1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="right")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))
ggsave("../Figure/Figure4_S3_H33_TE_markers_label.pdf")

#######################MTE
#------------------normal--------------
MTE_c0_name <- data.frame(name=WhichCells(nbt.H33,0))
class1 <- sapply( strsplit(as.character(MTE_c0_name$name), "_"), "[[", 5 )
MTE_c0_name.M <- MTE_c0_name[which(class1 %in% c("M")), ]
MTE_c0_tpm <- S3_tpm.new_H33[ ,which(colnames(S3_tpm.new_H33) %in% MTE_c0_name.M)]
colnames(MTE_c0_tpm)
#------------------Unnormal--------------
MTE_c3_name <- data.frame(name=WhichCells(nbt.H33,3))
clasS3 <- sapply( strsplit(as.character(MTE_c3_name$name), "_"), "[[", 5 )
MTE_c3_name.M <- MTE_c3_name[which(clasS3 %in% c("M")), ]
MTE_c3_tpm <- S3_tpm.new_H33[ ,which(colnames(S3_tpm.new_H33) %in% MTE_c3_name.M[12:82])]
colnames(MTE_c3_tpm)
normal_unnormal_MTE_tpm <- cbind(MTE_c0_tpm,MTE_c3_tpm)
normal_unnormal_MTE_tpm <- normal_unnormal_MTE_tpm[ ,-which(colnames(normal_unnormal_MTE_tpm) %in% c("S3.PTE_S3_hTE_E18_M_17","S3.PE_S3_hTE_E23_M_2"))]
normal_unnormal_MTE_tpm.log <- log2(normal_unnormal_MTE_tpm/10+1)
######################FIGURE6E
SIG.gene.chart <- data.frame(sample=colnames(normal_unnormal_MTE_tpm.log),BCS1L=t(normal_unnormal_MTE_tpm.log["BCS1L", ]),COX7B2=t(normal_unnormal_MTE_tpm.log["COX7B2", ]),GSTO1=t(normal_unnormal_MTE_tpm.log["GSTO1", ]),PRPSAP1=t(normal_unnormal_MTE_tpm.log["PRPSAP1", ]),CD53=t(normal_unnormal_MTE_tpm.log["CD53", ]),OXA1L=t(normal_unnormal_MTE_tpm.log["OXA1L", ]),GAS6=t(normal_unnormal_MTE_tpm.log["GAS6", ]),RPS38=t(normal_unnormal_MTE_tpm.log["RPS38", ]),RPS39=t(normal_unnormal_MTE_tpm.log["RPS39", ]),RPLP1=t(normal_unnormal_MTE_tpm.log["RPLP1", ]),Group=sapply( strsplit(as.character(colnames(normal_unnormal_MTE_tpm.log)), "_"), "[[", 1 ))
SIG.gene.chart$Group <- factor(SIG.gene.chart$Group,levels = c("S3.MTE","H33"))
library(ggplot2)
library(Rmisc) 
tgc1 <- summarySE(SIG.gene.chart, measurevar="BCS1L", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc1$Group <- factor(tgc1$Group)
p1 <- ggplot(tgc1, aes(x=Group, y=BCS1L, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=BCS1L-se, ymax=BCS1L+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+theme(panel.grid=element_blank())+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc2 <- summarySE(SIG.gene.chart, measurevar="COX7B2", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc2$Group <- factor(tgc2$Group)
p2 <- ggplot(tgc2, aes(x=Group, y=COX7B2, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=COX7B2-se, ymax=COX7B2+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc3 <- summarySE(SIG.gene.chart, measurevar="GSTO1", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc3$Group <- factor(tgc3$Group)
p3 <- ggplot(tgc3, aes(x=Group, y=GSTO1, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=GSTO1-se, ymax=GSTO1+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc4 <- summarySE(SIG.gene.chart, measurevar="PRPSAP1", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc4$Group <- factor(tgc4$Group)
p4 <- ggplot(tgc4, aes(x=Group, y=PRPSAP1, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=PRPSAP1-se, ymax=PRPSAP1+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc5 <- summarySE(SIG.gene.chart, measurevar="CD53", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc5$Group <- factor(tgc5$Group)
p5 <- ggplot(tgc5, aes(x=Group, y=CD53, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=CD53-se, ymax=CD53+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

library(easyGgplot2)
pdf("G:/Qiao lab/liudandan/human/Figure/H33_downgene.pdf",width = 5.4,height = 2)
ggplot2.multiplot(p1,p2,p3,p4,p5,cols=5)
dev.off()



tgc6 <- summarySE(SIG.gene.chart, measurevar="OXA1L", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc6$Group <- factor(tgc6$Group)
p6 <- ggplot(tgc6, aes(x=Group, y=OXA1L, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=OXA1L-se, ymax=OXA1L+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc7 <- summarySE(SIG.gene.chart, measurevar="GAS6", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc7$Group <- factor(tgc7$Group)
p7 <- ggplot(tgc7, aes(x=Group, y=GAS6, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=GAS6-se, ymax=GAS6+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc8 <- summarySE(SIG.gene.chart, measurevar="RPS38", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc8$Group <- factor(tgc8$Group)
p8 <- ggplot(tgc8, aes(x=Group, y=RPS38, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=RPS38-se, ymax=RPS38+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc9 <- summarySE(SIG.gene.chart, measurevar="RPS39", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc9$Group <- factor(tgc9$Group)
p9 <- ggplot(tgc9, aes(x=Group, y=RPS39, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=RPS39-se, ymax=RPS39+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))


tgc10 <- summarySE(SIG.gene.chart, measurevar="RPLP1", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc10$Group <- factor(tgc10$Group)
p10 <- ggplot(tgc10, aes(x=Group, y=RPLP1, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=RPLP1-se, ymax=RPLP1+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#bb47f5","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

library(easyGgplot2)
pdf("G:/Qiao lab/liudandan/human/Figure/PGD_H33_gene.pdf",width = 5.4,height = 2.2)
ggplot2.multiplot(p1,p6,p7,p8,cols=4)
dev.off()


ggplot2.multiplot(p8,p9,p10,cols=3)


############################################
library(Seurat)
library(dplyr)
nbt.H33.MTE <- CreateSeuratObject(raw.data =normal_unnormal_MTE_tpm, min.cells = 3, min.genes = 2000,names.field = 1,names.delim = "_")
nbt.H33.MTE <- NormalizeData(object = nbt.H33.MTE, normalization.method = "LogNormalize",scale.factor = 10000)
nbt.H33.MTE <- FindVariableGenes(object = nbt.H33.MTE, mean.function = ExpMean, dispersion.function = LogVMR,
                                 x.low.cutoff = 0.4, y.cutoff = 0.4)
nbt.H33.MTE.all.markers=FindAllMarkers(nbt.H33.MTE,test.use = "roc")
#nbt.H33.MTE.markers.use=subset(nbt.H33.MTE.all.markers,avg_diff>0&power>0.4)
write.table(nbt.H33.MTE.markers.use,"../Table/H33_MTE_markers.xls",quote = F,col.names = T,row.names = F,sep = "\t")
#nbt.H33.MTE.markers.use1 <- nbt.H33.MTE.markers.use[which(nbt.H33.MTE.markers.use$pct.1 > 0.8 & nbt.H33.MTE.markers.use$pct.2 < 0.5), ]
nbt.abnormal.MTE.markers.use.up=subset(nbt.H33.MTE.all.markers,avg_diff>0&power>0.4)
nbt.abnormal.MTE.markers.use.up.abnormal=nbt.abnormal.MTE.markers.use.up[which(nbt.abnormal.MTE.markers.use.up$cluster %in% c("H33")), ]

nbt.abnormal.MTE.markers.use.down=subset(nbt.H33.MTE.all.markers,avg_diff<0&power>0.4)
nbt.abnormal.MTE.markers.use.down.abnormal=nbt.abnormal.MTE.markers.use.down[which(nbt.abnormal.MTE.markers.use.down$cluster %in% c("H33")), ]

nbt.abnormal.MTE.markers.use.abnormal <- rbind(nbt.abnormal.MTE.markers.use.up.abnormal,nbt.abnormal.MTE.markers.use.down.abnormal)
write.table(nbt.abnormal.MTE.markers.use.abnormal,"../Table/E26_MTE_markers.xls",quote = F,col.names = T,row.names = F,sep = "\t")

SEQ <- as.vector(rownames(nbt.abnormal.MTE.markers.use.abnormal))
normal_unnormal_MTE_tpm_DEG <- normal_unnormal_MTE_tpm[SEQ, ]
mean <- rowMeans(normal_unnormal_MTE_tpm_DEG)
sd_value <- apply(normal_unnormal_MTE_tpm_DEG,1,sd)

z_score <- (normal_unnormal_MTE_tpm_DEG-mean)/sd_value
z_score[z_score>=1] <- 1
z_score[z_score<= -1] <- -1
subj <- sapply( strsplit(as.character(colnames(z_score)), "\t"), "[[", 1 )
aaka2 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1))
#aaka21 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 3))
row.names(aaka2)<-subj
#colnames(aaka2)<-c("Time")
#row.names(aaka21)<-subj
#colnames(aaka21)<-c("Time")
#aaka22 = cbind(aaka2,aaka21)

#aka5=list(Time=c("Early"="firebrick","Middle"="darkgreen","Late"="navy"))#Type=c("GM"="magenta",'RI'="olivedrab1","CM"="skyblue1","LH"="orange","IMM"="aquamarine","SOX4"="seagreen4","Blood"="mediumpurple","PD"="deepskyblue","ED"="tomato","BM"="turquoise1","CD"="khaki3")"firebrick","darkgreen","navy"
library(RColorBrewer)
library("pheatmap")
pheatmap( z_score,annotation_col =aaka2,#annotation_colors = aka5,gaps_row=(6),
          #color =  colorRampPalette(c("dodgerblue1", "black", "firebrick1"))(100),
          #color = colorRampPalette(c("cyan", "black", "firebrick1"))(100), 
          color = colorRampPalette(rev(brewer.pal(n = 5, name ="RdBu")))(50),
          #breaks = c(seq(-2,-1.1,length=25),seq(-1,1,length=50),seq(1.1,2,length=25)),
          cluster_row = F,cluster_col =F,show_rownames = F, 
          show_colnames =F, scale = "none",legend =T,border_color =T,
          #cellwidth = 0.6,cellheight = 8,fontsize_row=8,
          fontsize=8,fontsize_col = 5)
#######################PTE
PTE_c1_name <- data.frame(name=WhichCells(nbt.H33,1))
PTE_tpm <- S3_tpm.new_H33[ ,which(colnames(S3_tpm.new_H33) %in% PTE_c1_name$name[5:117])]
colnames(PTE_tpm)
PTE_class <- sapply( strsplit(as.character(colnames(PTE_tpm)), "_"), "[[", 5)
PTE_tpm.use <- PTE_tpm[ ,which(PTE_class %in% c("P"))]
colnames(PTE_tpm.use)
PTE_tpm.use.log <- log2(PTE_tpm.use/10+1)
library(Seurat)
library(dplyr)
nbt.H33.PTE <- CreateSeuratObject(raw.data =PTE_tpm, min.cells = 3, min.genes = 2000,names.field = 1,names.delim = "_")
nbt.H33.PTE <- NormalizeData(object = nbt.H33.PTE, normalization.method = "LogNormalize",scale.factor = 10000)
nbt.H33.PTE <- FindVariableGenes(object = nbt.H33.PTE, mean.function = ExpMean, dispersion.function = LogVMR,
                                 x.low.cutoff = 0.4, y.cutoff = 0.4)
nbt.H33.PTE.all.markers=FindAllMarkers(nbt.H33.PTE,test.use = "roc")
nbt.H33.PTE.markers.use=subset(nbt.H33.PTE.all.markers,avg_diff>0&power>0.4)
write.table(nbt.H33.PTE.markers.use,"../Table/H33_PTE_markers.xls",quote = F,col.names = T,row.names = F,sep = "\t")


######################FIGURE s7b
SIG.gene.chart1 <- data.frame(sample=colnames(PTE_tpm.use.log),BCS1L=t(PTE_tpm.use.log["BCS1L", ]),COX7B2=t(PTE_tpm.use.log["COX7B2", ]),GSTO1=t(PTE_tpm.use.log["GSTO1", ]),PRPSAP1=t(PTE_tpm.use.log["PRPSAP1", ]),CD53=t(PTE_tpm.use.log["CD53", ]),OXA1L=t(PTE_tpm.use.log["OXA1L", ]),GAS6=t(PTE_tpm.use.log["GAS6", ]),RPS38=t(PTE_tpm.use.log["RPS38", ]),Group=sapply( strsplit(as.character(colnames(PTE_tpm.use.log)), "_"), "[[", 1 ))
SIG.gene.chart1$Group <- factor(SIG.gene.chart1$Group,levels = c("S3.PTE","H33"))
library(ggplot2)
library(Rmisc) 
tgc1 <- summarySE(SIG.gene.chart1, measurevar="BCS1L", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc1$Group <- factor(tgc1$Group)
p1 <- ggplot(tgc1, aes(x=Group, y=BCS1L, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=BCS1L-se, ymax=BCS1L+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+theme(panel.grid=element_blank())+scale_fill_manual(values = c("#ADE7EC","#ff1493"))+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc2 <- summarySE(SIG.gene.chart1, measurevar="COX7B2", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc2$Group <- factor(tgc2$Group)
p2 <- ggplot(tgc2, aes(x=Group, y=COX7B2, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=COX7B2-se, ymax=COX7B2+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#ADE7EC","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc3 <- summarySE(SIG.gene.chart1, measurevar="GSTO1", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc3$Group <- factor(tgc3$Group)
p3 <- ggplot(tgc3, aes(x=Group, y=GSTO1, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=GSTO1-se, ymax=GSTO1+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#ADE7EC","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc4 <- summarySE(SIG.gene.chart1, measurevar="PRPSAP1", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc4$Group <- factor(tgc4$Group)
p4 <- ggplot(tgc4, aes(x=Group, y=PRPSAP1, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=PRPSAP1-se, ymax=PRPSAP1+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#ADE7EC","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

tgc5 <- summarySE(SIG.gene.chart1, measurevar="CD53", groupvars=c("Group"))
# Error bars represent standard error of the mean
tgc5$Group <- factor(tgc5$Group)
p5 <- ggplot(tgc5, aes(x=Group, y=CD53, fill=Group)) + 
  geom_bar(position=position_dodge(), stat="identity",width=0.5) +
  geom_errorbar(aes(ymin=CD53-se, ymax=CD53+se),
                width=.2, # 设置误差线的宽度 
                position=position_dodge(.9))+theme_bw()+scale_fill_manual(values = c("#ADE7EC","#ff1493"))+theme(panel.grid=element_blank())+theme(title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),axis.title.x= element_text(size=9,colour = "white"),legend.position = "none",  axis.text.x = element_text(size=9,colour="white",angle=45))

library(easyGgplot2)
pdf("G:/Qiao lab/liudandan/human/Figure/H33.PTE_downgene.pdf",width = 5.4,height = 2)
ggplot2.multiplot(p1,p2,p3,p4,p5,cols=5)
dev.off()

##############unused
###########MERGE S3.ISH & ISK
isk_ish_tpm <- cbind(C2_tpm,ISK_TPM)
isk_ish_tpm.log <- log2(isk_ish_tpm/10+1)
ligand_Data <- isk_ish_tpm.log[c("ENG","MAGEA3","FUT3","METTL7A","CDKN2A","FBLN2","ECM1","PTGDS","TNNT3","PHLDA1","PTPRS"), ]

sample <- data.frame(sample=colnames(ligand_Data),Type=sapply( strsplit(as.character(colnames(ligand_Data)), "_"), "[[", 1))
head(sample)               
library(ggplot2)
library(reshape2)
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
pdf(paste(week,name,"ligand_vio.pdf",sep = "_"),width = 2.8,height = 2*nrow(ligand_Data))
tmp_Exp <-cbind(as.data.frame( t(ligand_Data)),sample)
melt_data <- melt(tmp_Exp,id.vars=c("sample","Type"),variable.name = "Gene",value.name = "Exp")

head(melt_data)
p <- ggplot(data=melt_data, aes(x=Type,y=Exp))+geom_violin(aes(fill =Type))+theme_bw()+stat_summary(fun.y = median,geom = "point",fill="black",shape=21,size=1.5)+#geom_line()+
  theme(axis.line = element_line(size = 0.5, colour = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.title.y= element_text(size=9,colour = "white"),
        title= element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5),
        legend.title=element_text(size=9, face="bold"),
        legend.text= element_text(size=9),
        axis.title.x= element_text(size=9,colour = "white"),
        axis.text.x = element_text(size=9,colour="white"),
        strip.text.x = element_text(face="italic",size=9),strip.background = element_blank())+
  scale_fill_manual(values = c("#9660A8","#CFCEC8"))+
  #scale_color_manual(values = c("purple","pink","springgreen"))+
  theme(panel.grid=element_blank())+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))+facet_wrap(~Gene,ncol=4,scales="free")

print(p)

######################UNused
nbt.H33_C0_C3.all.markers=FindMarkers(nbt.H33,ident.1 = 3,ident.2 = 0,test.use = "roc")
UP.nbt.H33.out_C0_C3=subset(nbt.H33_C0_C3.all.markers,avg_diff>0&power>0.4)
UP.nbt.H33.out_C0_C3$Gene <- rownames(UP.nbt.H33.out_C0_C3)
DOWN.nbt.H33.out_C0_C3=subset(nbt.H33_C0_C3.all.markers,avg_diff<0&power>0.4)
DOWN.nbt.H33.out_C0_C3$Gene <- rownames(DOWN.nbt.H33.out_C0_C3)
write.table(UP.nbt.H33.out_C0_C3,"../Table/H33_up_markers.xls",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(DOWN.nbt.H33.out_C0_C3,"../Table/H33_down_markers.xls",quote = F,col.names = T,row.names = F,sep = "\t")

