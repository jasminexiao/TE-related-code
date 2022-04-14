setwd("E:/Qiaolab/liudandan/human/ALL_data")
TPM <- read.table("qc_tpm(1127).txt",as.is = TRUE)
#TPM.log <- log2(TPM/10+1)
Stage.class <- sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 1 )
S1.tpm <- TPM[ ,which(Stage.class %in% c("S1"))]
S1.tpm.log <- log2(S1.tpm/10+1)
#Embryo.class <- sapply( strsplit(as.character(colnames(S3.tpm)), "_"), "[[", 3 )
#S2.tpm.no.E25 <- S2.tpm[ ,-which(Embryo.class %in% c("E25"))]
#S2.tpm.no.E25.log <- log2(S2.tpm.no.E25/10+1)
library(Seurat)
library(dplyr)
nbt <- CreateSeuratObject(S1.tpm, project = "S1", min.cells = 3, min.features = 200)
nbt
nbt[["percent.mt"]] <- PercentageFeatureSet(nbt, pattern = "^MT-")
VlnPlot(nbt, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
nbt <- NormalizeData(nbt, normalization.method = "LogNormalize", scale.factor = 10000)
nbt <- NormalizeData(nbt)
nbt <- FindVariableFeatures(nbt, selection.method = "vst", nfeatures = 2000)

# plot variable features with and without labels
VariableFeaturePlot(nbt)
nbt <- ScaleData(nbt)
nbt <- RunPCA(nbt, features = VariableFeatures(object = nbt))
DimPlot(nbt, reduction = "pca")
DimHeatmap(nbt, dims = 1:15, cells = 500, balanced = TRUE)
nbt <- JackStraw(nbt, num.replicate = 100)
nbt <- ScoreJackStraw(nbt, dims = 1:20)
JackStrawPlot(nbt, dims = 1:20)
ElbowPlot(nbt)
nbt <- FindNeighbors(nbt, dims = 1:11)
nbt <- FindClusters(nbt, resolution = 0.2)
nbt <- RunUMAP(nbt, dims = 1:11)
DimPlot(nbt, reduction = "umap",label = T)


nbt.out=data.frame(nbt@reductions$umap@cell.embeddings)
nbt.out$SAMPLE=rownames(nbt.out)
nbt.out$Stage=sapply( strsplit(as.character(rownames(nbt.out)), "_"), "[[", 1 )
nbt.out$Embryo=sapply( strsplit(as.character(rownames(nbt.out)), "_"), "[[", 3 )

nbt.out$Location=sapply( strsplit(as.character(rownames(nbt.out)), "_"), "[[", 4 )

nbt.out$cell_type <-NA
table(nbt.out$Embryo)
#mycolors=c("blue","brown","purple","darkgreen","RoyalBlue1","grey","DarkGoldenrod3","springgreen","deeppink1","firebrick1","cyan","deepskyblue4","pink","lightgoldenrod3","hotpink","steelblue1","PaleGreen","indianred3","lightblue","mediumpurple1")
mycolors <- c("pink","gray","aquamarine2","green","red","cyan","purple","navy","darkorange2","cornflowerblue","deeppink","darkorchid1","darkturquoise","deepskyblue4",
              "dimgray","greenyellow","khaki1","darkolivegreen1",
              "black","darkorchid1","darkred","darksalmon","darkslateblue","darkturquoise","deepskyblue",
              "forestgreen")#,"gold4","black","maroon","mediumturquoise"
#nbt.out$Embryo <- factor(nbt.out$Embryo,levels = c("E6","E9","E14","E19","E20","E21","E22","E25"))
library(ggplot2)
p<-ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= Embryo,shape=Location))
p<-p+geom_point(size=1.5)+theme_bw()+
  scale_color_manual(values = c("#96C4E6","#194D5B","#6C15AB","#FAAE4E","#5579ED","#9C6628","#E560CD","#24B300"))+
  theme(panel.grid=element_blank())+scale_shape_manual(values=c(17,15),limits=c("M","P"))#+scale_x_continuous(limits = c(-40,40))+scale_y_continuous(limits = c(-40,40))+scale_color_manual(values = c("cyan","darkorange2","deepskyblue4","black","darkorchid1","darkred","darksalmon","forestgreen"),limits=c("E6","E9","E14","E19","E20","E21","E22","E25"))
pdf("../Figure/Revise_S1_UMAP_embryo.pdf",width = 5,height = 3.5)
p
dev.off()
S1.tpm.log1 <- S1.tpm.log[ ,as.vector(rownames(nbt.out))]
marker.plot <- function(x){
  nbt.out[ ,x]=t(S1.tpm.log1[x,])
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
pdf("../Figure/Figure4a.1.pdf",width = 13,height = 6)
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


pdf("../Figure/Figure4a.2.pdf",width = 13,height = 3)
ggplot2.multiplot(SOX2,NANOG,POU5F1,cols=3)
dev.off()

GATA4 <- marker.plot("GATA4")
PDGFRA <- marker.plot("PDGFRA")
FOXA2 <- marker.plot("FOXA2")
pdf("../Figure/Figure4a.3.pdf",width = 13,height = 3)
ggplot2.multiplot(GATA4,PDGFRA,FOXA2,cols=3)
dev.off()


pdf("../Figure/Figure4a.4.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= Embryo,shape=Location))+geom_point(size=1)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("cyan","darkorange2","deepskyblue4","black","darkorchid1","darkred","darksalmon","forestgreen"),limits=c("E6","E9","E14","E19","E20","E21","E22","E25"))+scale_shape_manual(values=c(17,15),limits=c("M","P"))
dev.off()

nbt.out[which(nbt.out$UMAP_2< -3), ]$cell_type=c("S1.ICM")
nbt.out[-which(nbt.out$UMAP_2< -3), ]$cell_type=c("S1.TE")

pdf("../Figure/revise_Human_S1_umap.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3.5)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("#ed4e1d","#92F725"))+scale_shape_manual(values=c(17,15),limits=c("M","P"))
dev.off()
write.table(nbt.out,"../Table/Figure2a_S1_coordinate.xls",quote = F,col.names = T,row.names = F,sep = "\t")
#############FIG S4a
TE_chart <- data.frame(CDX2=t(S1.tpm.log1["CDX2",]))
TE_chart[ ,"GATA2"]=t(S1.tpm.log1["GATA2",])
TE_chart[ ,"GATA3"]=t(S1.tpm.log1["GATA3",])
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

PE_chart <- data.frame(GATA4=t(S1.tpm.log1["GATA4",]))
PE_chart[ ,"PDGFRA"]=t(S1.tpm.log1["PDGFRA",])
PE_chart[ ,"FOXA2"]=t(S1.tpm.log1["FOXA2",])
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
  scale_colour_gradient2(low="blue",high="red",midpoint = 1,name= "log2(TPM/10+1)")+
  labs(title="PE")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

EPI_chart <- data.frame(SOX2=t(S1.tpm.log1["SOX2",]))
EPI_chart[ ,"NANOG"]=t(S1.tpm.log1["NANOG",])
EPI_chart[ ,"POU5F1"]=t(S1.tpm.log1["POU5F1",])
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
pdf("../Figure/revise_Human_S1_markers_exp.pdf",width = 6,height = 2)
ggplot2.multiplot(TE,EPI,PE,cols=3)
dev.off()


####################autocluster
nbt.cluster <- FindClusters(nbt, reduction.type = "pca", dims.use = 1:7, resolution = 0.1, print.output = F)
UMAPPlot(nbt.cluster)

##########ICM subtype
S1_tpm.select <- S1.tpm[c("SOX2","NANOG","POU5F1","GATA4","PDGFRA","FOXA2","CDX2","KRT18","GATA2","GATA3","CLDN10","FUT4"), ]
S1_tpm.select.log <- log2(S1_tpm.select/10+1)

library('RColorBrewer')
library("pheatmap")
pheatmap( S1_tpm.select.log,
          color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(50), 
          #breaks = c(seq(-2,-1.1,length=15),seq(-1,1,length=70),seq(1.1,2,length=15)),
          cluster_row = T,cluster_col =T,show_rownames = T, 
          show_colnames =F, scale = "none",legend =T,border_color =T,
          fontsize_row=8,
          fontsize=8,fontsize_col = 5
)

dst2_euc<-dist(t(S1_tpm.select.log))
hc<-hclust(dst2_euc,"complete")
plot(hc,hang=-1,col="darkgreen",cex=0.3)
rect.hclust(hc, k=2, border="red")

cl_2 <- rect.hclust(hc, k=2, which = c(1,2), border=c('red','green'))
names(cl_2[[1]])

nbt.out[which(nbt.out$SAMPLE %in% (names(cl_2[[1]]))), ]$cell_type <- c("S1.EPI")
C2_tpm.PE <- S1.tpm[ ,which(colnames(S1.tpm) %in% nbt.out[(nbt.out$cell_type %in% c("S1.EPI")), ]$SAMPLE)]
colnames(C2_tpm.PE) <- paste("S1.EPI",colnames(C2_tpm.PE),sep = "_")
write.table(C2_tpm.PE,"S1.EPI_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")

nbt.out[which(nbt.out$SAMPLE %in% (names(cl_2[[2]]))), ]$cell_type <- c("S1.TE")
C2_tpm.EPI <- S1.tpm[ ,which(colnames(S1.tpm) %in% nbt.out[(nbt.out$cell_type %in% c("S1.TE")), ]$SAMPLE)]
colnames(C2_tpm.EPI) <- paste("S1.TE",colnames(C2_tpm.EPI),sep = "_")
write.table(C2_tpm.EPI,"S1.TE_tpm.txt",quote = F,col.names = T,row.names = T,sep = "\t")
nbt.out$cell_type <- gsub("EPI","S1.ICM",nbt.out$cell_type)
pdf("../Figure/Figure4a.5.pdf",width = 6.5,height = 5)
ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color= cell_type,shape=Location))+geom_point(size=3.5)+theme_bw()+theme(panel.grid=element_blank())+scale_color_manual(values = c("#ed4e1d","#92F725"))+scale_shape_manual(values=c(17,15),limits=c("M","P"))
dev.off()
write.table(nbt.out,"Figure3a_S1_coordinate.xls",quote = F,col.names = T,row.names = F,sep = "\t")
#############S1.TE
S1.TE_name=data.frame(name=nbt.out[which(nbt.out$cell_type %in% c("S1.TE")), ]$SAMPLE)
S1.EPI_name=data.frame(name=nbt.out[which(nbt.out$cell_type %in% c("EPI")), ]$SAMPLE)
write.table(S1.TE_name,"S1.TE_name.txt",quote = F,col.names = T,row.names = F)
write.table(S1.EPI_name,"S1.EPI_name.txt",quote = F,col.names = T,row.names = F)

#################CELLPHONE
S1_TPM <- cbind(C2_tpm.PE,C2_tpm.EPI)
#ensambl_id <- read.table("G:/tang lab/huyuqiong/01.RetinaProject/cyd_analysis/ensembl_symbol_geneid.txt",header = T,sep = "\t")
ensambl_id <- read.table("G:/Qiao lab/liudandan/human/Table/cellphone/GeneSymbol_Ensembl_ID.txt",header = T,sep = "\t")
S1_TPM.sel <- S1_TPM[which(rownames(S1_TPM) %in% ensambl_id$gene_name), ]
seq <- as.vector(ensambl_id$gene_name)
S1_TPM.sel <- S1_TPM.sel[seq, ]
rownames(S1_TPM.sel) <- ensambl_id$ensembl
S1_TPM.sel.log <- log2(S1_TPM.sel/10+1)
S1_TPM.sel.log[1:10,1:3]
write.table(S1_TPM.sel.log,"S1_TPM.sel.Ensembl.txt",quote = F,col.names = T,row.names = T,sep="\t")

S1_chart <- data.frame(Cell=colnames(S1_TPM),cell_type=sapply( strsplit(as.character(colnames(S1_TPM)), "_"), "[[", 1))
write.table(S1_chart,"S1_chart.txt",quote = F,col.names = T,row.names = F,sep="\t")


library(mygene)
mygenes=c("TP53","AGTR1")

queryMany(c("TP53","AGTR1"), scopes="symbol", fields=c("ensembl.gene","entrezgene"), species="human")
out <- queryMany(rownames(S1_TPM), scopes="symbol", fields=c("ensembl.gene","entrezgene"), species="human")
colnames(out)

out.use <- out[ ,c(1,4)]

out.use2 <- data.frame(symbol=out.use@listData$query,ensembl=unlist(out.use@listData$ensembl))
symbol=data.frame(symbol=out.use@listData$query)
ensembl=data.frame(unlist(out.use@listData$ensembl))
class(ensembl)

data.list<-list()
for (i in 1:2){
  data.list[[i]]=out.use@listData$ensembl[i]$gene
}

test1 <- data.frame(out.use@listData$ensembl[[16]]$gene)
length(out.use@listData$ensembl[[1]]$gene)




data.list<-list()
for (i in 1:25370){
  data.list[[i]]=length(out.use@listData$ensembl[[i]]$gene)
}

data.list1 <- data.list[which(data.list[ ,3] %in% c("2","0")), ]

test1 <- unlist(data.list[lapply(data.list, length)==1])


write.table(print(out.use),"out_id.txt",quote = F,col.names = T,row.names = F,sep = "\t")
sum(out$ensembl)
mylist2 = out.use[-which(sapply(out.use$ensembl, is.null))]
out.use@listData
ensembl=out.use@listData$query

#########revise umap markers plot
nbt.out[ ,"ATP6V0D2"]=t(S1.tpm.log1["ATP6V0D2",])
ATP6V0D2 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=ATP6V0D2))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="ATP6V0D2")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"CYP26A1"]=t(S1.tpm.log1["CYP26A1",])
CYP26A1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=CYP26A1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="CYP26A1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"HAND1"]=t(S1.tpm.log1["HAND1",])
HAND1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=HAND1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.3,name= "log2(TPM/10+1)")+
  labs(title="HAND1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"SDC1"]=t(S1.tpm.log1["SDC1",])
SDC1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=SDC1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 1.3,name= "log2(TPM/10+1)")+
  labs(title="SDC1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"GCM1"]=t(S1.tpm.log1["GCM1",])
GCM1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=GCM1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 1,name= "log2(TPM/10+1)")+
  labs(title="GCM1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

nbt.out[ ,"PGF"]=t(S1.tpm.log1["PGF",])
PGF <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=PGF))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 2.7,name= "log2(TPM/10+1)")+
  labs(title="PGF")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))



nbt.out[ ,"ERVFRD_1"]=t(S1.tpm.log1["ERVFRD-1",])
ERVFRD_1 <- ggplot(nbt.out,aes(x=UMAP_1,y=UMAP_2,color=ERVFRD_1))+
  geom_point(size=0.35)+theme(axis.line = element_line(size = 0.5, colour = "black"),
                              axis.ticks.length = unit(0.15, "cm"),
                              panel.background = element_rect(fill = "white", colour = "white"),
                              axis.title.y= element_text(size=9),
                              title= element_text(size=9),
                              legend.text= element_text(size=9),
                              axis.text.x = element_text(size=9))+
  #scale_colour_gradient(high="red",low="grey",space="Lab",name= "log2(TPM/10+1)") + 
  scale_colour_gradient2(low="blue",high="red",midpoint = 0.8,name= "log2(TPM/10+1)")+
  labs(title="ERVFRD-1")+theme_bw()+theme(panel.grid=element_blank(),legend.position="none")+
  theme(plot.title = element_text(size=12, color="black", face="italic", vjust=0.5, hjust=0.5))+
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5))

pdf("../Figure/revise_Human_S1_newfindmarkers_exp.pdf",width = 6,height = 6)
ggplot2.multiplot(ATP6V0D2,CYP26A1,HAND1,SDC1,GCM1,PGF,ERVFRD_1,cols=3)
dev.off()
