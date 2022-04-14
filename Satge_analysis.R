setwd("F:/Qiaolab/liudandan/human/ALL_data")
summary <- read.table("qc_summary(1127).txt",header = T,sep = "\t")
summary$class <- paste(summary$Stage,summary$Embryo,summary$Location,sep="_")
table(summary$class)
TPM=read.table("qc_tpm(1127).txt",header = T,sep = "\t")
TPM.class1 <-sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 1 )
TPM.class2 <-sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 3 )
TPM.class3 <-sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 4 )
TPM.class <- paste(TPM.class1,TPM.class2,TPM.class3,sep="_")
table(TPM.class)
TPM.S1.E1.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E14_M"), ])$Name)]
TPM.S1.E1.M.mean <- data.frame(S1_hTE_E14_M=rowMeans(TPM.S1.E1.M))

TPM.S1.E13.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E14_P"), ])$Name)]
TPM.S1.E13.M.mean <- data.frame(S1_hTE_E14_P=rowMeans(TPM.S1.E13.M))

TPM.S1.E13.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E19_P"), ])$Name)]
TPM.S1.E13.P.mean <- data.frame(S1_hTE_E19_P=rowMeans(TPM.S1.E13.P))

TPM.S1.E18.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E19_M"), ])$Name)]
TPM.S1.E18.M.mean <- data.frame(S1_hTE_E19_M=rowMeans(TPM.S1.E18.M))

TPM.S1.E18.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E20_M"), ])$Name)]
TPM.S1.E18.P.mean <- data.frame(S1_hTE_E20_M=rowMeans(TPM.S1.E18.P))

TPM.S1.E23.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E20_P"), ])$Name)]
TPM.S1.E23.M.mean <- data.frame(S1_hTE_E20_P=rowMeans(TPM.S1.E23.M))

TPM.S1.E23.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E21_M"), ])$Name)]
TPM.S1.E23.P.mean <- data.frame(S1_hTE_E21_M=rowMeans(TPM.S1.E23.P))

TPM.S1.E24.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E22_M"), ])$Name)]
TPM.S1.E24.M.mean <- data.frame(S1_hTE_E22_M=rowMeans(TPM.S1.E24.M))

TPM.S1.E24.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E22_P"), ])$Name)]
TPM.S1.E24.P.mean <- data.frame(S1_hTE_E22_P=rowMeans(TPM.S1.E24.P))

#TPM.S1.E25.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E22_P"), ])$name)]
#TPM.S1.E25.M.mean <- data.frame(S1_hTE_E25_M=rowMeans(TPM.S1.E25.M))

TPM.S1.E26.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E25_M"), ])$Name)]
TPM.S1.E26.M.mean <- data.frame(S1_hTE_E25_M=rowMeans(TPM.S1.E26.M))

TPM.S1.E26.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E25_P"), ])$Name)]
TPM.S1.E26.P.mean <- data.frame(S1_hTE_E25_P=rowMeans(TPM.S1.E26.P))

TPM.S1.E30.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E9_M"), ])$Name)]
TPM.S1.E30.M.mean <- data.frame(S1_hTE_E9_M=rowMeans(TPM.S1.E30.M))

TPM.S1.E30.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E9_P"), ])$Name)]
TPM.S1.E30.P.mean <- data.frame(S1_hTE_E9_P=rowMeans(TPM.S1.E30.P))

TPM.S2.E10.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E1_M"), ])$Name)]
TPM.S2.E10.M.mean <- data.frame(S2_hTE_E1_M=rowMeans(TPM.S2.E10.M))

TPM.S2.E10.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E10_P"), ])$Name)]
TPM.S2.E10.P.mean <- data.frame(S2_hTE_E10_P=rowMeans(TPM.S2.E10.P))

TPM.S2.E11.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E11_M"), ])$Name)]
TPM.S2.E11.M.mean <- data.frame(S2_hTE_E11_M=rowMeans(TPM.S2.E11.M))

TPM.S2.E11.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E11_P"), ])$Name)]
TPM.S2.E11.P.mean <- data.frame(S2_hTE_E11_P=rowMeans(TPM.S2.E11.P))

TPM.S2.E12.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E12_M"), ])$Name)]
TPM.S2.E12.M.mean <- data.frame(S2_hTE_E12_M=rowMeans(TPM.S2.E12.M))

TPM.S2.E12.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E12_P"), ])$Name)]
TPM.S2.E12.P.mean <- data.frame(S2_hTE_E12_P=rowMeans(TPM.S2.E12.P))

TPM.S2.E14.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E13_M"), ])$Name)]
TPM.S2.E14.P.mean <- data.frame(S2_hTE_E13_M=rowMeans(TPM.S2.E14.P))

TPM.S2.E15.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E13_P"), ])$Name)]
TPM.S2.E15.M.mean <- data.frame(S2_hTE_E13_P=rowMeans(TPM.S2.E15.M))

TPM.S2.E15.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E2_M"), ])$Name)]
TPM.S2.E15.P.mean <- data.frame(S2_hTE_E2_M=rowMeans(TPM.S2.E15.P))

TPM.S2.E16.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E2_P"), ])$Name)]
TPM.S2.E16.M.mean <- data.frame(S2_hTE_E2_P=TPM.S2.E16.M)

TPM.S2.E16.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E24_P"), ])$Name)]
TPM.S2.E16.P.mean <- data.frame(S2_hTE_E24_P=rowMeans(TPM.S2.E16.P))


TPM.S2.E28.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E3_M"), ])$Name)]
TPM.S2.E28.P.mean <- data.frame(S2_hTE_E3_M=rowMeans(TPM.S2.E28.P))

TPM.S2.E29.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E3_P"), ])$Name)]
TPM.S2.E29.M.mean <- data.frame(S2_hTE_E3_P=rowMeans(TPM.S2.E29.M))

TPM.S2.E29.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E4_M"), ])$Name)]
TPM.S2.E29.P.mean <- data.frame(S2_hTE_E4_M=rowMeans(TPM.S2.E29.P))

TPM.S2.E6.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E4_P"), ])$Name)]
TPM.S2.E6.M.mean <- data.frame(S2_hTE_E4_P=rowMeans(TPM.S2.E6.M))

TPM.S2.E6.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E5_M"), ])$Name)]
TPM.S2.E6.P.mean <- data.frame(S2_hTE_E5_M=rowMeans(TPM.S2.E6.P))

TPM.S2.E7.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E5_P"), ])$Name)]
TPM.S2.E7.M.mean <- data.frame(S2_hTE_E5_P=rowMeans(TPM.S2.E7.M))

TPM.S2.E7.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E6_M"), ])$Name)]
TPM.S2.E7.P.mean <- data.frame(S1_hTE_E6_M=rowMeans(TPM.S2.E7.P))

TPM.S2.E8.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E6_P"), ])$Name)]
TPM.S2.E8.M.mean <- data.frame(S1_hTE_E6_P=rowMeans(TPM.S2.E8.M))

TPM.S2.E8.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E7_M"), ])$Name)]
TPM.S2.E8.P.mean <- data.frame(S2_hTE_E7_M=rowMeans(TPM.S2.E8.P))

TPM.S2.E9.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E7_P"), ])$Name)]
TPM.S2.E9.M.mean <- data.frame(S2_hTE_E7_P=rowMeans(TPM.S2.E9.M))

TPM.S2.E9.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E8_M"), ])$Name)]
TPM.S2.E9.P.mean <- data.frame(S2_hTE_E8_M=rowMeans(TPM.S2.E9.P))

TPM.S3.E19.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S2_E8_P"), ])$Name)]
TPM.S3.E19.M.mean <- data.frame(S2_hTE_E8_P=rowMeans(TPM.S3.E19.M))

TPM.S3.E19.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E15_M"), ])$Name)]
TPM.S3.E19.P.mean <- data.frame(S3_hTE_E15_M=rowMeans(TPM.S3.E19.P))

TPM.S3.E21.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E15_P"), ])$Name)]
TPM.S3.E21.M.mean <- data.frame(S3_hTE_E15_P=rowMeans(TPM.S3.E21.M))

TPM.S3.E21.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E16_ISH"), ])$Name)]
TPM.S3.E21.P.mean <- data.frame(S3_hTE_E16_ISH=rowMeans(TPM.S3.E21.P))

TPM.S3.E22.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E16_P"), ])$Name)]
TPM.S3.E22.M.mean <- data.frame(S3_hTE_E16_P=rowMeans(TPM.S3.E22.M))

TPM.S3.E22.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E17_M"), ])$Name)]
TPM.S3.E22.P.mean <- data.frame(S3_hTE_E17_M=rowMeans(TPM.S3.E22.P))

TPM.S3.E20.ISH <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E17_P"), ])$Name)]
TPM.S3.E20.ISH.mean <- data.frame(S3_hTE_E17_P=rowMeans(TPM.S3.E20.ISH))

TPM.S3.E20.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E18_M"), ])$Name)]
TPM.S3.E20.P.mean <- data.frame(S3_hTE_E18_M=rowMeans(TPM.S3.E20.P))


TPM.S3.E27.ISH <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E18_P"), ])$Name)]
TPM.S3.E27.ISH.mean <- data.frame(S3_hTE_E18_P=rowMeans(TPM.S3.E27.ISH))

TPM.S3.E27.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E23_ISH"), ])$Name)]
TPM.S3.E27.P.mean <- data.frame(S3_hTE_E23_ISH=rowMeans(TPM.S3.E27.P))

TPM.S3.E27.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E23_M"), ])$Name)]
TPM.S3.E27.M.mean <- data.frame(S3_hTE_E23_M=rowMeans(TPM.S3.E27.M))

TPM.S1.E31.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E23_P"), ])$Name)]
TPM.S1.E31.M.mean <- data.frame(S3_hTE_E23_P=rowMeans(TPM.S1.E31.M))

#TPM.ISK <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "ISK_ISK_ISK"), ])$name)]
#TPM.ISK.mean <- data.frame(ISK=rowMeans(TPM.ISK))


#TPM.S1.E31.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S1_E31_P"), ])$name)]
#TPM.S1.E31.P.mean <- data.frame(S1_hTE_E31_P=rowMeans(TPM.S1.E31.P))

#TPM.S3.E32.M <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E32_M"), ])$name)]
#TPM.S3.E32.M.mean <- data.frame(S3_hTE_E32_M=rowMeans(TPM.S3.E32.M))

#TPM.S3.E32.P <- TPM[ ,which(colnames(TPM) %in% (summary[which(summary$class %in% "S3_E32_P"), ])$name)]
#TPM.S3.E32.P.mean <- data.frame(S3_hTE_E32_P=rowMeans(TPM.S3.E32.P))

filename <- data.frame(name=table(summary$class))
filename[ ,1] <- gsub("_",".",filename[ ,1]) 
filename$rename <- paste("TPM",filename[ ,1],"mean",sep = ".")


filedata <- cbind(TPM.S1.E1.M.mean,TPM.S1.E13.M.mean,TPM.S1.E13.P.mean,TPM.S1.E18.M.mean,TPM.S1.E18.P.mean,TPM.S1.E23.M.mean,TPM.S1.E23.P.mean,TPM.S1.E24.M.mean,TPM.S1.E24.P.mean,TPM.S1.E26.M.mean,TPM.S1.E26.P.mean,TPM.S1.E30.M.mean,TPM.S1.E30.P.mean,TPM.S2.E10.M.mean,TPM.S2.E10.P.mean,TPM.S2.E11.M.mean,TPM.S2.E11.P.mean,TPM.S2.E12.M.mean,TPM.S2.E12.P.mean,TPM.S2.E14.P.mean,TPM.S2.E15.M.mean,TPM.S2.E15.P.mean,TPM.S2.E16.M.mean,TPM.S2.E16.P.mean,TPM.S2.E28.P.mean,TPM.S2.E29.M.mean,TPM.S2.E29.P.mean,TPM.S2.E6.M.mean,TPM.S2.E6.P.mean,TPM.S2.E7.M.mean,TPM.S2.E7.P.mean,TPM.S2.E8.M.mean,TPM.S2.E8.P.mean,TPM.S2.E9.M.mean,TPM.S2.E9.P.mean,TPM.S3.E19.M.mean,TPM.S3.E19.P.mean,TPM.S3.E20.ISH.mean,TPM.S3.E20.P.mean,TPM.S3.E21.M.mean,TPM.S3.E21.P.mean,TPM.S3.E22.M.mean,TPM.S3.E22.P.mean,TPM.S3.E27.ISH.mean,TPM.S3.E27.M.mean,TPM.S3.E27.P.mean,TPM.S1.E31.M.mean)

filedata.filt <- filedata[which(rowMeans(filedata)>1), ]


hc.candidates <- hclust(as.dist(1-abs(cor(filedata.filt,method="spearman"))), method="complete")# ward.D2 

pdf("../Figure/Figure2a.pdf",width = 12,height = 4.5)
plot(hc.candidates,hang=-1,col="black",cex=1)
dev.off()

#############STAGE1-3 marker
library(Seurat)
library(dplyr)
TPM.class <- sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 1 )
#TPM.no.ISK <- TPM[ ,-which(TPM.class %in% c("ISK"))]
nbt.stage <- CreateSeuratObject(TPM, min.cells = 3, min.genes = 2000, project = "NBT")
nbt.stage <- NormalizeData(object = nbt.stage, normalization.method = "LogNormalize",scale.factor = 10000)
nbt.stage <- FindVariableGenes(object = nbt.stage, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.6, y.cutoff = 0.5)#, x.high.cutoff = 3
length(x = nbt.stage@var.genes)
#nbt.stage.markers1 <- FindAllMarkers(object = nbt.stage, only.pos = TRUE, min.pct = 0.25,  
                                    #thresh.use=1.5)
nbt.stage.markers <- FindAllMarkers(object = nbt.stage, only.pos = TRUE, min.pct = 0.25,  thresh.use=2)
#nbt.stage.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
nbt.stage.all.markers=FindAllMarkers(nbt.stage,test.use = "roc")
nbt.stage.markers.use=subset(nbt.stage.all.markers,avg_diff>0&power>0.4)

write.table(nbt.stage.markers,"nbt.stage.markers.use.xls",quote = F,col.names = T,row.names = F,sep = "\t")
write.table(nbt.stage.markers.use,"nbt.stage.markers.use.old.xls",quote = F,col.names = T,row.names = F,sep = "\t")
nbt.stage.markers.use.gene <- as.vector(nbt.stage.markers.use$gene)
TPM.no.ISK.class <-sapply( strsplit(as.character(colnames(TPM)), "_"), "[[", 1 )
TPM.no.ISK.S1 <- TPM[ ,which(TPM.no.ISK.class %in% c("S1"))]
TPM.no.ISK.S2 <- TPM[ ,which(TPM.no.ISK.class %in% c("S2"))]
TPM.no.ISK.S3 <- TPM[ ,which(TPM.no.ISK.class %in% c("S3"))]
TPM.Stage <- cbind(TPM.no.ISK.S1,TPM.no.ISK.S2,TPM.no.ISK.S3)
nbt.stage.markers.use.tpm <- TPM.Stage[nbt.stage.markers.use.gene,]
mean <- rowMeans(nbt.stage.markers.use.tpm)
sd_value <- apply(nbt.stage.markers.use.tpm,1,sd)

z_score <- (nbt.stage.markers.use.tpm-mean)/sd_value
z_score[z_score>=1.5] <- 1.5
z_score[z_score<= -1.5] <- -1.5

ubj1 <- sapply( strsplit(as.character(colnames(z_score)), "/t"), "[[", 1 )
aaka22 = data.frame(sapply( strsplit(as.character(colnames(z_score)), "_"), "[[", 1 ))

row.names(aaka22)<-ubj1
#colnames(aaka3)<-c("cluster")
colnames(aaka22)<-c("Class")

aka4=list(Class=c("S1"="#37ACCB","S2"="#FB9100","S3"="#B0C915"))

library('RColorBrewer')
library("pheatmap")
pheatmap( z_score,annotation_col =aaka22,annotation_colors = aka4,
          color = colorRampPalette(rev(brewer.pal(n = 11, name ="RdYlBu")))(50), 
          #breaks = c(seq(-2,-1.1,length=15),seq(-1,1,length=70),seq(1.1,2,length=15)),
          cluster_row = F,cluster_col =F,show_rownames = F, 
          show_colnames =F, scale = "none",legend =T,border_color =T,
          cellwidth = 0.2,cellheight = 0.5,fontsize_row=8,
          fontsize=8,fontsize_col = 5,file="../new_plot/stage_marker_heatmap2.pdf"
)




#############STAGE1-3 specific gene function
library(cummeRbund)
TPM.S1.mean <- data.frame(Stage1=rowMeans(TPM.no.ISK.S1))
TPM.S2.mean <- data.frame(Stage2=rowMeans(TPM.no.ISK.S2))
TPM.S3.mean <- data.frame(Stage3=rowMeans(TPM.no.ISK.S3))
TPM.Stage.mean <- cbind(TPM.S1.mean,TPM.S2.mean,TPM.S3.mean)
in_data=TPM.Stage.mean[rowMeans(TPM.Stage.mean)>1,]
genesCluster<-function (object, k, logMode = T, pseudocount = 1, ...) 
{
  .local <- function (object, k, logMode = T, method = "none", 
                      pseudocount = 1, ...) 
  {
    require(cluster)
    m <- as.data.frame(object)
    m <- m[rowSums(m) > 0, ]
    if (logMode) {
      m <- log10(m + pseudocount)
    }
    if (!is.function(method)) {
      method = function(mat) {
        JSdist(makeprobs(t(m)))
      }
    }
    n <- method(m)
    clusters <- pam(n, k, ...)
    class(clusters) <- "list"
    clusters$fpkm <- m
    clusters
  }
  .local(object, k, logMode, pseudocount, ...)
}


cluster_result1<-genesCluster(in_data,k=9)
gene_cluster1<-csClusterPlot(cluster_result1)
pdf("../Figure/Figure2e.pdf",width=5,height=5)
gene_cluster1+scale_color_manual(values = c("#66CCCC","#99CC66","#FFFF66","#666699","#CCFF66","#FFCC99","#FF9900","#FF99CC","#FF9999"))
dev.off()

gene_groups<-as.data.frame(cluster_result1$clustering)
gene_groups$genes<-rownames(gene_groups)
colnames(gene_groups)=c("cluster","genes")
gene_cluster.list<-gene_groups[order(gene_groups$cluster),]
write.table(gene_cluster.list,"../Table/Stage_gene_list.txt",row.names=F,col.names=TRUE,sep="\t",quote=F)


