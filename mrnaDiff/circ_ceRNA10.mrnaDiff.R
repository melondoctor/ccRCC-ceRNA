#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")


#引用包
library(limma)
library(pheatmap)

inputFile="mrnaMatrix.txt"        #输入文件
conFile="sample1.txt"             #对照组样品
treatFile="sample2.txt"           #实验组样品
logFCfilter=1                     #logFC过滤阈值
adj.P.Val.Filter=0.05             #矫正后p值阈值
setwd("C:\\Users\\melon doctor\\Desktop\\biotech\\log\\10.mrnaDiff")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(inputFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=log2(data+1)       #如果表达数值很大，需要对数据取log2，可以这行前面#号删掉
data=normalizeBetweenArrays(data)

#读取样品信息
sample1=read.table(conFile,sep="\t",header=F,check.names=F)
sample2=read.table(treatFile,sep="\t",header=F,check.names=F)
conData=data[,as.vector(sample1[,1])]
treatData=data[,as.vector(sample2[,1])]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#差异分析
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut,file="mrna.all.xls",sep="\t",quote=F,col.names=F)

#输出矫正后的表达量
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file="mrna.normalize.txt",sep="\t",quote=F,col.names=F)

#输出差异结果
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file="mrna.diff.xls",sep="\t",quote=F,col.names=F)
write.table(diffSigOut,file="mrna.diff.txt",sep="\t",quote=F,col.names=F)

#绘制差异基因热图
geneNum=20
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep("C",conNum),rep("T",treatNum))
names(Type)=colnames(data)
Type=as.data.frame(Type)
pdf(file="mrna.heatmap.pdf",height=8,width=9)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 10,
         fontsize_row=7,
         fontsize_col=10)
dev.off()