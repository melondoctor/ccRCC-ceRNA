#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

pvalueFilter=0.05         #p值过滤条件
qvalueFilter=1            #矫正后的p值过滤条件
setwd("C:\\Users\\melon doctor\\Desktop\\biotech\\log\\22.GO")       #设置工作目录

rt=read.table("id.txt",sep="\t",header=T,check.names=F)     #读取id.txt文件
rt=rt[is.na(rt[,"entrezID"])==F,]                           #去除基因id为NA的基因
gene=rt$entrezID

#定义颜色
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

for(i in c("BP","CC","MF")){
	#GO富集分析
	kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont=i, readable =T)
	GO=as.data.frame(kk)
    GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
	write.table(GO,file=paste0(i,".txt"),sep="\t",quote=F,row.names = F)          #保存富集结果
	
	#定义显示Term数目
	showNum=30
	if(nrow(GO)<30){
		showNum=nrow(GO)
	}

	#柱状图
	pdf(file=paste0(i,".barplot.pdf"),width = 11,height = 7)
	bar=barplot(kk, drop = TRUE, showCategory =showNum,color = colorSel)
	print(bar)
	dev.off()
		
	#气泡图
	pdf(file=paste0(i,".bubble.pdf"),width = 11,height = 7)
	bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio", color = colorSel)
	print(bub)
	dev.off()
}