#install.packages("ggpubr")
#install.packages("reshape2")

#引用包
library(ggpubr)
library(reshape2)

setwd("C:\\Users\\melon doctor\\Desktop\\biotech\\log\\20.ceRNAboxplot")     #设置工作目录

#定义boxplot函数
bioBoxplot=function(expFile=null,geneFile=null,diffFile=null,outFile=null,height=null,width=null){
	#读取输入文件，并对输入文件整理
	rt=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)
	gene=read.table(geneFile,sep="\t",header=F,check.names=F)
	data=rt[as.vector(gene[,1]),]
	Type=gsub("(.*?)\\_(.*)","\\2",colnames(data))
	Type=ifelse(Type=="con","C","T")
	colnames(data)=gsub("(.*?)\\_(.*)","\\1",colnames(data))
	
	#读取差异文件
	diff=read.table(diffFile,sep="\t",header=T,check.names=F,row.names=1)
	diff=diff[as.vector(gene[,1]),]
	adjPval=as.character(diff$adj.P.Val)
	Sig=ifelse(adjPval<0.001,"***",ifelse(adjPval<0.01,"**",ifelse(adjPval<0.05,"*","")))
	
	data=t(data)
	data=cbind(as.data.frame(data),Type)

	data=melt(data,id.vars=c("Type"))
	colnames(data)=c("Type","Gene","Expression")
	yMax=max(data$Expression)    #定义显著性的y坐标值
	
	p=ggboxplot(data, x="Gene", y="Expression", color = "Type", orientation = "horizontal",
	     ylab="", add = "none", xlab="",width=0.8,
	     palette = c("blue","red"))+rotate_x_text(60)+
	     geom_text(data=diff,aes(label=Sig,x=row.names(diff), y=yMax), position=position_dodge(1), vjust=0)
	pdf(file=outFile,width=width,height=height)       #输出图片文件
	print(p)
    dev.off()
}


#绘制circRNA箱线图
bioBoxplot(expFile="circ.normalize.txt",
           geneFile="ceRNA.circList.txt",
           diffFile="circ.diff.txt",
           outFile="circ.boxplot.pdf",
           height=3, width=5)

#绘制miRNA箱线图
bioBoxplot(expFile="mirna.normalize.txt",
           geneFile="ceRNA.mirnaList.txt",
           diffFile="mirna.diff.txt",
           outFile="mirna.boxplot.pdf",
           height=4, width=5)

#绘制mRNA箱线图
bioBoxplot(expFile="mrna.normalize.txt",
           geneFile="ceRNA.mrnaList.txt",
           diffFile="mrna.diff.txt",
           outFile="mrna.boxplot.pdf",
           height=10, width=7)