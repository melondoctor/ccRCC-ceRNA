#install.packages("pheatmap")


library(pheatmap)       #引用包
setwd("C:\\Users\\melon doctor\\Desktop\\biotech\\log\\19.ceRNAheatmap")     #设置工作目录

#定义热图函数
bioHeatmap=function(expFile=null,geneFile=null,outFile=null,height=null,width=null){
	#读取输入文件，并对输入文件整理
	rt=read.table(expFile,sep="\t",header=T,check.names=F,row.names=1)
	gene=read.table(geneFile,sep="\t",header=F,check.names=F)
	data=rt[as.vector(gene[,1]),]
	Type=gsub("(.*?)\\_(.*)","\\2",colnames(data))
	Type=ifelse(Type=="con","C","T")
	colnames(data)=gsub("(.*?)\\_(.*)","\\1",colnames(data))
	
	#绘制热图
	names(Type)=colnames(data)
	Type=as.data.frame(Type)
	pdf(file=outFile,height=height,width=width)
	p=pheatmap(data, 
	         annotation=Type, 
	         color = colorRampPalette(c("blue", "white", "red"))(50),
	         cluster_cols =F,
	         show_colnames = F,
	         scale="row",
	         fontsize = 8,
	         fontsize_row=7,
	         fontsize_col=8)
	print(p)
	dev.off()
}

#绘制circRNA热图
bioHeatmap(expFile="circ.normalize.txt",
           geneFile="ceRNA.circList.txt",
           outFile="circ.heatmap.pdf",
           height=4, width=6)

#绘制miRNA热图
bioHeatmap(expFile="mirna.normalize.txt",
           geneFile="ceRNA.mirnaList.txt",
           outFile="mirna.heatmap.pdf",
           height=4, width=6)

#绘制mRNA热图
bioHeatmap(expFile="mrna.normalize.txt",
           geneFile="ceRNA.mrnaList.txt",
           outFile="mrna.heatmap.pdf",
           height=8, width=6)