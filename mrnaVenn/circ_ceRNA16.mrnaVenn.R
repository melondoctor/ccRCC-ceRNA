#install.packages("venn")


library(venn)                        #引用包
outFile="intersect.mRNA.txt"         #输出文件名称
setwd("C:\\Users\\melon doctor\\Desktop\\biotech\\log\\16.mrnaVenn")    #设置工作目录
geneList=list()

#读取差异mRNA
rt=read.table("mrna.diff.txt",sep="\t",header=T,check.names=F)
geneNames=as.vector(rt[,1])                #提取基因名称
geneNames=gsub("^ | $","",geneNames)       #去掉基因首尾的空格
uniqGene=unique(geneNames)                 #基因取unique
geneList[["diff mRNA"]]=uniqGene

#读取靶基因
rt=read.table("target.txt",sep="\t",header=T,check.names=F)
geneNames=as.vector(rt[,2])                #提取靶基因名称
geneNames=gsub("^ | $","",geneNames)       #去掉靶基因首尾的空格
uniqGene=unique(geneNames)                 #靶基因A取unique
geneList[["miRNA target"]]=uniqGene

#绘制venn图
mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="mRNA.venn.pdf",width=5,height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

#保存交集基因
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)