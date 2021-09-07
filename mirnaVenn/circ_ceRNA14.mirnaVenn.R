#install.packages("venn")


library(venn)                         #引用包
outFile="intersect.miRNA.txt"         #输出文件名称
setwd("C:\\Users\\melon doctor\\Desktop\\biotech\\log\\14.mirnaVenn")    #设置工作目录
geneList=list()

#读取差异miRNA
rt=read.table("mirna.diff.txt",sep="\t",header=T,check.names=F)
geneNames=as.vector(rt[,1])                #提取miRNA名称
geneNames=gsub("^ | $","",geneNames)       #去掉miRNA首尾的空格
uniqGene=unique(geneNames)                 #miRNA取unique
geneList[["diff miRNA"]]=uniqGene

#读取MRE
rt=read.table("MRE.txt",sep="\t",header=T,check.names=F)
geneNames=as.vector(rt[,2])                #提取miRNA名称
geneNames=gsub("^ | $","",geneNames)       #去掉miRNA首尾的空格
uniqGene=unique(geneNames)                 #miRNA取unique
geneList[["MRE"]]=uniqGene

#绘制venn图
mycol=c("#029149","#E0367A","#5D90BA","#431A3D","#FFD121","#D8D155","#223D6C","#D20A13","#088247","#11AA4D","#7A142C","#5D90BA","#64495D","#7CC767")
pdf(file="miRNA.venn.pdf",width=5,height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

#保存交集miRNA
intersectGenes=Reduce(intersect,geneList)
write.table(file=outFile,intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)