## This scripts used to compare PCA result between SNP and AS
par(mfrow=c(1,2))
AS_PCA <- read.table("/mnt/01/hungying/cassy.folder/AS_GWAS/AS_GWAS_PCA.txt", head=TRUE)
#AS_PCA <- read.table("/mnt/01/hungying/cassy.folder/AS_GWAS/AS_all/GAPIT.PCA.csv", head=TRUE)
SNP_PCA <- read.table("/mnt/01/hungying/cassy.folder/SAM_GAPIT_CMLM1/GAPIT.PCA.csv",header =T, sep=",") #Using SNP as PCA source
line_info <- read.csv("/mnt/01/hungying/AS_GWAS/line_info.csv",header =T, sep="," )
AS_table <- merge(AS_PCA, line_info, by.x="X.Trait.", by.y="Inbred.line", all.x=T, all.y=F)
SNP_table <- merge(SNP_PCA[,c(1,2,3)], line_info, by.x="taxa", by.y="Inbred.line", all.x=T, all.y=F)

plot(SNP_table$PC1[SNP_table$Pop.structure=="unclassified"], SNP_table$PC2[SNP_table$Pop.structure=="unclassified"], col="gray", pch=20, cex=0.8,
     xlim=c(-250,600), ylim=c(-500, 500),xlab="PC1 (5.9%)",ylab="PC2 (3.3%)", main="PCA_SNP")
points(SNP_table$PC1[SNP_table$Pop.structure=="stiff stalk"], SNP_table$PC2[SNP_table$Pop.structure=="stiff stalk"], col="red", pch=20, cex=0.8)
points(SNP_table$PC1[SNP_table$Pop.structure=="non-stiff stalk"], SNP_table$PC2[SNP_table$Pop.structure=="non-stiff stalk"], col="blue", pch=20, cex=0.8)
points(SNP_table$PC1[SNP_table$Pop.structure=="tropical"], SNP_table$PC2[SNP_table$Pop.structure=="tropical"], col="purple", pch=20, cex=0.8)
points(SNP_table$PC1[SNP_table$Pop.structure=="sweet corn"], SNP_table$PC2[SNP_table$Pop.structure=="sweet corn"], col="yellow", pch=20, cex=0.8)
points(SNP_table$PC1[SNP_table$Pop.structure=="popcorn"], SNP_table$PC2[SNP_table$Pop.structure=="popcorn"], col="green", pch=20, cex=0.8)
points(SNP_table$PC1[SNP_table$Pop.structure=="landrace"], SNP_table$PC2[SNP_table$Pop.structure=="landrace"], col="yellow", pch=20, cex=0.8)
legend(par("usr")[2]*0.95, par("usr")[3]*0.95, xjust=1, yjust=0, legend=c("unclassified", "stiff stalk", "non-stiff stalk  ", "tropical", "popcorn"),
       col=c("gray","red","blue", "purple", "green"), bg="white", pch=20, cex=0.6)

plot(AS_table$PC1[AS_table$Pop.structure=="unclassified"], AS_table$PC2[AS_table$Pop.structure=="unclassified"], col="gray", pch=20, cex=0.8,
     xlim=c(-150,110), ylim=c(-100, 100), xlab="PC1 (5.4%)",ylab="PC2 (2.4%)", main="PCA_AS")
points(AS_table$PC1[AS_table$Pop.structure=="stiff stalk"], AS_table$PC2[AS_table$Pop.structure=="stiff stalk"], col="red", pch=20, cex=0.8)
points(AS_table$PC1[AS_table$Pop.structure=="non-stiff stalk"], AS_table$PC2[AS_table$Pop.structure=="non-stiff stalk"], col="blue", pch=20, cex=0.8)
points(AS_table$PC1[AS_table$Pop.structure=="tropical"], AS_table$PC2[AS_table$Pop.structure=="tropical"], col="purple", pch=20, cex=0.8)
points(AS_table$PC1[AS_table$Pop.structure=="sweet corn"], AS_table$PC2[AS_table$Pop.structure=="sweet corn"], col="yellow", pch=20, cex=0.8)
points(AS_table$PC1[AS_table$Pop.structure=="popcorn"], AS_table$PC2[AS_table$Pop.structure=="popcorn"], col="green", pch=20, cex=0.8)
points(AS_table$PC1[AS_table$Pop.structure=="landrace"], AS_table$PC2[AS_table$Pop.structure=="landrace"], col="brown", pch=20, cex=0.8)
legend(par("usr")[2]*0.95, par("usr")[3]*0.95, xjust=1, yjust=0, legend=c("unclassified", "stiff stalk", "non-stiff stalk  ", "tropical", "popcorn"),
       col=c("gray","red","blue", "purple", "green"), bg="white", pch=20, cex=0.6)





