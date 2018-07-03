#Install packages (Do this section only for new installation of R)
#-------------------------------------------------------------------------------
source("http://www.bioconductor.org/biocLite.R") 
biocLite("multtest")
install.packages("gplots")
#Step 0: Import library and GAPIT functions run this section each time to start R)
#######################################################################################
library('MASS') # required for ginv
library(multtest)
library(gplots)
library(scatterplot3d)
library(compiler) #required for cmpfun
source("http://www.zzlab.net/GAPIT/emma.txt")
source("/mnt/01/hungying/cassy.folder/gapit_functions.txt")

## Testing AS_GWAS SNP PCA control ##
geno <- read.delim("/mnt/01/hungying/AS_GWAS/Share/MLM-AS-GWAS/Event.list.071615.GAPit0075.hmp", head = FALSE)
phenoG <- read.delim("//mnt/01/hungying/AS_GWAS/Share/MLM-AS-GWAS/GAPIT_pheno.txt", header = T, sep="\t") 
myCV <- read.table("/mnt/01/hungying/cassy.folder/SAM_GAPIT_CMLM1/GAPIT.PCA.csv",header =T, sep=",") #Using SNP as PCA source
setwd("/mnt/01/hungying/AS_GWAS/Share/MLM-AS-GWAS/")
myGAPIT_SUPER <- GAPIT(
  Y=phenoG[,c(1,23)],  		
  G=geno,				
  #CV=myCV,
  #PCA.total=1,				
  sangwich.top="MLM", #options are GLM,MLM,CMLM, FaST and SUPER 
  sangwich.bottom="SUPER", #options are GLM,MLM,CMLM, FaST and SUPER 
  group.from = 100,
  group.to = 400,
  LD=0.1
)

