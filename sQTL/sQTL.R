## eQTL_for AS_GWAS event
## Analysis Setting
library(plyr)
library(MatrixEQTL)
## Location of the package with the data files.
base.dir = "/mnt/01/hungying/cassy.folder/AS_GWAS/"
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

SNP_file_name = paste(base.dir, "SNP_matrixQTL_SAMalt_m2.txt", sep="");
snps_location_file_name = paste(base.dir, "SNPpos_matrixQTL_SAM.txt", sep="");
# Gene expression file name
expression_file_name = paste(base.dir, "AS_count_snp_m2_maf05.txt", sep="");
gene_location_file_name = paste(base.dir, "SAMAltloc.txt", sep="");
# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste(base.dir, "PCA_alt_m2.txt", sep=""); # 20 PCs

# Output file name
output_file_name_cis = tempfile(pattern="Cis.out");
output_file_name_tra = tempfile(pattern="Tra.out");
# Only associations significant at this level will be saved
pvOutputThreshold_cis = 0.05/(860999);
pvOutputThreshold_tra = 0.05/(860999);
# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
#errorCovariance = read.table("Sample_Data/errorCovariance.txt");
# Distance for local gene-SNP pairs
cisDist = 1e6;
# Output file name
output_file_name_cis = tempfile(pattern="Cis.out");
output_file_name_tra = tempfile(pattern="Tra.out");
## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t"; # the TAB character
snps$fileOmitCharacters = "99"; # denote missing values;
snps$fileSkipRows = 1; # one row of column labels
snps$fileSkipColumns = 1; # one column of row labels
snps$fileSliceSize = 100000; # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);
## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t"; # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1; # one row of column labels
gene$fileSkipColumns = 1; # one column of row labels
gene$fileSliceSize = 5000; # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t"; # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1; # one row of column labels
cvrt$fileSkipColumns = 1; # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}
## Minor allele frequency filtering
maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

# Look at the distribution of MAF
hist(maf[maf<0.1],seq(0,0.1,length.out=100))
cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.1);
snps$RowReorder(maf>0.05);
cat('SNPs before filtering:',nrow(snps))

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

## Create an artificial dataset and plot the histogram and Q-Q plot of all p-values
# a Q-Q plot
meq = Matrix_eQTL_engine(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt, 
  output_file_name = output_file_name_cis,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = modelLINEAR, 
  errorCovariance = numeric(), 
  verbose = TRUE,
  pvalue.hist = "qqplot")

plot(meq, pch = 16, cex = 0.7)

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);


## Results:  #par(resetPar()) 
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected local eQTLs:', '\n');
show(me$cis$eqtls)
cat('Detected distant eQTLs:', '\n');
show(me$trans$eqtls)
## Make the histogram of local and distant p-values
plot(meq)
#load("/mnt/01/hungying/AS_GWAS/as_eQTL_analysis.RData")
head(meq$all$eqtls$snps)
class(meq$all$eqtls$snps)
length(meq$all$eqtls$snps)
length(meq$all$eqtls$pvalue)