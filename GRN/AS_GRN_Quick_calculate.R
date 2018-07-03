### Quick calculation
#source("~/GENIE3_R/GENIE3_python/Gene_AS_GRN_comparison.R")

### GOseq enrichment test function
Go_goseq <- function(gene_list){
  library(goseq)
  Godb <- read.delim("/mnt/01/hungying/Sequence/go_ensembl_zea_mays.gaf", header = T, sep = "\t")
  rawreads <- read.table("/mnt/01/hungying/RD_GWAS/ref/NSF-SAM.380.individuals.FGSv5b.60.readcounts.txt", header=T, sep="\t") #loading Raw count
  rawreads$go.term <- 0
  rawreads$go.term[rawreads$Gene %in% gene_list] <- 1
  Goterm <- rawreads[,c("Gene","go.term")]
  sum(rawreads$go.term)
  gene.vector <- Goterm$go.term
  names(gene.vector) <- Goterm$Gene
  rawreads$sum <- rowSums(rawreads[,6:386])
  gene.len <- read.delim("/mnt/01/hungying/Sequence/transcript_length", sep="\t", header=F)
  colnames(gene.len)<-c("ID","length")
  rownames(gene.len) <- as.character(gene.len$ID)
  test.geneLen <- as.integer(gene.len$length)
  names(test.geneLen) <- rownames(gene.len)
  pwf.counts <- nullp(gene.vector, bias.data=test.geneLen, plot.fit=T)
  go <- goseq(pwf.counts, id=rownames(rawreads), gene2cat=Godb[,c("db_object_symbol","term_accession")],method="Sampling",repcnt = 100000)
  go <- go[go$over_represented_pvalue < 0.05,]
  return(go)
}

Mapman_Func_Enrich <- function(gene_list){
  i = 1
  source("/mnt/01/hungying/cassy.folder/mapman.enrichment.v2.3.R")
  mapman.annotation.file = "/mnt/01/hungying/cassy.folder/Zm_B73_5b_FGS_cds_2011.overrepresented.test.txt"
  rawreads <- read.table("/mnt/01/hungying/RD_GWAS/ref/NSF-SAM.380.individuals.FGSv5b.60.readcounts.txt", header=T, sep="\t") #loading Raw count
  rawreads$go.term <- 0
  rawreads$go.term[rawreads$Gene %in% gene_list] <- 1 
  Goterm <- rawreads[,c("Gene","go.term")]
  colnames(Goterm) <- c("genes","ts.cate")
  write.table(Goterm, "~/Sequence/temp", sep="\t", col.names = T, row.names = F, quote = F)
  result <- ORT2("~/Sequence/temp",mapman.annotation.file,i)
  return(result)
}

Po_goseq <- function(gene_list){
  library(goseq)
  Godb <- read.delim("/mnt/01/hungying/Sequence/PO_gene2cat.tab", header = T, sep = "\t")
  colnames(Godb) <- c("ID","PO")
  rawreads <- read.table("/mnt/01/hungying/RD_GWAS/ref/NSF-SAM.380.individuals.FGSv5b.60.readcounts.txt", header=T, sep="\t") #loading Raw count
  rawreads$go.term <- 0
  rawreads$go.term[rawreads$Gene %in% gene_list] <- 1
  Goterm <- rawreads[,c("Gene","go.term")]
  sum(rawreads$go.term)
  gene.vector <- Goterm$go.term
  names(gene.vector) <- Goterm$Gene
  rawreads$sum <- rowSums(rawreads[,6:386])
  gene.len <- read.delim("/mnt/01/hungying/Sequence/transcript_length", sep="\t", header=F)
  colnames(gene.len)<-c("ID","length")
  rownames(gene.len) <- as.character(gene.len$ID)
  test.geneLen <- as.integer(gene.len$length)
  names(test.geneLen) <- rownames(gene.len)
  pwf.counts <- nullp(gene.vector, bias.data=test.geneLen, plot.fit=T)
  PO <- goseq(pwf.counts, id=rownames(rawreads), gene2cat=Godb,repcnt = 100000,method="Sampling")
  PO <- PO[PO$over_represented_pvalue < 0.05,]
  return(PO)
}
des <- read.delim("/mnt/01/hungying/RD_GWAS/ref/Zmays_284_6a.annotation_info.txt", header=F)
names(des) <- c("PhytoID", "geneid", "txid", "proid", "PFAM", "Panther", "KOG", "KEGGec",
                "KEGGorth", "GO", "ara_hit_id", "ara_hit_symbol", "ara_hit_defline",
                "rice_hit_id", "rice_hit_symbol", "rice_hit_defile")
maizeGDB <- read.csv("/mnt/01/hungying/Sequence/Zea_mays.AGPv3.gene.gff3", sep="\t", header=F)
maizeGDB$V10 <- gsub("^ID.+;biotype", "",maizeGDB$V10)
maizeGDB$V10 <- gsub(";logic_name.+$", "",maizeGDB$V10)
Get_anno <- function(list_before,target){
  gene_Reg_anno <- merge(list_before, des, by.x=target, by.y="geneid", all=F)
  gene_Reg_anno <- merge(gene_Reg_anno[,c(1,2,3,4,5,13,14,15,18)], maizeGDB[,c(1,10)], by.x=target, by.y="V1")
  #gene_Reg_anno <- gene_Reg_anno[!duplicated(gene_Reg_anno$Tar),]
  gene_Reg_anno <- gene_Reg_anno[order(gene_Reg_anno$Score, decreasing=TRUE),]
  return(gene_Reg_anno)
}

### Analysis network output  (GRMZM2G079823;intron_retention;26201313;26202129) (GRMZM2G079823;intron_retention;26201313;26202825)
ASgene_net <- read.table("/mnt/01/hungying/GENIE3_R/GENIE3_python/SAM_AS_net_noASAS_Top05.txt", sep = "\t", header = F)
colnames(ASgene_net) <- c("Reg", "Tar", "Score")
AS_signals <- read.delim("/mnt/01/hungying/AS_GWAS/AS_GRN/Full_AS_list.txt")
for(i in 1:58){
  #i = 1
  AS_name <- as.character(AS_signals$x[i])
  AS_gene <- as.character(sub(";.+$","",AS_name))
  #Upstream analysis
  gene_Tar1 <- ASgene_net[ASgene_net$Tar == AS_gene, ]
  AS_Tar1 <- ASgene_net[ASgene_net$Tar == AS_name, ]
  gene_Tar1_anno <- Get_anno(gene_Tar1, "Reg")
  AS_Tar1_anno <- Get_anno(AS_Tar1, "Reg")
  write.table(gene_Tar1_anno, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_gene_anno/", AS_gene, ".txt", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(AS_Tar1_anno, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_gene_anno/", AS_name, ".txt", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  #Overlap
  up <- matrix(0, ncol = 5, nrow = 1)
  up <- data.frame(up)
  names(up) <- c("Gene", "AS", "Common", "Unique_gene_Gene%", "Unique_gene_AS%")
  gene_Tar1_anno <- gene_Tar1_anno[!duplicated(gene_Tar1_anno$Reg),]
  AS_Tar1_anno <- AS_Tar1_anno[!duplicated(AS_Tar1_anno$Reg),]
  up$Gene <- length(gene_Tar1_anno$Reg)
  up$AS <- length(AS_Tar1_anno$Reg)
  up$Common <- sum(AS_Tar1_anno$Reg %in% gene_Tar1_anno$Reg)
  up$`Unique_gene_Gene%` <- (1-(up$Common / up$Gene)) * 100
  up$`Unique_gene_AS%` <- (1-(up$Common / up$AS)) * 100
  write.table(up, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_overlap/", AS_name, ".txt", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  #GO ,PO and Mapman
  AS_Tar1_go <- Go_goseq(AS_Tar1_anno$Reg)
  AS_Tar1_Map <- Mapman_Func_Enrich(AS_Tar1_anno$Reg)
  AS_Tar1_po <- Po_goseq(AS_Tar1_anno$Reg)
  gene_Tar1_go <- Go_goseq(gene_Tar1_anno$Reg)
  gene_Tar1_Map <- Mapman_Func_Enrich(gene_Tar1_anno$Reg)
  gene_Tar1_po <- Po_goseq(gene_Tar1_anno$Reg)
  write.table(AS_Tar1_go, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_go/", AS_name, ".go", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(AS_Tar1_Map, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_map/", AS_name, ".map", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(AS_Tar1_po, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_po/", AS_name, ".po", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(gene_Tar1_go, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_go/", AS_gene, ".go", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(gene_Tar1_Map, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_map/", AS_gene, ".map", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(gene_Tar1_po, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/up_po/", AS_gene, ".po", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  #Downstream
  gene_Reg1 <- ASgene_net[ASgene_net$Reg == AS_gene, ]
  AS_Reg1 <- ASgene_net[ASgene_net$Reg == AS_name, ]
  gene_Reg1_anno <- Get_anno(gene_Reg1, "Tar")
  AS_Reg1_anno <- Get_anno(AS_Reg1, "Tar")
  write.table(gene_Reg1_anno, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_gene_anno/", AS_gene, ".txt", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(AS_Reg1_anno, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_gene_anno/", AS_name, ".txt", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  #Overlap
  down <- matrix(0, ncol = 5, nrow = 1)
  down <- data.frame(down)
  names(down) <- c("Gene", "AS", "Common", "Unique_gene_Gene%", "Unique_gene_AS%")
  gene_Reg1_anno <- gene_Reg1_anno[!duplicated(gene_Reg1_anno$Tar),]
  AS_Reg1_anno <- AS_Reg1_anno[!duplicated(AS_Reg1_anno$Tar),]
  down$Gene <- length(gene_Reg1_anno$Tar)
  down$AS <- length(AS_Reg1_anno$Tar)
  down$Common <- sum(AS_Reg1_anno$Tar %in% gene_Reg1_anno$Tar)
  down$`Unique_gene_Gene%` <- (1-(down$Common / down$Gene)) * 100
  down$`Unique_gene_AS%` <- (1-(down$Common / down$AS)) * 100
  write.table(down, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_overlap/", AS_name, ".txt", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  #Go, PO, Mapman
  AS_Reg1_go <- Go_goseq(AS_Reg1_anno$Tar)
  AS_Reg1_Map <- Mapman_Func_Enrich(AS_Reg1_anno$Tar)
  AS_Reg1_po <- Po_goseq(AS_Reg1_anno$Tar)
  gene_Reg1_go <- Go_goseq(gene_Reg1_anno$Tar)
  gene_Reg1_Map <- Mapman_Func_Enrich(gene_Reg1_anno$Tar)
  gene_Reg1_po <- Po_goseq(gene_Reg1_anno$Tar)
  write.table(AS_Reg1_go, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_go/", AS_name, ".go", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(AS_Reg1_Map, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_map/", AS_name, ".map", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(AS_Reg1_po, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_po/", AS_name, ".po", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(gene_Reg1_go, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_go/", AS_gene, ".go", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(gene_Reg1_Map, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_map/", AS_gene, ".map", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
  write.table(gene_Reg1_po, paste("/mnt/01/hungying/AS_GWAS/AS_GRN/down_po/", AS_gene, ".po", sep=""), 
              quote = F, sep = "\t", col.names = T,row.names = F)
}

