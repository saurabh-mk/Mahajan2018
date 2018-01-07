## author: Saurabh Mahajan
## created: Dec 30, 2016
## description: We have two sources of data for rRNA CN- IMG database and rrnDB.
## There is some discrepancy between the two. Since rrnDB is curated and more reliable,
## I will replace rRNA CN data from IMG with rRNA CN from rrnDB, wherever possible.
## depends on: IMG DB data (downloaded May 5, 2016); rrnDB 5.1

## get data from IMG DB
IMG_data <- read.table(file = "traits_data/all_genomes_IMG_data.txt", header = T, stringsAsFactors = F)
rrnaCN_IMG <- round(IMG_data$rRNA.Count.....assembled/3, digits = 2)
names(rrnaCN_IMG) <- rownames(IMG_data)

##read the available rRNA DB data
rrnDB_data <- read.delim(file = "traits_data/rrnDB-5.1.tsv", stringsAsFactors = F)

## create lists of just the refseq/genbank genome ids for taxa of interest, and for rRNA DB data
## this will serve as a common identifier between the two datasets
GCF_num_from_IMG <- sapply(rownames(IMG_data), function(x) strsplit(x = x, split = "_")[[1]][2])
GCF_num_from_IMG2 <- sapply(GCF_num_from_IMG, function(x) strsplit(x = x, split = "[.]")[[1]][1])
GCF_from_rrnDB <- rrnDB_data$Data.source.record.id
GCF_num_from_rrnDB <- sapply(GCF_from_rrnDB, function(x) strsplit(x = x, split = "_")[[1]][2])
GCF_num_from_rrnDB2 <- sapply(GCF_num_from_rrnDB, function(x) strsplit(x = x, split = "[.]")[[1]][1])

## find GCF numbers of taxa of interest for which data is avaialble in rRNA DB
common_GCF_num <- intersect(GCF_num_from_IMG2, GCF_num_from_rrnDB2)

## get rows in rrnaDB dataset corresponding to taxa that we require
row_oi <- sapply(GCF_num_from_IMG2, function(x) which(GCF_num_from_rrnDB2==x)[1])

## get rRNA CN data and write to file
rrnaCN_rrnDB <- rrnDB_data[row_oi, "X16S.gene.count"]
names(rrnaCN_rrnDB) <- rownames(IMG_data)

##wherever available, replace rRNA CN in IMG database with 16s numbers from rrnDB
rrnaCN_modified <- sapply(rownames(IMG_data), function(x) {ifelse(is.na(rrnaCN_rrnDB[x]), rrnaCN_IMG[x], rrnaCN_rrnDB[x])})
rrnaCN_final <- sapply(rrnaCN_modified, function(x) ifelse(x<1, 1, x))

rRNA_CN_data <- cbind.data.frame(rownames(IMG_data), IMG_data$Genome.Name.Sample.Name, rrnaCN_rrnDB, rrnaCN_IMG, rrnaCN_final)
write.table(x = rRNA_CN_data, file = "traits_data/rRNA_CN.txt", quote = T, row.names = F, col.names = c("genome_id","Genome_Name","rRNA_CN_rrnDB","rRNA_CN_IMG", "rRNA_CN_final"))