#!/usr/bin/ Rscript

#author: Saurabh Mahajan
#date: Oct 16, 2016
#description: this parses the output of ENCprime into tables which can be read by subsequent analysis scripts;
#outputs to same directory where the input is stored
#input arguments: directory path where genome specific ENCprime output is saved; genome to process

args_oi <- commandArgs(trailingOnly = TRUE)
dir_oi <- args_oi[1]
genome_oi <- args_oi[2]

enc_data <- readLines(con = paste0(dir_oi, "/", genome_oi, "_cds_from_genomic.fna.encp"))
enc_values <- unname(sapply(enc_data, function(dataline) strsplit(x = dataline, split = "]: ")[[1]][2]))
enc_names <- unname(sapply(enc_data, function(dataline) strsplit(x = dataline, split = "]: ")[[1]][1]))
enc_matrix <- do.call(rbind, strsplit(x = enc_values, split = " "))
rownames(enc_matrix) <- enc_names
colnames(enc_matrix) <- c("Nc", "Ncp","scChi","sumChi","df", "p","B_KM", "n_codons")
write.table(x = enc_matrix, file = paste0(dir_oi, "/", genome_oi, "_cds_from_genomic.fna.enc_matrix"), quote = T, sep = "\t", row.names = T)

codfreq_data <- readLines(con = paste0(dir_oi, "/", genome_oi, "_cds_from_genomic.fna.codfreq"))
codfreq_values <- unname(sapply(codfreq_data, function(dataline) strsplit(x = dataline, split = "]> ")[[1]][2]))
codfreq_names <- unname(sapply(codfreq_data, function(dataline) strsplit(x = dataline, split = "]> ")[[1]][1]))
codfreq_matrix <- do.call(rbind, strsplit(x = codfreq_values, split = " "))
rownames(codfreq_matrix) <- codfreq_names
codon_colnames <- c("TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG")
colnames(codfreq_matrix) <- strsplit(x = codon_colnames, split = " ")[[1]]
write.table(x = codfreq_matrix, file = paste0(dir_oi, "/", genome_oi,"_cds_from_genomic.fna.codfreq_matrix"), quote = T, sep = "\t", row.names = T)

codcnt_data <- readLines(con = paste0(dir_oi, "/", genome_oi, "_cds_from_genomic.fna.codcnt"))
codcnt_values <- unname(sapply(codcnt_data, function(dataline) strsplit(x = dataline, split = "]>")[[1]][2]))
codcnt_names <- unname(sapply(codcnt_data, function(dataline) strsplit(x = dataline, split = "]>")[[1]][1]))
codcnt_matrix <- do.call(rbind, strsplit(x = codcnt_values, split = " "))
rownames(codcnt_matrix) <- codcnt_names
codon_colnames <- c("TTT TTC TTA TTG TCT TCC TCA TCG TAT TAC TAA TAG TGT TGC TGA TGG CTT CTC CTA CTG CCT CCC CCA CCG CAT CAC CAA CAG CGT CGC CGA CGG ATT ATC ATA ATG ACT ACC ACA ACG AAT AAC AAA AAG AGT AGC AGA AGG GTT GTC GTA GTG GCT GCC GCA GCG GAT GAC GAA GAG GGT GGC GGA GGG")
colnames(codcnt_matrix) <- strsplit(x = codon_colnames, split = " ")[[1]]
write.table(x = codcnt_matrix, file = paste0(dir_oi, "/", genome_oi, "_cds_from_genomic.fna.codcnt_matrix"), quote = T, sep = "\t", row.names = T)

acgtfreq_data <- readLines(con = paste0(dir_oi, "/", genome_oi,"_cds_from_genomic.fna.acgtfreq"))
acgtfreq_values <- unname(sapply(acgtfreq_data, function(dataline) strsplit(x = dataline, split = "]> ")[[1]][2]))
acgtfreq_names <- unname(sapply(acgtfreq_data, function(dataline) strsplit(x = dataline, split = "]> ")[[1]][1]))
acgtfreq_matrix <- do.call(rbind, strsplit(x = acgtfreq_values, split = " "))
rownames(acgtfreq_matrix) <- acgtfreq_names
colnames(acgtfreq_matrix) <- c("fA", "fC", "fG", "fT")
write.table(x = acgtfreq_matrix, file = paste0(dir_oi, "/", genome_oi, "_cds_from_genomic.fna.acgtfreq_matrix"), quote = T, sep = "\t", row.names = T)