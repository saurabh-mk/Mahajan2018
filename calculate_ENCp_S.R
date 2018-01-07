#title: "Calculate CUB metrics ENC' and S and amino acid usage"
#author: Saurabh Mahajan

##Choose parameters which decide which versions are calculated
ignore_ribo <- TRUE #if TRUE- defines HEG according to ENC' and ignores ribosomal+ genes
ENCp_percentile <- 0.01 #for HEG definition, in percentile

##Basic set up
tree_oi <- read.tree(file = "trees/bacteria_pruned_IMG_ids.newick")
genome_list <- read.csv(file = "genome_lists/all_curated_genomes_list.txt", header = F, stringsAsFactors = F)

codon_table <- read.csv(file = "codon_aa_map.txt", header = T, stringsAsFactors=F)
degen2_aas <- unique(codon_table[codon_table$degen==2, "aa"])
degen4_aas <- unique(codon_table[codon_table$degen==4, "aa"])
degen6_aas <- unique(codon_table[codon_table$degen==6, "aa"])
aa_unique <- setdiff(unique(codon_table$aa), c("Stop", "Met", "Trp"))
codons_aa <- sapply(aa_unique, function(x) codon_table[codon_table$aa==x, "codon"])
pos_list <- 1:3
base_list <- c("A", "C", "G", "T")

#dataframes to store genome specific and amino acid specific data
aawise_ENCp_HEG_all <- matrix(data = NA, nrow = length(genome_list$V1), ncol = length(aa_unique), dimnames = list(genome_list$V1, aa_unique))
aawise_ENCp_rest_all <- matrix(data = NA, nrow = length(genome_list$V1), ncol = length(aa_unique), dimnames = list(genome_list$V1, aa_unique))
ENCp_all <- matrix(data = NA, nrow = length(genome_list$V1), ncol = 3, dimnames = list(genome_list$V1, c("ENCp_HEG", "ENCp_rest", "deltaENCp")))
aawise_S_all <- matrix(data = NA, nrow = length(genome_list$V1), ncol = length(degen2_aas)+2, dimnames = list(genome_list$V1, c(degen2_aas, "Ile", "S_avg")))
aaUsage_full <- matrix(data=NA, nrow=length(genome_list$V1), ncol=length(aa_unique), dimnames=list(genome_list$V1, aa_unique))
aaUsage_HEG <- matrix(data=NA, nrow=length(genome_list$V1), ncol=length(aa_unique), dimnames=list(genome_list$V1, aa_unique))
aaUsage_rest <- matrix(data=NA, nrow=length(genome_list$V1), ncol=length(aa_unique), dimnames=list(genome_list$V1, aa_unique))

genomes_found <- c()
##CUB calculations
for(genome_oi in genome_list$V1){ #for each genome
  tryCatch({
    codCnt_matrix <- read.table(file = paste0("example_codon_data/",genome_oi,"/", genome_oi,"_cds_from_genomic.fna.codcnt_matrix"), header = T, sep = "\t", row.names = 1)
    ENCp_matrix <- read.table(file = paste0("example_codon_data/",genome_oi,"/", genome_oi,"_cds_from_genomic.fna.enc_matrix"), header = T, sep = "\t", row.names = 1)
    
    ENCp_threshold <- quantile(x = ENCp_matrix$Ncp, probs = ENCp_percentile, na.rm = T)
    
    #identify HEG genes with varying criteria
    #genes in lowest x percentile and with >99 codons
    HEG_rows_def_ENCp <- rownames(ENCp_matrix)[which(ENCp_matrix$Ncp<ENCp_threshold & ENCp_matrix$n_codons>99)]
    #genes annotated as ribosmoal genes, RNA polymerases, and elongation factors
    ribosomal_rows <- rownames(codCnt_matrix)[grepl("ribosomal protein", rownames(codCnt_matrix))]
    EF_rows <- rownames(codCnt_matrix)[grepl("elongation factor Tu", rownames(codCnt_matrix))]
    RNApol_rows <- rownames(codCnt_matrix)[grepl("\\[protein=DNA-directed RNA polymerase\\]", rownames(codCnt_matrix)) | (grepl("DNA-directed RNA polymerase", rownames(codCnt_matrix)) & grepl("subunit A'", rownames(codCnt_matrix)))  | (grepl("DNA-directed RNA polymerase", rownames(codCnt_matrix)) & grepl("subunit B", rownames(codCnt_matrix)))]
    HEG_rows_def_ribo <- c(ribosomal_rows, EF_rows, RNApol_rows)
    
    #definition of HEG
    if(ignore_ribo==TRUE){
      HEG_rows <- setdiff(HEG_rows_def_ENCp, HEG_rows_def_ribo)
      remaining_rows <- setdiff(rownames(codCnt_matrix), c(HEG_rows_def_ribo, HEG_rows))
    } else{
      HEG_rows <- HEG_rows_def_ribo
      remaining_rows <- setdiff(rownames(codCnt_matrix), HEG_rows_def_ribo)
    }
    
    # #calculate position specific base composition to derive expected frequencies;
    # #**base composition is derived from codon counts, pooled independently for HEG and remaining genes**
    # #for HEG
    baseCnt_HEG <- matrix(data = NA, nrow = 4, ncol = 3, dimnames = list(base_list, 1:3))
    baseFreq_HEG <- matrix(data = NA, nrow = 4, ncol = 3, dimnames = list(base_list, 1:3))
    codCnt_HEG <- colSums(x = codCnt_matrix[HEG_rows,], na.rm = T) #pool codon counts for all HEG
    codon_all <- colnames(codCnt_matrix)
    for(pos_oi in 1:3){
      for(base_oi in 1:4){
        baseCnt_HEG[base_oi, pos_oi] <- sum(codCnt_HEG[sapply(codon_all, function(x) {strsplit(x, split = "")[[1]][pos_oi]==base_list[base_oi]})])
      }
    }
    #correct for stop codons
    baseCnt_HEG["T",1] <- baseCnt_HEG["T",1] - sum(codCnt_HEG[c("TGA", "TAA", "TAG")])
    baseCnt_HEG["A",2] <- baseCnt_HEG["A",2] - sum(codCnt_HEG[c("TAA", "TAG")])
    baseCnt_HEG["G",2] <- baseCnt_HEG["G",2] - sum(codCnt_HEG[c("TGA")])
    baseCnt_HEG["A",3] <- baseCnt_HEG["A",3] - sum(codCnt_HEG[c("TGA", "TAA")])
    baseCnt_HEG["G",3] <- baseCnt_HEG["G",3] - sum(codCnt_HEG[c("TAG")])

    #calculate frequncies
    baseFreq_HEG[,1] <- round(baseCnt_HEG[,1]/colSums(baseCnt_HEG)[1], digits = 3)
    baseFreq_HEG[,2] <- round(baseCnt_HEG[,2]/colSums(baseCnt_HEG)[2], digits = 3)
    baseFreq_HEG[,3] <- round(baseCnt_HEG[,3]/colSums(baseCnt_HEG)[3], digits = 3)

    #for remaining genes
    baseCnt_rest <- matrix(data = NA, nrow = 4, ncol = 3, dimnames = list(base_list, 1:3))
    baseFreq_rest <- matrix(data = NA, nrow = 4, ncol = 3, dimnames = list(base_list, 1:3))
    codCnt_rest <- colSums(x = codCnt_matrix[remaining_rows,], na.rm = T) #pool codon counts for all remaining genes
    for(pos_oi in 1:3){
      for(base_oi in 1:4){
        baseCnt_rest[base_oi, pos_oi] <- sum(codCnt_rest[sapply(codon_all, function(x) {strsplit(x, split = "")[[1]][pos_oi]==base_list[base_oi]})])
      }
    }

    baseCnt_rest["T",1] <- baseCnt_rest["T",1] - sum(codCnt_rest[c("TGA", "TAA", "TAG")])
    baseCnt_rest["A",2] <- baseCnt_rest["A",2] - sum(codCnt_rest[c("TAA", "TAG")])
    baseCnt_rest["G",2] <- baseCnt_rest["G",2] - sum(codCnt_rest[c("TGA")])
    baseCnt_rest["A",3] <- baseCnt_rest["A",3] - sum(codCnt_rest[c("TGA", "TAA")])
    baseCnt_rest["G",3] <- baseCnt_rest["G",3] - sum(codCnt_rest[c("TAG")])

    baseFreq_rest[,1] <- round(baseCnt_rest[,1]/colSums(baseCnt_rest)[1], digits = 3)
    baseFreq_rest[,2] <- round(baseCnt_rest[,2]/colSums(baseCnt_rest)[2], digits = 3)
    baseFreq_rest[,3] <- round(baseCnt_rest[,3]/colSums(baseCnt_rest)[3], digits = 3)


    #calculate expected codon frequencies as products of position specific base frequencies
    calc_exp_codon_freq <- function(codon_oi, baseFreqs){
      b1 <- strsplit(codon_oi, "")[[1]][1]
      b2 <- strsplit(codon_oi, "")[[1]][2]
      b3 <- strsplit(codon_oi, "")[[1]][3]
      return(round(baseFreqs[b1,1]*baseFreqs[b2,2]*baseFreqs[b3,3], digits = 5))
    }

    codon_sense <- setdiff(codon_all, c("TGA", "TAA", "TAG", "ATG", "TGG"))
    exp_codFreq_HEG <- sapply(codon_sense, function(x) calc_exp_codon_freq(x,baseFreq_HEG))
    exp_codFreq_rest <- sapply(codon_sense, function(x) calc_exp_codon_freq(x,baseFreq_rest))

    #calculate amino acid specific ENC' following Novembre (2001)
    chiSq_HEG <- sapply(aa_unique, function(x) 0)
    F_HEG <- sapply(aa_unique, function(x) 0)
    chiSq_rest <- sapply(aa_unique, function(x) 0)
    F_rest <- sapply(aa_unique, function(x) 0)

    codNormFreq_HEG <- sapply(codon_sense, function(x) NA) #actual
    codNormFreq_rest <- sapply(codon_sense, function(x) NA)
    exp_codNormFreq_HEG <- sapply(codon_sense, function(x) NA) #expected
    exp_codNormFreq_rest <- sapply(codon_sense, function(x) NA)

    for(aa_oi in aa_unique){ #for each amino acid
      codons_oi <- codon_table[codon_table$aa==aa_oi, "codon"]
      degen_oi <- codon_table[codon_table$aa==aa_oi, "degen"][1]
      #get normalized codon frequencies within each amino acid
      sum_codCnt_HEG <- sum(codCnt_HEG[codons_oi])
      sum_codCnt_rest <- sum(codCnt_rest[codons_oi])
      sum_exp_codFreq_HEG <- sum(exp_codFreq_HEG[codons_oi])
      sum_exp_codFreq_rest <- sum(exp_codFreq_rest[codons_oi])
      for(codon_oi in codons_oi){
        codNormFreq_HEG[codon_oi] <- codCnt_HEG[codon_oi]/sum_codCnt_HEG
        codNormFreq_rest[codon_oi] <- codCnt_rest[codon_oi]/sum_codCnt_rest
        exp_codNormFreq_HEG[codon_oi] <- exp_codFreq_HEG[codon_oi]/sum_exp_codFreq_HEG
        exp_codNormFreq_rest[codon_oi] <- exp_codFreq_rest[codon_oi]/sum_exp_codFreq_rest
      }
      #get chi squared and F values
      for(codon_oi in codons_oi){
        chiSq_HEG[aa_oi] <- chiSq_HEG[aa_oi] + (codNormFreq_HEG[codon_oi] - exp_codNormFreq_HEG[codon_oi])^2/(exp_codNormFreq_HEG[codon_oi])
        chiSq_rest[aa_oi] <- chiSq_rest[aa_oi] + (codNormFreq_rest[codon_oi] - exp_codNormFreq_rest[codon_oi])^2/(exp_codNormFreq_rest[codon_oi])
      }
      chiSq_HEG[aa_oi] <- round(chiSq_HEG[aa_oi]*sum(codCnt_HEG[codons_oi]), digits = 4)
      chiSq_rest[aa_oi] <- round(chiSq_rest[aa_oi]*sum(codCnt_rest[codons_oi]), digits = 4)
      F_HEG[aa_oi] <- round((chiSq_HEG[aa_oi] + sum_codCnt_HEG - degen_oi)/(degen_oi*(sum_codCnt_HEG-1)), digits = 4)
      F_rest[aa_oi] <- round((chiSq_rest[aa_oi] + sum_codCnt_rest - degen_oi)/(degen_oi*(sum_codCnt_rest-1)), digits = 4)
    }

    aawise_ENCp_HEG_all[genome_oi,] <- round(1/F_HEG, digits = 4)
    aawise_ENCp_rest_all[genome_oi,] <- round(1/F_rest, digits = 4)

    #average across amino acid degeneracy classes, following Wright (1990), Novembre (2001)
    ENCp_HEG <- sum(2,9/mean(F_HEG[degen2_aas]),1/F_HEG["Ile"],5/mean(F_HEG[degen4_aas]),3/mean(F_HEG[degen6_aas]), na.rm = T)
    ENCp_rest <- sum(2,9/mean(F_rest[degen2_aas]),1/F_rest["Ile"],5/mean(F_rest[degen4_aas]),3/mean(F_rest[degen6_aas]), na.rm = T)

    ENCp_HEG <- round(ENCp_HEG, digits = 4)
    ENCp_rest <- round(ENCp_rest, digits = 4)
    ENCp_diff <- round((ENCp_rest-ENCp_HEG)/ENCp_rest,digits=4)
    ENCp_all[genome_oi, ] <- c(ENCp_HEG, ENCp_rest, ENCp_diff)

    ##calculate S
    S_CUB  <- sapply(c(degen2_aas, "Ile", "avg"), function(x) 0)
    for(aa_oi in c(degen2_aas, "Ile")){
      codons_oi <- codon_table[codon_table$aa==aa_oi, "codon"]
      S_CUB[aa_oi] <- log((codCnt_HEG[codons_oi[2]]/codCnt_HEG[codons_oi[1]])/(codCnt_rest[codons_oi[2]]/codCnt_rest[codons_oi[1]]))
    }

    #calculate weighted average of S, following Sharp et al (2005)
    aas_oi <- c("Phe", "Ile", "Tyr", "Asn")
    S_weights <- sapply(aas_oi, function(x) sum(codCnt_HEG[codon_table[codon_table$aa==x, "codon"]][1:2]))
    S_CUB["avg"] <- weighted.mean(x = S_CUB[aas_oi], w = S_weights, na.rm = T)
    
    aawise_S_all[genome_oi,] <- round(S_CUB, digits = 4)
    
  	# #calculate amino acid usage
  	aaCnt_HEG <- sapply(aa_unique, function(x) sum(codCnt_HEG[codons_aa[[x]]], na.rm = T))
  	aaUsage_HEG[genome_oi, ] <- t(round(aaCnt_HEG/sum(aaCnt_HEG), digits=4))
  	aaCnt_rest <- sapply(aa_unique, function(x) sum(codCnt_rest[codons_aa[[x]]], na.rm = T))
  	aaUsage_rest[genome_oi, ] <- t(round(aaCnt_rest/sum(aaCnt_rest), digits=4))
  	aaCnt_full <- aaCnt_HEG+aaCnt_rest
  	aaUsage_full[genome_oi, ] <- t(round(aaCnt_full/sum(aaCnt_full),digits=4))
  	
  	genomes_found <- c(genomes_found, genome_oi)
    }, error=function(e){
            print(e)
    })
}

if(ignore_ribo==TRUE){
	file_suffix <- paste0("_p", ENCp_percentile)
}else{
	file_suffix <- ""
}
write.table(aawise_ENCp_HEG_all, file = paste0("traits_data/ENCp_HEG_aawise", file_suffix,".txt"))
write.table(aawise_ENCp_rest_all, file = paste0("traits_data/ENCp_rest_aawise", file_suffix, ".txt"))
write.table(ENCp_all, file = paste0("traits_data/ENCp", file_suffix, ".txt"))

aawise_ENCp_diff <- (aawise_ENCp_rest_all-aawise_ENCp_HEG_all)/aawise_ENCp_rest_all
aawise_ENCp_diff <- round(aawise_ENCp_diff, digits = 4)
write.table(aawise_ENCp_diff, file = paste0("traits_data/deltaENCp_aawise", file_suffix, ".txt"))

write.table(aawise_S_all, file = paste0("traits_data/S_aawise", file_suffix, ".txt"))

write.table(aaUsage_HEG, file = paste0("traits_data/aaUsage_HEG", file_suffix, ".txt"))
write.table(aaUsage_rest, file = paste0("traits_data/aaUsage_others", file_suffix, ".txt"))
write.table(aaUsage_full, file = paste0("traits_data/aaUsage_full", file_suffix, ".txt"))