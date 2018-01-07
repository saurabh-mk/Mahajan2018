# Mahajan2018
#### Code and data corresponding to a paper on amino acid specific tRNA gene copy numbers and codon usage bias.

This repository contains-
1. raw data
1. code to process raw data to data directly used to generate figures
1. data used to generate figures, and
1. code to generate figures from processed data

* **raw data**
    1. [lists](https://github.com/saurabh-mk/Mahajan2018/tree/master/genome_lists) of all (and clade specific) genomes used in the study
    1. [phylogeny](https://github.com/saurabh-mk/Mahajan2018/tree/master/trees) used in the study
    1. details of all genomes from IMG database ([here](https://github.com/saurabh-mk/Mahajan2018/blob/master/traits_data/all_genomes_IMG_data.txt))
    1. rrnDB v5.1 (a database of rRNA copy numbers)([here](https://github.com/saurabh-mk/Mahajan2018/blob/master/traits_data/rrnDB-5.1.tsv))
    1. tRNA gene copy number data ([here](https://github.com/saurabh-mk/Mahajan2018/blob/master/traits_data/tRNA_GCN.txt))
    1. codon specific translation time for _E. coli_ and _B. subtilis_ ([here](https://github.com/saurabh-mk/Mahajan2018/blob/master/traits_data/decoding_times.tab))
    1. a codon table file ([here](https://github.com/saurabh-mk/Mahajan2018/blob/master/codon_aa_map.txt))
    
* **code to process raw data**
    1. a [script](https://github.com/saurabh-mk/Mahajan2018/blob/master/codon_count_calculations_script.sh) to calculate _genewise_ codon counts, base composition, and codon usage metrics like ENC', for multiple genomes
       1. accessory python script to split multifasta files
       1. accessory R script to tabulate codon usage data 
       1. sample output
    1. a [script](https://github.com/saurabh-mk/Mahajan2018/blob/master/calculate_ENCp_S.R) to calculate _genomewise_ amino acid specific CUB metrics (&#x394;ENC', S), amino acid usage for all genomes
    1. a [script](https://github.com/saurabh-mk/Mahajan2018/blob/master/obtain_rRNA_CN.R) to obtain 'consensus' rRNA copy numbers from IMG data and rrnDB
    
* **data used to generate figures**
In all data files, first row contains headers (one less than columns); first column contains genome ids (refseq/genbak)
    1. genomewise rRNA copy numbers
    1. genomewise anticodon specific tRNA gene copy numbers
    1. genomewise amino acid specific CUB: amino acid averaged &#x394;ENC', amino acid specific &#x394;ENC' and S
    1. genomewise amino acid specific CUB: amino acid averaged &#x394;ENC', amino acid specific &#x394;ENC' and S (using alternate definition of highly expressed genes)
    1. genomewise amino acid usage in highly expressed genes
    
* **code used to generate figures**
    1. Figure 1, S1, S2, and S3
    1. Figure 2, 3, and 4
    1. Figure 5, 6, S4, S5, S6, and S7
    1. Figure 7
    1. Figure S8

If you have any trouble, please write to me at saurabh.mk@gmail.com.
