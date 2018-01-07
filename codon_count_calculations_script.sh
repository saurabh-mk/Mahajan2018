## Author: Saurabh Mahajan
## Created: Dec 30, 2017
## Description: Calculates genewise codon counts and ENC' using ENCprime program (Novembre 2001). SeqCount fails on large fasta files, so this script breaks it into chunks of user defined number of sequences and runs SeqCount, ENCprime, and pools back data into single files.
## Input arguments: file with list of genome ids (same as cds file labels); batch size i.e. how may cds to put in each sub-file
## Depends on: ENCprime (https://github.com/NovembreLab/ENCprime); split_multifasta.py script; an R script (tabulate_ENCprime_output.R) to parse data into neat tables. Before compiling, the ENCprime.c file had to be modified to increase length of output filename on line 127 from 40 to 200.
## Run as: "parallel -j# -a example_genome_data/required_genomes_list.txt ./codon_count_calculation_script.sh ::: 5"; # is number of processors available, 5 is batch size

genome=$1
batch_size=$2

ENCprime_dir="" #main directory for ENCprime
data_dir="" #where cds_from_genomic.fna files are stored
output_dir="" #where output data is to be stored
code_dir="" #where split_multifasta.py and tabulate_ENCprime_output.R files are saved

#create specific output directory for the genome
mkdir $output_dir/$genome
sp_output_dir=$output_dir/$genome

#count the number of cds in the cds file, each starts with a ">"
n_cds=$(grep -c ">" $data_dir/${genome}_cds_from_genomic.fna)

# split the cds_from_genomic_file into chunks
python $code_dir/split_multifasta.py $genome $data_dir $sp_output_dir $batch_size
remaining_records=$(($n_cds % $batch_size))

#run SeqCount program of ENCprime
if [ "$remaining_records" -ne 0 ]; then # if number of cds is not multiple of batch size
	n_files=$(( ($n_cds / $batch_size) + 1))
	for batch_no in $(seq $(($n_files - 1))); do
		echo "running SeqCount for batch no " $batch_no of $n_files
		$ENCprime_dir/bin/SeqCount -c $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna $batch_size
		$ENCprime_dir/bin/SeqCount -n $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna $batch_size
		#rm $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna #remove the file, its job is done
	done

	echo "running SeqCount for batch no " $n_files of $n_files
	records_end=$(($n_cds % $batch_size)) #number of remaining records, need to pass to SeqCounts

	$ENCprime_dir/bin/SeqCount -c $sp_output_dir/${genome}_cds_from_genomic_$n_files.fna $records_end
	$ENCprime_dir/bin/SeqCount -n $sp_output_dir/${genome}_cds_from_genomic_$n_files.fna $records_end
	#rm $sp_output_dir/${genome}_cds_from_genomic_$n_files.fna #remove the file, its job is done

else # when number of cds is multiple of batch size
	n_files=$(($n_cds / $batch_size))
	for batch_no in $(seq $n_files); do
		echo "running SeqCount for batch no " $batch_no of $n_files
		$ENCprime_dir/bin/SeqCount -c $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna $batch_size
		$ENCprime_dir/bin/SeqCount -n $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna $batch_size
		#rm $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna #remove the file, its job is done
	done
fi

##run the ENCprime program
for batch_no in $(seq $n_files); do
	echo "running ENCprime for batch no " $batch_no of $n_files
	$ENCprime_dir/bin/ENCprime $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.codcnt $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.acgtfreq 1 $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.encp 0 -q

	#pasting data for all sub-files together; and removing batch files
	tail -n +2 -q $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.encp | head -n -1 >> $sp_output_dir/${genome}_cds_from_genomic.fna.encp
	#rm $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.encp
	tail -n +2 -q $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.acgtfreq | head -n -1 >> $sp_output_dir/${genome}_cds_from_genomic.fna.acgtfreq
	#rm $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.acgtfreq
	tail -n +2 -q $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.acgtcnt | head -n -1 >> $sp_output_dir/${genome}_cds_from_genomic.fna.acgtcnt
	#rm $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.acgtcnt
	tail -n +4 -q $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.codfreq | head -n -1 >> $sp_output_dir/${genome}_cds_from_genomic.fna.codfreq
	#rm $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.codfreq
	tail -n +4 -q $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.codcnt | head -n -1 >> $sp_output_dir/${genome}_cds_from_genomic.fna.codcnt
	#rm $sp_output_dir/${genome}_cds_from_genomic_$batch_no.fna.codcnt

done

#tabulate the output to make it readable by subsequent programs
Rscript --vanilla tabulate_ENCprime_output.R ${sp_output_dir} ${genome}
