#!/bin/bash

# examples
# ./batch_correction_qsub.sh CRC_otu limma none
# ./batch_correction_qsub.sh CRC_thomas_otu DomainCorrect none
# ./batch_correction_qsub.sh CRC_thomas_otu bmc none
#./batch_correction_qsub.sh CRC_thomas_otu ComBat_with_biocovariates_with_seqbatch logscale
#./batch_correction_qsub.sh CRC_otu raw none
#./batch_correction_qsub.sh CRC_k raw none
#./batch_correction_qsub.sh CRC_k7 raw clrscale

#./batch_correction_qsub.sh AGP_complete_otu raw none
#./batch_correction_qsub.sh AGP_complete_otu raw clr_scale
# ./batch_correction_qsub.sh AGP_complete_otu limma none
# ./batch_correction_qsub.sh AGP_complete_otu ComBat logscale
# ./batch_correction_qsub.sh AGP_complete_otu bmc none
# ./batch_correction_qsub.sh AGP) ComBat_with_biocovariates_with_seqbatch logscale
# ./batch_correction_qsub.sh AGP_complete_otu raw none
# ./batch_correction_qsub.sh AGP_complete_otu raw clr_scale
# for SV
# # ./batch_correction_qsub.sh AGP_complete_otu raw none
# ./batch_correction_qsub.sh CRC_otu raw none
# ./batch_correction_qsub.sh CRC_thomas_otu raw none# 
# ./batch_correction_qsub.sh AGP_max_k7 raw none # 5327228
# ./batch_correction_qsub.sh CRC_k7 raw none  # 5325053 
# ./batch_correction_qsub.sh Thomas_k7 raw none
# ./batch_correction_qsub.sh Hispanic_k7 raw none

# ./batch_correction_qsub.sh CRC_k7 DomainCorrect none
# ./batch_correction_qsub.sh CRC_k7 DomainCorrect clr_scale
# ./batch_correction_qsub.sh AGP_max_k7 DomainCorrect clr_scale # 5327228
# ./batch_correction_qsub.sh AGP_max_k7 DomainCorrect none # 5327228

dataset_input=$1
method_input=$2
trans_input=$3

echo $dataset_input

if [ "$dataset_input" == "CRC_thomas_otu" ]
then
	echo $dataset_input
	for method in $method_input; 
		do for phen in bin_crc_normal; 
			do for tran in $trans_input; 
				do for sv in 10; 
					do for k in 6; 
						do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; 
					done; 
				done; 
			done; 
		done; 
	done
fi

if [[ "$dataset_input" == *"Thomas"* ]]
then
	echo $dataset_input
	for method in $method_input; 
		do for phen in bin_crc_normal; 
			do for tran in $trans_input; 
				do for sv in 10; 
					do for k in 7; 
						do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; 
					done; 
				done; 
			done; 
		done; 
	done
fi

if [ "$dataset_input" == "CRC_otu" ]
then
	echo $dataset_input
	for method in $method_input; 
		do for phen in bin_crc_normal; 
			do for tran in $trans_input; 
				do for sv in 1; 
					do for k in 6; 
						do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 0 -arg 0; 
					done; 
				done; 
			done; 
		done; 
	done
fi

if [[ "$dataset_input" == *"CRC_k"* ]]
then
	echo $dataset_input
	for method in $method_input; 
		do for phen in bin_crc_normal; 
			do for tran in $trans_input; 
				do for sv in 10; 
					do for k in 7; 
						do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 0 -arg 0; 
					done; 
				done; 
			done; 
		done; 
	done
fi




if [[ "$dataset_input" == *"AGP"* ]]
then
	echo "AGP"
	for method in $method_input; 
		do for phen in "bin_antibiotic_last_year"; #"bmi_corrected"; # #"bmi_corrected"; #bin_antibiotic_last_year
			do for tran in $trans_input; 
				do for sv in 1; 
					do
					if [[ "$dataset_input" == *"_k"* ]]
					then
						for k in 7; 
						do
							echo "KMER"
							/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 15 -t 24 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg "Instrument" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; 
						done;
					else 
						echo "OTU"
						/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 15 -t 24 -hp -v 3.6.0 -arg otu -arg 5 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_complete -arg "$method" -arg $sv -arg "Instrument" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; 
					fi
					done;
				done; 
			done; 
		done; 
fi


if [[ "$dataset_input" == *"Hispanic"* ]]
then
	echo "hispanic"
	for method in $method_input; 
		do for phen in "bin_antibiotic"; #"bmi_corrected"; #bin_antibiotic_last_year
			do for tran in $trans_input; 
				do for sv in 1; 
					do
					if [[ "$dataset_input" == *"_k"* ]]
					then
						for k in 7; 
						do
							echo "KMER"
							/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 15 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg "mastermix_lot..exp." -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; 
						done;
					else 
						echo "OTU"
						/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 15 -t 24 -v 3.6.0 -arg otu -arg 5 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg "mastermix_lot..exp." -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; 
					fi
					done;
				done; 
			done; 
		done; 
fi




