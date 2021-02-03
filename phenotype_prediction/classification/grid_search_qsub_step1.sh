#!/bin/bash

# example 
# ./grid_search_qsub_step1.sh 1 CRC_thomas_otu ComBat logscale ComBatlogscale_LODO
# ./grid_search_qsub_step1.sh 2701 CRC_thomas_otu ComBat logscale ComBatlogscale
# ./grid_search_qsub_step1.sh 5401 CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale_LODO
# ./grid_search_qsub_step1.sh 8101 CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale
# ./grid_search_qsub_step1.sh 1 CRC_thomas_otu raw none dataaug_LODO
# ./grid_search_qsub_step1.sh 2701 CRC_thomas_otu limma none limma_LODO
# ./grid_search_qsub_step1.sh 5401 CRC_thomas_otu bmc none bmc_LODO




# ./grid_search_qsub_step1.sh 10801 Thomas_k7 ComBat logscale ComBatlogscale_LODO
# ./grid_search_qsub_step1.sh 13501 Thomas_k7 ComBat logscale ComBatlogscale

# ./grid_search_qsub_step1.sh 10801 CRC_otu bmc none bmc
# ./grid_search_qsub_step1.sh 13501 CRC_otu limma none limma
# ./grid_search_qsub_step1.sh 16201 CRC_otu ComBat logscale ComBatlogscale

# ./grid_search_qsub_step1.sh 5401 CRC_thomas_otu limma none limma
# ./grid_search_qsub_step1.sh 8101 CRC_thomas_otu bmc none bmc


# ./grid_search_qsub_step1.sh 16201 CRC_otu raw none raw
# ./grid_search_qsub_step1.sh 2701 CRC_otu raw none raw_LODO
# ./grid_search_qsub_step1.sh 1 CRC_otu raw clrscale minervaclrscale
# ./grid_search_qsub_step1.sh 2701 CRC_otu raw clrscale minervaclrscale_LODO

# ./grid_search_qsub_step1.sh 1 AGP_complete_otu raw none raw
# ./grid_search_qsub_step1.sh 5401 AGP_complete_otu raw clrscale minervaclrscale
# ./grid_search_qsub_step1.sh 5401 AGP_complete_otu raw clrscale minervaclrscale LODO

#Data augment

# ./grid_search_qsub_step1.sh 1 CRC_k7 raw none dataaug
# ./grid_search_qsub_step1.sh 2701 Thomas_k7 raw none dataaug
# ./grid_search_qsub_step1.sh 5401 AGP_max_k7 raw none dataaug 

# ./grid_search_qsub_step1.sh 8101 CRC_otu raw clr_scale dataaug_clrscale - RAN
# ./grid_search_qsub_step1.sh 10801 CRC_thomas_otu raw clr_scale dataaug_clrscale - RUN
# ./grid_search_qsub_step1.sh 1 AGP_complete_otu raw clr_scale minervaclrscale
# ./grid_search_qsub_step1.sh 2700 AGP_complete_otu bmc none bmc
# ./grid_search_qsub_step1.sh 5400 AGP_complete_otu ComBat logscale ComBat
# ./grid_search_qsub_step1.sh 8101 AGP_complete_otu limma none limma




#python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch CRC_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 dataaug 3 1 1 CRC 100 entropy 1 0.1 5 1 0 0 study 1
#Domain correct
# ./grid_search_qsub_step1.sh 5401 CRC_otu raw none domaincorr
# ./grid_search_qsub_step1.sh 8101 CRC_thomas_otu raw none domaincorr
# ./grid_search_qsub_step1.sh 13501 AGP_complete_otu raw none domaincorr

# ./grid_search_qsub_step1.sh 1 CRC_k7 raw none domaincorr
# ./grid_search_qsub_step1.sh 2700 Thomas_k7 raw none domaincorr
# ./grid_search_qsub_step1.sh 1 AGP_max_k6 raw none domaincorr

#python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch CRC_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 domaincorr 4 1 1 CRC 100 entropy 1 0.1 5 1 0 0 study 1
#qsub -cwd -V -N dataaug -l h_data=7G,time=24:00:00 -b y -t 10848:13500 "./run_array_paraMINERVA_test_train_grid.sh" (first part took 18 hours to do 45 jobs sad)
#qsub -cwd -V -N domaincorr -l h_data=7G,time=24:00:00 -b y -t 13502:16200 "./run_array_paraMINERVA_test_train_grid.sh" (first part took 18 hours to do 45 jobs sad)

first_count_input=$1
dataset_input=$2
method_input=$3

trans_input=$4
name_input=$5

if [[ "$name_input" == *"LODO"* ]]; then
	use_lodo=1
else
	use_lodo=0
fi

if [[ "$name_input" == *"minerva"* ]]; then
	minerva_input=1
elif [[ "$name_input" == *"dataaug"* ]]; then
	minerva_input=3
elif [[ "$name_input" == *"domaincorr"* ]]; then
	minerva_input=4
else
	minerva_input=0
fi

COUNTER=$first_count_input-1
echo $COUNTER
for trainit in {0..49}; 
	do for nest in 100 1000 1500; 
		do for crit in entropy gini; 
			do for leaf in 1 5 10; 
				do for feat in 0.1 0.3 0.5; 
					do 
						COUNTER=$((COUNTER + 1)); 
						echo $COUNTER; 
						if [[ "$dataset_input" == *"thomas"* ]]; then
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit $use_lodo dataset_name 1" > inputs/data_$COUNTER.in; 	
						fi
						if [[ "$dataset_input" == *"Thomas"* ]]; then
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit $use_lodo study 1" > inputs/data_$COUNTER.in; 
						fi
						if [[ "$dataset_input" == *"CRC_otu"* ]]; then
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit $use_lodo study 1" > inputs/data_$COUNTER.in; 
						fi

						if [[ "$dataset_input" == *"CRC_k"* ]]; then
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input 0 1 1 $nest $crit $leaf $feat 5 1 $trainit $use_lodo study 1" > inputs/data_$COUNTER.in; 
						fi

						if [[ "$dataset_input" == *"AGP_complete"* ]]; then
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_antibiotic_last_year 0 0 10 10 $name_input $minerva_input 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit $use_lodo Instrument 1" > inputs/data_$COUNTER.in; 
						fi

						if [[ "$dataset_input" == *"AGP_max"* ]]; then
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_antibiotic_last_year 0 0 10 10 $name_input $minerva_input 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit $use_lodo Instrument 1" > inputs/data_$COUNTER.in; 
						fi


					done; 
				done; 
			done; 
		done; 
	done;


if [[ "$dataset_input" == *"AGP_max"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=15G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"

elif [[ "$dataset_input" == *"AGP_complete"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=16G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"

else
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=5G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"
fi

#qsub -cwd -V -N bmc -l h_data=10G,time=24:00:00 -b y -t 2703:5399 "./run_array_paraMINERVA_test_train_grid.sh"
# qsub -cwd -V -N ComBatgrid -l h_data=16G,time=24:00:00 -b y -t 5507:8099 "./run_array_paraMINERVA_test_train_grid.sh"



