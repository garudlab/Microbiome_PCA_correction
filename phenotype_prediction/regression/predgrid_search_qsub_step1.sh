#!/bin/bash

# example 


#python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch CRC_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 dataaug 3 1 1 CRC 100 entropy 1 0.1 5 1 0 0 study 1
#Domain correct

# ./predgrid_search_qsub_step1.sh 1 AGP_complete_otu raw none raw reg 5234847

#./predgrid_search_qsub_step1.sh 1 AGP_complete_otu raw clr_scale minervaclrscale reg
# ./predgrid_search_qsub_step1.sh 51 AGP_complete_otu limma none limma reg
# ./predgrid_search_qsub_step1.sh 101 AGP_complete_otu bmc none bmc reg
#./predgrid_search_qsub_step1.sh 151 AGP_complete_otu ComBat logscale ComBatlogscale reg

# ./predgrid_search_qsub_step1.sh 51 AGP_complete_otu raw none domaincorr
# ./predgrid_search_qsub_step1.sh 101 AGP_max_k6 raw none domaincorr
# ./predgrid_search_qsub_step1.sh 152 AGP_max_k5 raw none domaincorr
# ./predgrid_search_qsub_step1.sh 152 AGP_max_k7 raw none raw

first_count_input=$1
dataset_input=$2
method_input=$3

trans_input=$4
name_input=$5

enet_input=$6

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
do
	if [[ "$enet_input" == *"enet"* ]]; then
		for l1ratio in 0 .1 .5 .7 .9 .95 .99 1; 
			do for alpha in 0.025 0.05 .125 .25 .5 1 2 4;
				do 
					COUNTER=$((COUNTER + 1)); 
					echo $COUNTER; 
				
					if [[ "$dataset_input" == *"AGP_complete"* ]]; then
						echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $trainit 1 $alpha $l1ratio Instrument" > pinputs/data_$COUNTER.in; 
					fi

					if [[ "$dataset_input" == *"AGP_max"* ]]; then
						echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $trainit 1 $alpha $l1ratio Instrument" > pinputs/data_$COUNTER.in; 
					fi
				done; 
			done; 
	else
		COUNTER=$((COUNTER + 1)); 
		echo $COUNTER; 
		if [[ "$dataset_input" == *"AGP_complete"* ]]; then
			echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $trainit 0 0 0 Instrument" > pinputs/data_$COUNTER.in; 
		fi

		if [[ "$dataset_input" == *"AGP_max"* ]]; then
			echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $trainit 0 0 0 Instrument" > pinputs/data_$COUNTER.in; 
		fi
	fi
done;



if [[ "$dataset_input" == *"AGP_max"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=16G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_prediction.sh"

elif [[ "$dataset_input" == *"AGP_complete"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=16G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_prediction.sh"

else
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=5G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_prediction.sh"
fi



