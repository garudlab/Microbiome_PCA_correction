#!/bin/bash

#examples

# ./predgrid_search_qsub_step2.sh AGP_max_k5 ComBat none ComBat val reg
# ./predgrid_search_qsub_step2.sh AGP_max_k5 bmc none bmc val reg

# ./predgrid_search_qsub_step2.sh AGP_complete_otu raw none raw val reg
# ./predgrid_search_qsub_step2.sh AGP_complete_otu ComBat logscale ComBatlogscale val reg
# ./predgrid_search_qsub_step2.sh AGP_complete_otu raw clr_scale minervaclrscale val reg
# ./predgrid_search_qsub_step2.sh AGP_complete_otu bmc none bmc val reg
# ./predgrid_search_qsub_step2.sh AGP_complete_otu limma none limma val reg

# ./predgrid_search_qsub_step2.sh AGP_complete_otu raw none domaincorr val reg
# ./predgrid_search_qsub_step2.sh AGP_max_k6 raw none domaincorr val reg
./predgrid_search_qsub_step2.sh AGP_max_k6 raw none raw val reg
./predgrid_search_qsub_step2.sh AGP_max_k5 raw none raw val reg
# ./predgrid_search_qsub_step2.sh AGP_max_k7 raw none raw val reg

live_run=0
dataset_input=$1
echo $dataset_input
method_input=$2

trans_input=$3
name_input=$4
val_input=$5
enet_input=$6

if [[ "$name_input" == *"LODO"* ]]; then
	use_lodo=1
else
	use_lodo=0
fi

if [[ "$val_input" == *"val"* ]]; then
	use_val=1
else
	use_val=0
fi

echo $name_input
if [[ "$name_input" == *"minerva"* ]]; then
	minerva_input=1
elif [[ "$name_input" == *"dataaug"* ]]; then
	minerva_input=3
elif [[ "$name_input" == *"domaincorr"* ]]; then
	minerva_input=4
else
	minerva_input=0
fi
		
if [[ "$enet_input" == *"enet"* ]]; then
	use_enet=1
else
	use_enet=0
fi

#otu tom

if [[ "$dataset_input" == *"AGP_complete"*  ]]; then
	qsub -cwd -V -N "$name_input"_asmb -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $use_enet Instrument"			
fi

if [[ "$dataset_input" == *"AGP_max"*  ]]; then
	qsub -cwd -V -N "$name_input"_asmb -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $use_enet Instrument"			
fi
