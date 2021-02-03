#!/bin/bash

#examples
# qsub -cwd -V -N CombatotulodoTOMvalasmb -l h_data=5G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu ComBatfilter_TRUE_trans_logscale BatchCorrected bin_crc_normal 0 0 10 10 ComBatlogscale_LODO 1 1 0 1 CRC 1 dataset_name 1"			
# qsub -cwd -V -N CombatotulodoTOMvalasmb -l h_data=5G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu ComBat_with_batch2filter_TRUE_trans_logscale BatchCorrected bin_crc_normal 0 0 10 10 ComBatBatch2logscale_LODO 1 1 0 1 CRC 1 dataset_name 1"			

# ./grid_search_qsub_step2.sh CRC_thomas_otu ComBat logscale ComBatlogscale_LODO val
# ./grid_search_qsub_step2.sh CRC_thomas_otu ComBat logscale ComBatlogscale val

# ./grid_search_qsub_step2.sh CRC_thomas_otu limma none limma_LODO val
# ./grid_search_qsub_step2.sh CRC_thomas_otu bmc none bmc_LODO val
# ./grid_search_qsub_step2.sh CRC_thomas_otu raw none dataaug_LODO val
#./grid_search_qsub_step2.sh CRC_thomas_otu raw clr_scale dataaug_clrscale val

# ./grid_search_qsub_step2.sh CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale_LODO val
# ./grid_search_qsub_step2.sh CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale val
# ./grid_search_qsub_step2.sh Thomas_k7 raw none raw val


# for k in 6 7; do for method in combat; do ./grid_search_qsub_step2.sh CRC_k"$k" $method none $method val; done; done

# ./grid_search_qsub_step2.sh Thomas_k7 raw clr_scale minervaclrscaleLODO val

#./grid_search_qsub_step2.sh "Thomas_k6&Thomas_k7" raw clr_scale minervaclrscale val  #live run
#./grid_search_qsub_step2.sh "Thomas_k6&Thomas_k7" raw none raw val  #live run
#./grid_search_qsub_step2.sh "AGP_complete_otu" raw none raw val
#./grid_search_qsub_step2.sh "AGP_complete_otu" raw clrscale minervaclrscale val

#./grid_search_qsub_step2.sh "CRC_otu" raw none raw val
#./grid_search_qsub_step2.sh "CRC_k6&CRC_k7" raw none raw val
# ./grid_search_qsub_step2.sh CRC_k6&CRC_k7 limma none raw val

#./grid_search_qsub_step2.sh "CRC_k6&CRC_k7" raw none raw val
#./grid_search_qsub_step2.sh "CRC_k7" ComBat none combat val


# ./grid_search_qsub_step2.sh CRC_otu raw none dataaug val
# ./grid_search_qsub_step2.sh CRC_thomas_otu limma none limma val
# ./grid_search_qsub_step2.sh CRC_thomas_otu bmc none bmc val
# ./grid_search_qsub_step2.sh CRC_k7 raw none domaincorr val #5406153 
# ./grid_search_qsub_step2.sh Thomas_k7 raw none domaincorr val
# ./grid_search_qsub_step2.sh AGP_complete_otu raw none dataaug val
# ./grid_search_qsub_step2.sh AGP_complete_otu raw none domaincorr val
# ./grid_search_qsub_step2.sh AGP_max_k7 raw none domaincorr

# no assemblies run yet

# ./grid_search_qsub_step2.sh CRC_k7 raw none dataaug val

# ./grid_search_qsub_step2.sh CRC_k7 raw none dataaug val
# ./grid_search_qsub_step2.sh Thomas_k7 raw none dataaug val


# ./grid_search_qsub_step2.sh CRC_k7 raw none raw val

# Run OTUs for Complete AGP
# Check that all the bc types are there - > running raw_clr_scale, 

# ./grid_search_qsub_step2.sh CRC_otu limma none limma val
# ./grid_search_qsub_step2.sh CRC_otu ComBat logscale ComBatlogscale val
# ./grid_search_qsub_step2.sh CRC_otu bmc none bmc val
# ./grid_search_qsub_step2.sh AGP_complete_otu ComBat logscale ComBat val
#./grid_search_qsub_step2.sh AGP_complete_otu raw clr_scale minervaclrscale val

# ./grid_search_qsub_step2.sh AGP_max_k7 raw none domaincorr val

live_run=0
dataset_input=$1
echo $dataset_input
method_input=$2

trans_input=$3
name_input=$4
val_input=$5

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
		


#otu tom
if [[ "$dataset_input" == *"thomas"* ]]; then
	qsub -cwd -V -N otulodoTOMvalasmb -l h_data=5G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input $use_val 0 1 CRC $use_lodo dataset_name 1"			
fi

if [[ "$dataset_input" == *"Thomas"* ]]; then
	if [ $live_run == 1 ]; then
		./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $use_minerva $use_val 1 1 CRC $use_lodo study 1
	else
		qsub -cwd -V -N klodoTOMvalasmb -l h_data=5G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch "$dataset_input" "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input $use_val 1 1 CRC $use_lodo study 1"			
	fi
fi

if [[ "$dataset_input" == *"CRC_k"*  ]]; then
	if [ $live_run == 1 ]; then
		./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $use_minerva $use_val 0 1 1 $use_lodo study 1		
	else
		qsub -cwd -V -N klodoTOMvalasmb -l h_data=5G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch "$dataset_input" "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input $use_val 0 1 1 $use_lodo study 1"			
	
	fi

fi

if [[ "$dataset_input" == *"CRC_otu"*  ]]; then
	qsub -cwd -V -N "$name_input"_asmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input $use_val 1 1 CRC $use_lodo study 1"			
fi

if [[ "$dataset_input" == *"AGP_"*  ]]; then
	qsub -cwd -V -N "$name_input"_asmb -l h_data=15G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_antibiotic_last_year 0 0 10 10 $name_input $minerva_input $use_val 1 1 Yes $use_lodo Instrument 1"			
fi

