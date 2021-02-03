#!/bin/bash
. /u/local/Modules/default/init/modules.sh
#source /u/home/b/briscoel/project-ngarud/miniconda2/bin/activate /u/home/b/briscoel/project-ngarud/miniconda2/envs/python3
module load python/anaconda3

while read -r arg_1 arg_2 arg_3 arg_4 arg_5 arg_6 arg_7 arg_8 arg_9 arg_10 arg_11 arg_12 arg_13 arg_14 arg_15 arg_16 arg_17 arg_18 arg_19 arg_20 arg_21; do
        python paraMINERVA_test_train_prediction.py $arg_1 $arg_2 $arg_3 $arg_4 $arg_5 $arg_6 $arg_7 $arg_8 $arg_9 $arg_10 $arg_11 $arg_12 $arg_13 $arg_14 $arg_15 $arg_16 $arg_17 $arg_18 $arg_19 $arg_20 $arg_21;
done < pinputs/data_$SGE_TASK_ID.in
