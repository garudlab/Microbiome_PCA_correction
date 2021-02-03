#!/bin/bash
. /u/local/Modules/default/init/modules.sh
#source /u/home/b/briscoel/project-ngarud/miniconda2/bin/activate /u/home/b/briscoel/project-ngarud/miniconda2/envs/python3
module load python/anaconda3
python MINERVA_test_train_prediction.py $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11} ${12} ${13} ${14}

