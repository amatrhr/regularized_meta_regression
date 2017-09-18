#!/bin/csh

#$ -M claurin@nd.edu	# Email address for job notification
#$ -m abe
#$ -t 1:60
#$ -N diss_sub# Specify job name
#$ -w e             	# Ask SGE to notify you, if you make a bad resource request
#$ -q *@@daccss


module load R
cd ~/dis_sim

mkdir ~/dis_sim/dsim_${SGE_TASK_ID}

cd ~/dis_sim/dsim_${SGE_TASK_ID}

	foreach subtask (`seq 25`)
		
		R CMD BATCH --no-save --no-restore "--args ${SGE_TASK_ID}_${subtask}" ../diss_submit_v3.R diss_submit_${SGE_TASK_ID}_${subtask}.txt

	end 
	cat `ls -1 ds*rep*.txt` > ../ds_t${SGE_TASK_ID}_s${subtask}.txt
cd ..
