#!/bin/sh
#SBATCH --time=100:05:00 --mem-per-cpu=4000
#SBATCH -p batch
#SBATCH -o toy_final2-%a.out
#SBATCH --array=1-20

. /etc/profile.d/modules.sh
module load matlab
matlab -r "run_MVT_allpara_emdbht($SLURM_ARRAY_TASK_ID,'V2_1_1')"
matlab -r "run_MVT_allpara_emdbht($SLURM_ARRAY_TASK_ID,'V2_1_2')"
matlab -r "run_MVT_allpara_emdbht($SLURM_ARRAY_TASK_ID,'V2_1_3')"
matlab -r "run_MVT_allpara_emdbht($SLURM_ARRAY_TASK_ID,'V2_2_2')"
matlab -r "run_MVT_allpara_emdbht($SLURM_ARRAY_TASK_ID,'V2_2_3')"
matlab -r "run_MVT_allpara_emdbht($SLURM_ARRAY_TASK_ID,'V2_3_3')"

