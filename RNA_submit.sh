#!/bin/bash

#SBATCH -p all 
#SBATCH --mail-user=chongm@email.sc.edu
#SBATCH --time=1-0
#SBATCH --job-name="RNA_Seq_Analysis"
#SBATCH --output=RNA_Seq_Cluster_lfdr/RNA_Seq_CV%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-100

module load R/3.3.0

Rscript --no-save --no-restore --verbose RNA_Seq_Count_CV.R ${SLURM_ARRAY_TASK_ID} 2>&1 



