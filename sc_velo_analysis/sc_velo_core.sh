#! /bin/bash

### Submit this Script with: sbatch <script.sh> ###

# Parameters for slurm (don't remove the # in front of #SBATCH!)
#  Use partition shortterm/debug/longterm:
#SBATCH --partition=shortterm
#  Use so many node:
#SBATCH --nodes=1
#  Request so many cores (hard constraint):
#SBATCH -c 4
#  Request so much of memory (hard constraint):
#SBATCH --mem=250GB
#  Find your job easier with a name:
#SBATCH --job-name=sc_velo
#set slurm file output nomenclature
#SBATCH --output "slurm-%x-%j.out"

PATH=$WORK/.omics/anaconda3/bin:$PATH #add the anaconda installation path to the bash path
source $WORK/.omics/anaconda3/etc/profile.d/conda.sh # some reason conda commands are not added by default

# Load your necessary modules:
conda activate scVelocity
module load nextflow/v22.04.1

# Move to SCRATCH all the relevant scripts.
cp * $SCRATCH/
cp -r ../src $SCRATCH/
cd $SCRATCH

# Submit the Nextflow Script:
nextflow run sc_velo.nf -params-file sc_velo_params.yaml --id ${SCRATCH/"/scratch/"/}

