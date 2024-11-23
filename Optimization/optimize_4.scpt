#!/bin/sh 
### General options
### -- specify queue -- 
#BSUB -q hpc
### -- set the job Name -- 
#BSUB -J Optimization_4
### -- ask for number of cores (default: 1) -- 
#BSUB -n 12
### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"
### -- specify that we need 2GB of memory per core/slot -- 
#BSUB -R "rusage[mem=4GB]"
### -- specify that we want the job to get killed if it exceeds 3 GB per core/slot -- 
#BSUB -M 5GB
### -- set walltime limit: hh:mm -- 
#BSUB -W 72:00 
### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
#BSUB -u s194113@dtu.dk
### -- send notification at start --
#BSUB -B
### -- send notification at completion--
#BSUB -N
### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o Administrator_4.out
#BSUB -e Administrator_4.err
# -- end of LSF options --

source ~pkwi/WORK/miniconda3/bin/activate
module load ~pkwi/McStas/mcstas/3.dtu/mcstas-module
module load gcc/12.3.0-binutils-2.40
module load mpi/4.1.5-gcc-12.3.0-binutils-2.40 

python3 optimization_4.py

