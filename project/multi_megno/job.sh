#! /bin/sh
#$ -N Multi_MEGNO
#$ -cwd
#$ -l h_rt=04:00:00
#$ -l h_vmem=32G
module load anaconda
source activate rebound_python
python3 problem.py
source deactivate