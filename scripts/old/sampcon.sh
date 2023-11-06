#$ -S /bin/bash
#$ -N mTORC_sampcon
#$ -o /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis
#$ -j y
#$ -l h_rt=2:00:00
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -q gpu.q
#$ -l hostname='qb3-id*'

module load CBI conda-stage
conda activate imp
module load Sali python3/pyrmsd/4.1.gita558b8a
module load cuda/7

mkdir /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis/analysis.413611/sampcon30
cd /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis/analysis.413611/sampcon30

python3 /wynton/home/sali/mhancock/imp-sampcon/pyext/src/exhaust.py --sysname mtorc2 --path /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis/analysis.413611 --rmfA /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis/analysis.413611/hdbscan/A_1_10000.rmf3 --rmfB /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis/analysis.413611/hdbscan/B_1_10000.rmf3 --scoreA /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis/analysis.413611/hdbscan/A_1_10000.txt --scoreB /wynton/home/sali/mhancock/mtorc2/exp_8/66/analysis/analysis.413611/hdbscan/B_1_10000.txt --selection /wynton/home/sali/mhancock/mtorc2/scripts/custom_ranges.txt --density /wynton/home/sali/mhancock/mtorc2/scripts/custom_ranges.txt --skip --cluster_threshold 30 --align -m cuda -c 4 -g 2.0 -gp