#$ -S /bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=48:00:00
#$ -l mem_free=10G
#$ -l scratch=10G
#$ -N MTORC_analysis
#$ -pe smp 40
#$ -l hostname='qb3-id*'

#echo "$TMPDIR"
eval "$(conda shell.bash hook)"
conda activate imp
echo $NSLOTS
python3 ../scripts/04_run_analysis_trajectories.py
#scp -r "$TMPDIR"/mtorc2/rigid_body_1/data/output/* /wynton/home/sali/mhancock/mtorc2/rigid_body_1/data/output

echo "JOB PERFORMANCE"
[[ -n "$TMPDIR" ]] && qstat -j "$JOB_ID"
