module load python3/pyrmsd/4.1.gita558b8a
module load imp-fast/2.14.0

export top_dir=/salilab/park3/rakesh/Projects/3-SPOTSComplex/3-IntModeling/ClosedFinal100
export analys_dir=$top_dir/analys/
export name=SPOTS

python /home/rakesh/WORK/Projects/SPOTSComplex/IMPModeling/AnalysisModules/imp-sampcon/pyext/src/exhaust.py --sysname $name --path $analys_dir --mode cuda --cores 4 --align --density density.txt --gridsize 5.0 --gnuplot --scoreA A_models_clust-1.txt --scoreB B_models_clust-1.txt --rmfA A_models_clust-1.rmf3 --rmfB B_models_clust-1.rmf3
