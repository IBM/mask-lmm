#!/bin/bash
#SBATCH -A pdrineas
# SBATCH -n 128
# SBATCH -n 64
#SBATCH -N 1
#SBATCH -t 10-00:00:00

# set name of job
# SBATCH --job-name=MASK-LMM
# SBATCH -o mask.out
# SBATCH -e mask.err

# mail alert at (b)eginning, (e)nd and (a)bortion of execution
# SBATCH --mail-type=ALL
# send mail to the following address
# SBATCH --mail-user=mcburch@purdue.edu

# use submission environment
#SBATCH --export=ALL
cd "$SLURM_SUBMIT_DIR" || exit $?
source /depot/pdrineas/data/gwas_scripts/software_paths.conf

module load anaconda

PATH_TO_PLINK2='/depot/pdrineas/data/gwas_software/updated_plink2/plink2'

sample_sketch_size=${3}
marker_sketch_size=${4}
block_size=${5}
data=${1}
pruned=${2}
pvalthresh="5.0E-8"

echo " "
echo "parameters:"
echo "sample sketch dimension: ${sample_sketch_size}"
echo "marker sketch dimension: ${marker_sketch_size}"
echo "block size: ${block_size}"
echo "dataset: ${data}"
echo "dataset used for GRM sketch: ${pruned}"
echo "p-value threshold: ${pvalthresh}"
echo " "

python3 /scratch/negishi/mcburch/masklmm/mask_real.py ${data} ${pruned} ${sample_sketch_size} ${marker_sketch_size} $block_size

sort -g -k 5 masklmm_output > temp
mv temp masklmm_output
sed -i '1i\SNP Chr ChrPos ChiSq PValue' masklmm_output

awk -v awkvar="$pvalthresh" '$5 < awkvar {print}' masklmm_output > masklmm_output_hits

numhits=$(wc -l < masklmm_output_hits)

echo " "
echo "number of assocs.: ${numhits}"
echo " "


# # How to run with varied settings

# mask="/scratch/negishi/mcburch/masklmm/mask_real.sh"
# dir_prefix="HYP-pr-"
# data="../HYP_ukb_imp_merged_v3"
# #data="../qc/pruned_0.8/ukbimp_ps_1"
# #pruned="../HYP_ukb_imp_merged_v3"
# pruned="../qc/pruned_0.8/ukbimp_ps_1"

# samples=( 0.1 0.2 0.3 )
# markers=( 0.1 0.2 0.3 0.4 0.5 0.6 0.7 )

# for i in "${samples[@]}"
# do
#         for j in "${markers[@]}"
#         do
#                 #echo "$i $j"
#                 rm -rf ${dir_prefix}${i}-${j}/; mkdir ${dir_prefix}${i}-${j}/
#                 cd ${dir_prefix}${i}-${j}/; sbatch -n 128 -J ${dir_prefix}${i}-${j} $mask $data $pruned $i $j 50000; cd ..
#         done
# done