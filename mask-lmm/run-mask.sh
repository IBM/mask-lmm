#!/bin/bash

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

python3 run-mask.py ${data} ${pruned} ${sample_sketch_size} ${marker_sketch_size} $block_size

sort -g -k 5 masklmm_output > temp
mv temp masklmm_output
sed -i '1i\SNP Chr ChrPos ChiSq PValue' masklmm_output

awk -v awkvar="$pvalthresh" '$5 < awkvar {print}' masklmm_output > masklmm_output_hits

numhits=$(wc -l < masklmm_output_hits)

echo " "
echo "number of assocs.: ${numhits}"
echo " "
