# Packages

import sys, time

# Inputs

data_path = str(sys.argv[1])
bed_fn = data_path+".bed"
pheno_fn = data_path+".phen_w_header"
cov_fn = data_path+".cov_w_header"
pruned_bed_fn = str(sys.argv[2])+".bed"

# Parameters

maxiters = 10
sample_sketch_size = float(sys.argv[3])
marker_sketch_size = float(sys.argv[4])
block_size = int(sys.argv[5])

# MaSkLMM pacakge

import MaSkLMM as MaSkLMM

# Running MaSkLMM

start = time.time()
newton = MaSkLMM.run(bed_fn, pruned_bed_fn, pheno_fn, cov_fn, sample_sketch_size = sample_sketch_size, marker_sketch_size = marker_sketch_size, maxiters = maxiters, block_size = block_size)
end = time.time()
mask_timing = end - start
print("execution time (s):", mask_timing)
print(" ")

# Checking convergence for outputs

if newton.converged == True:
    pass
else:
    if newton.root < 0 or abs(newton.root) > 1000:
        sys.exit(0)
    else:
        pass

# end of script


