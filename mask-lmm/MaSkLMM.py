import pandas as pd
import numpy as np
np.seterr(invalid='ignore')

from os import environ
N_THREADS = 8
environ['OMP_NUM_THREADS'] = str(N_THREADS)
environ['OPENBLAS_NUM_THREADS'] = str(N_THREADS)
environ['MKL_NUM_THREADS'] = str(N_THREADS)

import sys, math
from scipy import optimize, linalg
import scipy.stats as stats
from pysnptools.snpreader import Bed, Pheno
from pysnptools.util import intersect_apply

def get_data(bed_fn, pruned_bed_fn, pheno_fn, cov_fn):
    """
    Function reading in data and performing normalization

    :param bed_fn: Filename of of PLINK Bed-formatted files (ex: "/path/to/files/toy.bed")

    :param pheno_fn: Filename of of PLINK phenotype-formatted files (ex: "/path/to/files/toy.phen")

    :param cov_fn: Filename of of PLINK covariate-formatted files (ex: "/path/to/files/toy.cov")

    """

    ignore_in = None
    #store snp data from Bedfile
    snp_on_disk = Bed(bed_fn, count_A1=False)
    prune_on_disk = Bed(pruned_bed_fn, count_A1=False)
    #store covariates
    cov_on_disk = Pheno(cov_fn)
    #store phenotype
    pheno_on_disk = Pheno(pheno_fn)

    #intersect and sort IDs from bed, pheno, cov returning new data with same IDs
    ignore_out, bed_out, prune_out, pheno_out, cov_out = intersect_apply([ignore_in,
                                                   snp_on_disk,
                                                   prune_on_disk,
                                                   pheno_on_disk,
                                                   cov_on_disk])

    sample_id = bed_out.iid
    num_samples = bed_out.iid_count
    snp_id = bed_out.sid
    snp_chrom = bed_out.pos
    snp_chrom[0][1] = 0; snp_chrom[0][2] = 0; snp_ids_and_chrom = pd.DataFrame(np.hstack([snp_id.reshape(snp_id.shape[0], 1),snp_chrom]))
    snp_ids_and_chrom.columns = ['SNP','Chr','ChrPos', 'ChrPos2']
    num_snps = bed_out.sid_count
    num_covars = cov_out.sid_count

    #read pheno, cov values
    pheno = pheno_out.read(dtype = np.float32).standardize(num_threads = N_THREADS).val
    cov = cov_out.read(dtype = np.float32).standardize(num_threads = N_THREADS).val

    return bed_out, prune_out, pheno, cov, snp_ids_and_chrom, num_samples, num_snps, num_covars

def get_K(Z):
    """
    Function computing sketched relatedness matrix

    :param Z: Sketched genotype matrix

    """

    n, m = Z.shape
    
    return np.float32((1/m)*(Z @ Z.T))

def sample_sketch(M, sk_sz, type="clkwdf"):
    """
    Function performing sketching on the samples for the genotype and phenotype matrices

    :param M: Matrix / vector to be sketched (genotype and phenotype matrix in separate calls)

    :param sk_sz: Sketch dimension to be used on the samples

    :param type: Flag for different types of sketching. Currently, support count sketching and gaussian projections

    """

    if type != "clkwdf":
        # sample sketching for Z, y, X
        n, _ = M.shape
        sketch_rows = int(sk_sz*n)
        if sk_sz < 1.0:
            S = np.random.randn(sketch_rows, n) / np.sqrt(sketch_rows)
            M_sketched = S @ M
        else:
            return M
        
        return M_sketched
    
    else:
        # sample sketching for Z, y, X
        n, _ = M.shape
        sketch_rows = int(sk_sz*n)
        if sk_sz < 1.0:
            M_sketched = np.float32(linalg.clarkson_woodruff_transform(M, sketch_rows, seed = 13))
        else:
            return M
        
        return M_sketched

def marker_sketch(M, sk_sz, type="clkwdf"):
    """
    Function performing sketching on the markers for the genotype and phenotype matrices

    :param M: Matrix / vector to be sketched (genotype and phenotype matrix in separate calls)

    :param sk_sz: Sketch dimension to be used on the markers

    :param type: Flag for different types of sketching. Currently, support count sketching and gaussian projections

    """

    if type != "clkwdf":
        # sample sketching for Z, y, X
        _, m = M.shape
        sketch_cols = int(sk_sz*m)
        if sk_sz < 1.0:
            S = np.random.randn(m, sketch_cols) / np.sqrt(sketch_cols)
            M_sketched = M @ S
        else:
            return M  
              
        return M_sketched
    
    else:
        # marker sketching for GRM
        _, m = M.shape
        sketch_cols = int(sk_sz*m)
        if sk_sz < 1.0:
            M_sketched = np.float32(linalg.clarkson_woodruff_transform(M.T, sketch_cols, seed = 13).T)
        else:
            return M
        
        return M_sketched

def get_H_tau(tau_0, n, K):
    """
    Function computing H_tau (see paper for details)

    :param tau_0: Current value for tau estimate

    :param n: Number of samples

    :param K: Sketched GRM matrix

    """

    return np.float32(K + (tau_0*np.identity(n)))

def get_U_term(X):
    """
    Function computing U_term (see paper for details)

    :param X: Matrix of covariates computed after sketching

    """

    U_x, S, _ = np.linalg.svd(X, full_matrices=False)

    return np.float32(np.identity(X.shape[0]) - (U_x @ U_x.T))

def get_projected_matrix(U_term, H_tau):
    """
    Function computing projection term 

    :param U_term: see 'get_U_term()' for details
    
    :param H_tau: see 'get_H_tau()' for details

    """

    UHU = U_term @ H_tau @ U_term

    return np.float32(UHU)

def estimate_variance_comp(y, M_inv, n, c):
    """
    Function computing sigma g squared term (see paper for details)

    :param y: Sketched phenotype vector / matrix

    :param M_inv: Pseudoinverse term (see paper for details)

    :param n: Number of samples

    :param c: Number of covariates

    """

    num = (y.T @ M_inv @ y)[0][0]

    den = n - c

    return num/den

def get_p_inverse(P):
    """
    Function computing inverse of projection matrix

    :param P: Projection matrix

    """

    # Get the indices for the off-diagonal elements
    row_indices, col_indices = np.indices(P.shape)

    # Set all off-diagonal elements to 0
    P[row_indices != col_indices] = 0

    #P_inv = np.float32(np.linalg.pinv(P, hermitian = True))
    P_inv = np.copy(P)

    diagonal_elements = np.diag_indices_from(P_inv)

    P_inv[diagonal_elements] = 1 / P_inv[diagonal_elements]

    return P_inv

def get_lle(x, n, c, pheno, K, U_term):
    """
    Function computing log-likelihood function at estimated values (see paper for details)

    :param x: Current value for tau estimate

    :param n: Number of samples

    :param c: Number of covariates

    :param pheno: Sketched phenotype vector / matrix

    :param K: Sketched GRM matrix

    :param U_term: see 'get_U_term()' for details

    """
    
    H_tau = get_H_tau(x, n, K)
    
    P = get_projected_matrix(U_term, H_tau)

    P_inv = get_p_inverse(P)
    
    sigma_g = abs(estimate_variance_comp(pheno, P_inv, n, c))

    fact = (n - c)/2

    comp1 = -1*fact*np.log(2*math.pi)
    
    with np.errstate(divide='raise'):
        try:
            print("current value of sigma_g:", sigma_g)
            print(" ")
            comp2 = -1*fact*np.log(sigma_g)
        except FloatingPointError:
            #print("current value of sigma_g:", sigma_g)
            print("ERROR: FloatingPointError occurred within lle computation")
            sys.exit(0)

    # comp2 = -1*fact*np.log(sigma_g)
    sign_logdet, logdet = np.linalg.slogdet(P)

    # subtract logdet(P) or add logdet(P^-1) to compute lle correctly
    comp3 = (1/2)*logdet # equivalent to: -(1/2)*sign_logdet*logdet

    comp4 = -1*fact

    result = np.sum([comp1, comp2, comp3, comp4])

    return result

def get_pvals(snp_on_disk, Y, X, H_tau, sigma_e, sigma_g, block_size, sample_sketch_size, snp_ids_and_chrom):
    """
    Function computing test statistics

    :param snp_on_disk: SNP data

    :param Y: Sketched phenotype vector / matrix

    :param X: Matrix of covariates computed after sketching

    :param H_tau: see 'get_H_tau()' for details

    :param sigma_e: Estimate for sigma e squared

    :param sigma_g: Estimate for sigma g squared

    :param block_size: Size of the block

    :param sample_sketch_size: Sketch dimension to use for the number of samples in the input (given 'n' samples, the sketch will have 'n * sample_sketch_size' rows)

    :param marker_sketch_size: Sketch dimension to use for the number of markers when computing the GRM (given 'm' markers, the sketch will have 'm * marker_sketch_size' columns)

    :param snp_ids_and_chrom: Array of labels for rsIDs and corresponding chromosomes

    """
    
    # normalize phenotype
    Y = np.float32((Y - np.mean(Y)).flatten())
    _, c = X.shape

    # compute V
    V = sigma_g*H_tau

    # compute inverse of V
    # V_inv = np.float32(np.linalg.pinv(sigma_g*H_tau, hermitian = True))
    V_inv = get_p_inverse(V)
    
    # vectorized computation in blocks (need to place data in sample-sketched space to compute test statistic)
    num_blocks = math.floor(snp_on_disk.sid_count / block_size)
    if snp_on_disk.sid_count % block_size != 0:
        num_blocks += 1
        
    num_sample_blocks = math.floor(snp_on_disk.iid_count / block_size)
    if snp_on_disk.iid_count % block_size != 0:
        num_sample_blocks += 1

    i = 0
    starting_range = 0
    ending_range = starting_range + block_size - 1
    # generate p-values for each block of markers
    while i < num_blocks:
        j = 0
        start_idx_samples = 0
        end_idx_samples = start_idx_samples + block_size - 1
        Z = np.empty(shape=[0, ending_range - starting_range + 1])
        while j < num_sample_blocks:
            subset_on_disk = snp_on_disk[start_idx_samples:end_idx_samples+1:1 ,starting_range:ending_range+1:1]
            block = subset_on_disk.read(dtype = np.float32).standardize(num_threads = N_THREADS).val
            M_sketched = sample_sketch(block, sample_sketch_size)
            Z = np.concatenate((Z,M_sketched),axis=0)
            j += 1
            if j+1 != num_sample_blocks:
                start_idx_samples = end_idx_samples + 1
                end_idx_samples = start_idx_samples + block_size - 1
            else:
                start_idx_samples = end_idx_samples + 1
                end_idx_samples = snp_on_disk.iid_count - 1
        i += 1
        del block
        del M_sketched
            
        n, m = Z.shape
        num_vec = np.square(Z.T @ V_inv @ Y, dtype="float32")
        den_vec = np.float32(np.einsum('...i,...i->...', Z.T @ V_inv, Z.T))

        # compute chi-square statistics 
        chi2stats = np.absolute( np.array(num_vec / den_vec).reshape((m, 1)) )
        chi2stats = np.nan_to_num(chi2stats)
        adjchisq = 2 * np.array(chi2stats)
        adjpvals = np.float64(stats.f.sf(adjchisq, dfn = 1, dfd = n-(c+1))[:,0])

        # save p-values
        snp_ids_and_chrom_block = snp_ids_and_chrom.iloc[starting_range:ending_range + 1,:][['SNP','Chr','ChrPos2']]
        snp_ids_and_chrom_block['Chisq'] = np.array(adjchisq).reshape(adjchisq.shape[0], 1)
        snp_ids_and_chrom_block['PValue'] = np.array(adjpvals).reshape(adjpvals.shape[0], 1)
        snp_ids_and_chrom_block.to_csv('masklmm_output', mode='a', sep = '\t', index=False, header=False)
        
        if i+1 != num_blocks:
            starting_range = ending_range + 1
            ending_range = starting_range + block_size - 1
        else:
            starting_range = ending_range + 1
            ending_range = snp_on_disk.sid_count - 1
        
    del num_vec 
    del den_vec 

def generateSketch(block_num, block_size, marker_sketch_size, sample_sketch_size, snp_on_disk):
    """
    Function generating sketch in blocks.

    :param block_num: Block number

    :param block_size: Size of the block

    :param sample_sketch_size: Sketch dimension to use for the number of samples in the input (given 'n' samples, the sketch will have 'n * sample_sketch_size' rows)

    :param marker_sketch_size: Sketch dimension to use for the number of markers when computing the GRM (given 'm' markers, the sketch will have 'm * marker_sketch_size' columns)

    :param snp_on_disk: SNP data

    """
    
    # number of blocks for the markers given the block size
    num_marker_blocks = math.floor(snp_on_disk.sid_count / block_size)
    if snp_on_disk.sid_count % block_size != 0:
        num_marker_blocks += 1
    
    # number of blocks for the samples given the block size
    num_sample_blocks = math.floor(snp_on_disk.iid_count / block_size)
    if snp_on_disk.iid_count % block_size != 0:
        num_sample_blocks += 1
        
    # start / end range for marker block (using block id)
    start_idx_markers = (block_num) * block_size
    end_idx_markers = start_idx_markers + block_size - 1
    if block_num+1 == num_marker_blocks:
        end_idx_markers = snp_on_disk.sid_count - 1

    j = 0
    # start / end range for sample blocks
    start_idx_samples = 0
    end_idx_samples = start_idx_samples + block_size - 1
    sketchedBlock = np.empty(shape=[0, end_idx_markers - start_idx_markers + 1])
    # sketch samples in each block and concatenate
    while j < num_sample_blocks:
        subset_on_disk = snp_on_disk[start_idx_samples:end_idx_samples+1:1 ,start_idx_markers:end_idx_markers+1:1]
        block = subset_on_disk.read(dtype = np.float32).standardize(num_threads = N_THREADS).val
        M_sketched = sample_sketch(block, sample_sketch_size)
        sketchedBlock = np.concatenate((sketchedBlock,M_sketched),axis=0)
        j += 1
        if j+1 != num_sample_blocks:
            start_idx_samples = end_idx_samples + 1
            end_idx_samples = start_idx_samples + block_size - 1
        else:
            start_idx_samples = end_idx_samples + 1
            end_idx_samples = snp_on_disk.iid_count - 1
    
    # sketch markers in concatenated block
    M_sketched = marker_sketch(sketchedBlock, marker_sketch_size)
    
    start_index = (block_num) * (block_size * marker_sketch_size)
    end_index = start_index + (block_size * marker_sketch_size)
    if block_num+1 == num_marker_blocks:
        end_index = snp_on_disk.sid_count * marker_sketch_size
    
    # place sketched block in appropriate place of array
    normZ_s1_s2[:, int(start_index):int(end_index)] = M_sketched 
    
def MaSkLMM(snp_on_disk, prune_on_disk, pheno, cov, num_covars, sample_sketch_size, marker_sketch_size, snp_ids_and_chrom, maxiters, block_size):
    """
    Function performing matrix sketching-based linear mixed modeling for association studies.

    :param snp_on_disk: SNP data

    :param prune_on_disk: Pruned SNP data to be used in GRM approximation (can be original data)

    :param pheno: Normalized phenotype vector / matrix

    :param cov: Set to 'None' as top 2 principal components computed after sketching

    :param num_covars: Number of covariates

    :param sample_sketch_size: Sketch dimension to use for the number of samples in the input (given 'n' samples, the sketch will have 'n * sample_sketch_size' rows)

    :param marker_sketch_size: Sketch dimension to use for the number of markers when computing the GRM (given 'm' markers, the sketch will have 'm * marker_sketch_size' columns)

    :param snp_ids_and_chrom: Array of labels for rsIDs and corresponding chromosomes

    :param maxiters: Maximum number of iterations to be used in the Newton-Raphson estimation

    :param block_size: Size of the block

    """
    
    # identify number of column blocks using user-defined block size 
    num_marker_blocks = math.floor(prune_on_disk.sid_count / block_size)
    if prune_on_disk.sid_count % block_size != 0:
        num_marker_blocks += 1

    # need to define empty sketch as global so that all threads / processes have access
    global normZ_s1_s2 
    normZ_s1_s2 = np.empty(shape=[int(sample_sketch_size*prune_on_disk.iid_count), int(marker_sketch_size*prune_on_disk.sid_count)], dtype = np.float32)
    num_samples, _ = normZ_s1_s2.shape
    
    # list of block ids
    blocklist = list(range(num_marker_blocks))
    
    # generate sketched input in blocks
    for block in blocklist:
        generateSketch(block, block_size, marker_sketch_size, sample_sketch_size, prune_on_disk)

    # sketch phenotype and covariates              
    pheno = sample_sketch(pheno, sample_sketch_size)
    cov = sample_sketch(cov, sample_sketch_size)
    
    # compute GRM term ("constants")
    K = get_K(normZ_s1_s2)
    del normZ_s1_s2
    
    # compute projection term ("constants")
    U_term = get_U_term(cov)
        
    # scipy newton raphson (using secant to auto-compute the derivate)
    root, newton_output = optimize.newton(get_lle, 1.0, args=( num_samples,
                                                               num_covars,
                                                               pheno,
                                                               K,
                                                               U_term, ), rtol=1e-2, full_output = True, 
                                                               maxiter=maxiters, disp = False)
    
    # test for convergence and stability of the root (as to not waste time computing test statistic)
    print(newton_output)
    print(' ')
    if newton_output.converged == True:
        pass
    else:
        if abs(newton_output.root) > 5000:
            print("METHOD DID NOT CONVERGE; aborting...")
            return newton_output
        else:
            # if root hasn't diverged too much then attempt at solution
            pass
    
    # generate terms for test statistics (additional computations after converging on sigma_g)
    H_tau = get_H_tau(abs(root), num_samples, K)
    del K
    P = get_projected_matrix(U_term, H_tau)
    P_inv = get_p_inverse(P)
    del P 
    sigma_g = estimate_variance_comp(pheno, P_inv, num_samples, num_covars)
    del P_inv 

    # compute p-values 
    get_pvals(snp_on_disk, pheno, cov, H_tau, abs(root)*sigma_g, sigma_g, block_size, sample_sketch_size, snp_ids_and_chrom)
        
    return newton_output

def run(bed_fn, pruned_bed_fn, pheno_fn, cov_fn = None, sample_sketch_size = 0.5, marker_sketch_size = 0.5, maxiters = 1, block_size = 10000):
    """
    Function performing matrix sketching-based linear mixed modeling for association studies.

    :param bed_fn: Filename of of PLINK Bed-formatted files (ex: "/path/to/files/toy.bed")

    :param pruned_bed_fn: Filename of of PLINK Bed-formatted files (ex: "/path/to/files/toy.bed")

    :param pheno_fn: Filename of of PLINK phenotype-formatted files (ex: "/path/to/files/toy.phen")

    :param cov_fn: Filename of of PLINK covariate-formatted files (ex: "/path/to/files/toy.cov")

    :param sample_sketch_size: Sketch dimension to use for the number of samples in the input (given 'n' samples, the sketch will have 'n * sample_sketch_size' rows)

    :param marker_sketch_size: Sketch dimension to use for the number of markers when computing the GRM (given 'm' markers, the sketch will have 'm * marker_sketch_size' columns)

    :param maxiters: Maximum number of iterations to be used in the Newton-Raphson estimation

    :param block_size: Size of the block

    """

    # read in data
    snp_on_disk, prune_on_disk, pheno, cov, snp_ids_and_chrom, num_samples, num_snps, num_covars = get_data(bed_fn, pruned_bed_fn, pheno_fn, cov_fn)

    # run MaSkLMM
    newton = MaSkLMM(snp_on_disk, prune_on_disk, pheno, cov, 
                     num_covars, sample_sketch_size, marker_sketch_size, snp_ids_and_chrom, maxiters, block_size)

    # check for convergence 
    if newton.converged == True:
        pass
    else:
        if newton.root < 0 or abs(newton.root) > 1000:
            #print("METHOD DID NOT CONVERGE; aborting...")
            return None,newton
            #sys.exit(0)
        else:
            # if root hasn't diverged too much then attempt at solution
            pass

    return newton
