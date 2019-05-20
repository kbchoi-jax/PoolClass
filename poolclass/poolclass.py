# -*- coding: utf-8 -*-

"""Main module"""

from __future__ import print_function
from . import utils
from . import get_data
from past.builtins import xrange
import os
import time
import glob
import numpy as np
from scipy.sparse import lil_matrix, dok_matrix, csr_matrix
from scipy.special import digamma, polygamma
from scipy.stats import chi2
from numpy.linalg import inv
from statsmodels.discrete.discrete_model import NegativeBinomial
import loompy
import pandas as pd
from subprocess import call
try:
    import cPickle as pickle
except:
    import pickle
try:
    xrange
except NameError:
    xrange = range

LOG = utils.get_logger()


def __update_shape(sufficient, tol=0.000001, max_iters=100):
    shape = 0.5 / sufficient
    for cur_iter in range(max_iters):
        shape_prev = shape.copy()
        g = np.log(shape) - sufficient - digamma(shape)
        h = 1/shape - polygamma(1, shape)
        shape = 1 / (1/shape + g/(np.power(shape, 2) * h))
        abs_err = np.abs(shape - shape_prev)
        abs_err_max = np.max(abs_err)
        abs_err_arg = np.argmax(abs_err)
        if abs_err_max < tol:
            LOG.debug("__update_shape ran total of %d iterations (Max error=%.9f)" % (cur_iter+1, abs_err_max))
            break
        LOG.debug('Iter #%04d: Error = %.4f [argmax = %d]' % (cur_iter, abs_err_max, abs_err_arg))
    return shape


def __em4tge(cntmat, scaler, percentile, tol, max_iters):
    cntmat_scaled = cntmat.copy()
    cntmat_scaled.data = cntmat_scaled.data / scaler[cntmat_scaled.indices]
    # libsz_scaled = np.squeeze(np.asarray(cntmat_scaled.sum(axis=0)))
    ridx = np.repeat(np.arange(cntmat.shape[0]), np.diff(cntmat.indptr))
    #
    # Initial mu and phi
    #
    x_nnz = np.squeeze(np.asarray((cntmat > 0).sum(axis=1)))
    mean_x_scaled = np.squeeze(np.asarray(cntmat_scaled.sum(axis=1))) / x_nnz
    mean_x_scaled_square = np.squeeze(np.asarray(cntmat_scaled.power(2).sum(axis=1))) / x_nnz
    var_x_scaled = mean_x_scaled_square - np.power(mean_x_scaled, 2)
    phi = var_x_scaled / mean_x_scaled - 1
    phi[phi < 0] = 0.001
    mu = mean_x_scaled / phi
    p = cntmat.copy()
    lamb = cntmat.copy()
    log_lamb = cntmat.copy()

    for cur_iter in range(max_iters):
        #
        # Initialize an iteration
        #
        mu_prev = mu.copy()

        #
        # E-step
        #
        p.data = phi[ridx] / (scaler[cntmat.indices]*phi[ridx] + 1)
        x_plus_mu = cntmat.copy()
        x_plus_mu.data += mu[ridx]
        lamb.data = x_plus_mu.data * p.data
        mean_lamb = np.squeeze(np.asarray(lamb.sum(axis=1))) / x_nnz
        log_lamb.data = digamma(x_plus_mu.data) + np.log(p.data)
        mean_log_lamb = np.squeeze(np.asarray(log_lamb.sum(axis=1))) / x_nnz
        suff = np.log(mean_lamb) - mean_log_lamb

        #
        # M-step
        #
        mu = __update_shape(suff, tol, max_iters)
        phi = mean_lamb / mu

        #
        # Check termination
        #
        err = np.abs(mu - mu_prev)
        err_pct = np.percentile(err, percentile)
        num_converged = sum(err < tol)
        if err_pct < tol:
            break
        LOG.warn('Iter #%04d: %6s genes converged below the tolerance level of %.1E' % (cur_iter+1, num_converged, tol))
        LOG.debug('Median error=%.6f' % err_pct)
    if cur_iter+1 == max_iters:
        LOG.warn('Reached the maximum number of iterations')
    return lamb, mu, phi, err


def __em4tge_dense(cntmat, scaler, percentile, tol, max_iters):
    #
    # Initial mu and phi
    #
    cntmat_scaled = cntmat / scaler
    mean_x_scaled = cntmat_scaled.mean(axis=1)
    var_x_scaled = cntmat_scaled.var(axis=1)
    phi = var_x_scaled / mean_x_scaled - 1
    phi[phi < 0] = 0.001
    mu = mean_x_scaled / phi
    phi = np.squeeze(np.asarray(phi))
    mu = np.squeeze(np.asarray(mu))

    for cur_iter in range(max_iters):
        #
        # Initialize an iteration
        #
        mu_prev = mu.copy()

        #
        # E-step
        #
        p = (1 / (np.outer(scaler, phi) + 1) * phi).T
        x_plus_mu = np.add(cntmat, mu[:, np.newaxis])
        lamb = np.multiply(x_plus_mu, p)
        mean_lamb = lamb.mean(axis=1)
        log_lamb = digamma(x_plus_mu) + np.log(p)
        suff = np.log(mean_lamb) - log_lamb.mean(axis=1)

        #
        # M-step
        #
        mu = __update_shape(np.squeeze(np.asarray(suff)), tol, max_iters)
        phi = np.divide(np.squeeze(np.asarray(mean_lamb)), mu)

        #
        # Check termination
        #
        err = np.abs(mu - mu_prev)
        err_pct = np.percentile(err, percentile)
        num_converged = sum(err < tol)
        if err_pct < tol:
            break
        LOG.warn('Iter #%04d: %6s genes converged below the tolerance level of %.1E' % (cur_iter+1, num_converged, tol))
        LOG.debug('Median error=%.6f' % err_pct)
    if cur_iter+1 == max_iters:
        LOG.warn('Reached the maximum number of iterations')
    return lamb, mu, phi, err


def run_em(loomfile, model, common_scale, percentile, tol, max_iters):
    if model == 'normalizing' and model == 'pooling':
        raise RuntimeError('The model should be either \"normalizing\" or \"pooling\".')
    with loompy.connect(loomfile) as ds:
        num_genes, num_cells = ds.shape
        LOG.info('Loading data from %s' % loomfile)
        if model == 'normalizing':
            origmat = ds.sparse().tocsr()
        elif model == 'pooling':
            origmat = ds.layers[''][:, :]
        else:
            raise NotImplementedError('Only Gamma-Poisson model is available for TGE in run_em.')
        LOG.info('Processing data matrix')
        if 'Selected' in ds.ca.keys():  # Get selected cells
            csurv = np.where(ds.ca.Selected > 0)[0]
            cntmat = origmat[:, csurv]
        else:
            csurv = np.arange(num_cells)
            cntmat = origmat
        LOG.info('The number of selected cells: %d' % len(csurv))
        
        libsz = np.squeeze(np.asarray(cntmat.sum(axis=0)))
        scaler = libsz / common_scale
        if 'Selected' in ds.ra.keys():  # Get selected genes
            gsurv1 = ds.ra.Selected > 0
        else:
            gsurv1 = np.ones(num_genes)
        gsurv2 = np.squeeze(np.asarray((cntmat > 0).sum(axis=1) > 0))
        gsurv = np.where(np.logical_and(gsurv1, gsurv2))[0]
        LOG.info('The number of selected genes: %d' % len(gsurv))
        cntmat = cntmat[gsurv, :]
        LOG.info('Running EM algorithm for TGE')
        if model == 'normalizing':
            lambda_mat, mu, phi, err = __em4tge(cntmat, scaler, percentile, tol, max_iters)
        elif model == 'pooling':
            lambda_mat, mu, phi, err = __em4tge_dense(cntmat, scaler, percentile, tol, max_iters)
        else:
            raise NotImplementedError('Only Gamma-Poisson model is available for TGE in run_em.')
        LOG.info('There were %d genes that converged below the tolerance level of %.1E' % (sum(err < tol), tol))
        LOG.info('Saving results to %s' % loomfile)
        if model == 'normalizing':
            resmat = csr_matrix((origmat.shape))
            resmat.indptr = np.ones(resmat.indptr.shape, dtype='int') * lambda_mat.indptr[-1]
            resmat.indptr[0] = 0
            resmat_indptr_part = np.repeat(lambda_mat.indptr[1:-1], np.diff(gsurv))
            resmat.indptr[1:len(resmat_indptr_part)+1] = resmat_indptr_part
            resmat.indices = csurv[lambda_mat.indices]
            resmat.data = lambda_mat.data
        elif model == 'pooling':
            resmat = np.zeros(origmat.shape)
            resmat[gsurv, :][:, csurv] = lambda_mat
        else:
            raise NotImplementedError('Only Gamma-Poisson model is available for TGE in run_em.')
        ds.layers['lambda'] = resmat
        mu_res = dok_matrix((num_genes, 1), float)
        mu_res[gsurv] = mu[:, np.newaxis]
        ds.ra['mu'] = mu_res
        phi_res = dok_matrix((num_genes, 1), float)
        phi_res[gsurv] = phi[:, np.newaxis]
        ds.ra['phi'] = phi_res
        err_res = dok_matrix((num_genes, 1), float)
        err_res[gsurv] = err[:, np.newaxis]
        ds.ra['err'] = err_res
        g_selected = dok_matrix((num_genes, 1), float)
        g_selected[gsurv] = 1
    LOG.info("Finished EM for TGE")


def __ez_test(y, X, G=None, offset=None, exposure=None):
    if G is None:
        G = X
    n = len(y)
    ind = np.squeeze((y == 0)).astype(int)
    df = max(X.shape[1], G.shape[1])
    nb = NegativeBinomial(y, X, offset=offset, exposure=exposure).fit()
    k_hat = nb.params[-1]
    B_hat = nb.params[:-1]
    tau_h = np.exp(X.dot(B_hat))
    N = 1 + k_hat*tau_h
    L = 1/k_hat
    bb1 = tau_h / N
    I_bb = np.matmul(np.matmul(X.T, np.diag(bb1)), X)
    D = N**2 * (N**L - 1)
    aa = tau_h**2 / D
    I_aa = np.matmul(np.matmul(G.T, np.diag(aa)), G)
    I_ab = np.matmul(np.matmul(G.T, -np.diag(aa)), X)
    I_ak = np.matmul(G.T, tau_h * (k_hat**(-2) * N * np.log(N) - tau_h/k_hat) / D)
    I_bk = np.zeros(X.shape[1])
    K = k_hat**(-4) + 2*k_hat**(-3)
    I_kk = sum(2*np.log(N)/k_hat**3 +
           tau_h * (1 - 2/N) / k_hat**2 -
           tau_h**2 * (tau_h + L/N) / N**2 +
           (1 - N**(-L)) * K * polygamma(1, L+1) -
           K / n * sum((1-ind)*polygamma(1, y+L+1)))
    I_abk = np.hstack((I_ab, I_ak[:, np.newaxis]))
    I_bbkk = np.vstack((np.hstack((I_bb, I_bk[:, np.newaxis])), np.hstack((I_bk, I_kk))))
    try:
        HH = I_aa - np.matmul(np.matmul(I_abk, inv(I_bbkk)), I_abk.T)
    except np.linalg.LinAlgError:
        LOG.warn("Singular matrix error: %d/%d nonzero data points" % ((y>0).sum(), n))
        return (-1, -1)
    gg = ind * tau_h / N - (1-ind) * tau_h / (N**(L+1) - N)
    U = G.T.dot(gg)
    try:
        test_S = U.dot(inv(HH)).dot(U)
        pvalue = 1 - chi2.cdf(test_S, df)
        return test_S, pvalue
    except np.linalg.LinAlgError:
        LOG.warn("Singular matrix error: %d/%d nonzero data points" % ((y>0).sum(), n))
        return (-1, -1)


def extra_zero_test(npzfile, common_scale, outfile=None):
    if outfile is None:
        outfile = os.path.join(os.path.dirname(npzfile), 'poolclass.score_test.tsv')
    LOG.warn('Name of output file: %s' % outfile)
    fhout = open(outfile, 'w')
    fhout.write('#TargetID\tNum_Zeros\tNum_Samples\tChi_sq\tP-value\n')
    fh = np.load(npzfile)
    gsurv = np.where(fh['Selected'])[0]
    csurv = np.where(fh['CellType'] != 'NA')[0]
    ctype = fh['CellType'][csurv]
    dmat = fh['Counts'][:, csurv]
    dmat = dmat[gsurv, :]
    num_genes, num_cells = dmat.shape
    LOG.warn('The number of genes = %d' % num_genes)
    LOG.warn('The number of cells = %d' % num_cells)
    if common_scale > 0:
        LOG.warn('Common scale: %d' % common_scale)
        expo = fh['Size'][csurv] / common_scale
    else:
        LOG.warn('Offset/Exposure will not be used')
        expo = None
    #X = np.ones((num_cells, 1))
    X = pd.get_dummies(ctype)
    for g in range(num_genes):
        y = dmat[g][:, np.newaxis]
        nnz = (y>0).sum()
        if nnz > 0:
            chi2val, pvalue = __ez_test(y, X, exposure=expo)
            if chi2val != -1 and not np.isnan(chi2val) and pvalue != -1 and not np.isnan(pvalue):
                fhout.write('%s\t%d\t%d\t%.6f\t%.5g\n' % (fh['GeneID'][g], num_cells-nnz, num_cells, chi2val, pvalue))
    fhout.close()


def submit(loomfile, chunk, common_scale, outdir, email, queue, mem, walltime, systype, dryrun):
    LOG.warn('Count file: %s' % loomfile)
    LOG.warn('HPC system type: %s' % systype)
    if dryrun:
        LOG.warn('Showing submission script only')
    with loompy.connect(loomfile, 'r') as ds:
        gsurv = np.where(ds.ra.Selected)[0]
        num_gsurv = len(gsurv)
        num_genes, num_cells = ds.shape
        LOG.warn('The number of selected genes: %d' % num_gsurv)
        LOG.warn('The number of selected cells: %d' % num_cells)
        LOG.warn('%d jobs will be submitted' % int(np.ceil(num_gsurv/chunk)))
    processed = 0
    if systype == 'pbs':
        tot_layer = ''
        for idx_start in xrange(0, num_gsurv, chunk):
            idx_end = min(idx_start+chunk, num_gsurv-1)
            start = gsurv[idx_start]
            if idx_end < num_gsurv-1:
                end = gsurv[idx_end]
                genes = gsurv[idx_start:idx_end]
            else:  #idx_end == num_gsurv-1:
                end = num_genes
                genes = gsurv[idx_start:]
            LOG.info('Chunk start: %d, end %d' % (start, end))
            infile = os.path.join(outdir, '_chunk.%05d-%05d.npz' % (start, end))
            LOG.debug('Genes: %s' % ' '.join(genes.astype(str)))
            LOG.debug('Total %d genes submitted in this job' % len(genes))
            data_dict = dict()
            data_dict['shape'] = (len(genes), num_cells)
            with loompy.connect(loomfile, 'r') as ds:
                data_dict['GeneID'] = ds.ra.GeneID[genes]
                data_dict['Counts'] = ds.layers[tot_layer][genes, :]
                data_dict['Size'] = ds.ca.Size
                data_dict['CellType'] = ds.ca.CellType
                data_dict['Selected'] = np.ones(len(genes))  # select all
                np.savez_compressed(infile, **data_dict)
            outfile = os.path.join(outdir, 'poolclass.score_test.%05d-%05d.tsv' % (start, end))
            job_par = 'SCALE=%d,OUTFILE=%s,INFILE=%s' % (common_scale, outfile, infile)
            cmd = ['qsub']
            if email is not None:
                cmd += ['-M', email]
            if queue is not None:
                cmd += ['-q', queue]
            if mem > 0:
                cmd += ['-l', 'mem=%d' % mem]
            if walltime > 0:
                cmd += ['-l', 'walltime=%d:00:00' % walltime]
            cmd += ['-v', job_par]
            cmd += [os.path.join(os.path.dirname(os.environ['_']), 'run_score_test_on_cluster.sh')]
            if dryrun:
                print(" ".join(cmd))
            else:
                LOG.info(" ".join(cmd))
                call(cmd)
                time.sleep(1.0)
            processed += len(genes)
        LOG.debug('Total %d genes were submitted' % processed)
        LOG.warn('Job submission complete')
    elif 'lsf':
        raise NotImplementedError('LSF submission is not yet supported')
    else:
        raise RuntimeError('No plan to support other job scheduling system until we see many requests')
