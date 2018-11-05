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
from scipy.sparse import lil_matrix, dok_matrix, csr_matrix, hstack
from scipy.special import digamma, polygamma
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


def __em4tgx(cntmat, scaler, percentile, tol, max_iters):
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


def __em4tgx_dense(cntmat, scaler, percentile, tol, max_iters):
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


def score_test(cntfile):
    raise NotImplementedError('Coming soon!')