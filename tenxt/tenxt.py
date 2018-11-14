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


def extra_zero_test(cntfile):
    raise NotImplementedError('Coming soon!')