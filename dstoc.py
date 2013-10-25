# Tomer Levinboim, October 2013
import numpy as np


# An implementation of dykstra's projection algorithm onto the set of doubly-stochastic matrices.
def dykstra(X, tol):
    T = 1000
    P = np.zeros(X.shape)
    Q = np.zeros(X.shape)
    R = np.zeros(X.shape)
    X0 = X
    for t in xrange(T):
        X1 = proj_row(X0 + P)
        P = X0 + P - X1
        X2 = proj_col(X1 + Q)
        Q = X1 + Q - X2;
        X0 = proj_nonneg(X2 + R)
        R = X2 + R - X0
        if np.all(X2 > -tol):
            break
    return X0, t


def proj_row(X):
    N = X.shape[0]
    e = np.mat(np.ones((N, 1)))
    X = X - (X * e - e) * e.T / N
    return X


def proj_col(X):
    N = X.shape[0]
    e = np.mat(np.ones((N, 1)))
    X = X - e * (e.T * X - e.T) / N
    return X


def proj_nonneg(X):
    X[X < 0] = 0
    return X


def assert_projections():
    N = 5
    X = (np.random.randn(N, N))
    Zr = proj_row(X)
    Zc = proj_col(X)
    Z0 = proj_nonneg(X)
    e = np.mat(np.ones((N, 1)))
    assert np.linalg.norm(Zr * e - 1) < 1e-8
    assert np.linalg.norm(e.T * Zc - 1) < 1e-8
    assert np.all(Z0 >= 0)