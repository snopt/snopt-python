"""
Diet problem (LP) solved with SNOPTA
"""

import numpy as np
import scipy.sparse as sp
from optimize.snopt7 import SNOPT_solver


def dieta_fun(status,x,needF,needG,cu,iu,ru):
    # LP has no nonlinear terms in the objective
    F = []
    G = []
    return status, F, G


inf   = 1.0e20
snopt = SNOPT_solver()

snopt.setOption('Print file','dieta.out')
snopt.setOption('Minor print level',1)
snopt.setOption('Summary frequency',1)
snopt.setOption('Print frequency',1)

nF      = 4
ObjRow  = 4
n       = 6

# We provide the linear components of the Jacobian
# matrix as a dense matrix.
A       = np.array([ [110, 205, 160, 160, 420, 260],
                     [  4,  32,  13,   8,   4, 14],
                     [  2,  12,  54, 285,  22, 80],
                     [  3,  24,  13,   9,  20, 19]],float)

# Alternatively, the user can provide the Jacobian in two
# other formats:

# 1. in coordinate form via scipy's coordinate matrix class
# A = sp.coo_matrix(A)

# 2. in coordinate form explicitly
# iAfun = ... row indices
# jAvar = ... col indices
# A     = ... matrix values

x0      = np.ones(n)

xlow    = np.zeros(n)
xupp    = np.array([ 4, 3, 2, 8, 2, 2],float)

Flow    = np.array([ 2000, 55, 800, -inf ],float)
Fupp    = inf*np.ones(nF)

snopt.setOption('Verbose',True)
snopt.snopta(name='   dieta',xlow=xlow,xupp=xupp,x0=x0,
             Flow=Flow,Fupp=Fupp,ObjRow=ObjRow,A=A,
             usrfun=dieta_fun)
