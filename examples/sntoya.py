"""
This is the toy example problem in the SNOPT documentation.
     min         x1
     s.t   .   x1**2 + 4x2**2 <= 4
           (x1-2)**2 +  x2**2 <= 5
            x1 >=0

We define the function F(x):
   F(x)  = [     x1             ]
           [     x1**2 + 4x2**2 ]
           [ (x1-2)**2 +  x2**2 ]
with ObjRow = 1 to indicate the objective.

The Jacobian is:
  F'(x)  = [   1       0  ]   [ 1 0 ]    [   0       0  ]
           [   2x1    8x2 ] = [ 0 0 ]  + [   2x1    8x2 ]
           [ 2(x1-2)  2x2 ]   [ 0 0 ]    [ 2(x1-2)  2x2 ]
                             linear(A)     nonlinear (G)

"""

import numpy           as np
import scipy.sparse    as sp
from   optimize.snopt7 import SNOPT_solver


def sntoya_objF(status,x,needF,needG,cu,iu,ru):
    F = np.array([                      x[1], # objective row
                   x[0]**2        + 4.0*x[1]**2,
                  (x[0] - 2.0)**2 +     x[1]**2 ])
    return status, F


def sntoya_objFG(status,x,needF,needG,cu,iu,ru):
    F = np.array([                      x[1], # objective row
                   x[0]**2        + 4.0*x[1]**2,
                  (x[0] - 2.0)**2 +     x[1]**2 ])

    G = np.array([ 2*x[0], 8*x[1], 2*(x[0]-2), 2*x[1] ])
    return status, F, G



inf   = 1.0e20

snopt = SNOPT_solver()

snopt.setOption('Verbose',True)
snopt.setOption('Solution print',True)
snopt.setOption('Print file','sntoya.out')


# Either dtype works, but the names for x and F have to be of
# the correct length, else they are both ignored by SNOPT:
xNames  = np.array([ '      x0', '      x1' ])
FNames  = np.array([ '      F0', '      F1', '      F2' ],dtype='c')

x0      = np.array([ 1.0, 1.0 ])

xlow    = np.array([ 0.0, -inf])
xupp    = np.array([ inf,  inf])

Flow    = np.array([ -inf, -inf, -inf ])
Fupp    = np.array([  inf,  4.0,  5.0 ])

ObjRow  = 1

# We first solve the problem without providing derivative info
snopt.snopta(name=' sntoyaF',x0=x0,xlow=xlow,xupp=xupp,
             Flow=Flow,Fupp=Fupp,ObjRow=ObjRow,
             usrfun=sntoya_objF,xnames=xNames,Fnames=FNames)

# Now we set up the derivative structures...

# A and G provide the linear and nonlinear components of
# the Jacobian matrix, respectively. Here, G and A are given
# as dense matrices.
#
# For the nonlinear components, enter any nonzero value to
# indicate the location of the nonlinear deriatives (in this case, 2).
# A must be properly defined with the correct derivative values.

A = np.array([ [0, 1],
               [0, 0],
               [0, 0]])

G = np.array([ [0, 0],
               [2, 2],
               [2, 2]])

# Alternatively, A and G can be input in coordinate form via scipy's
# coordinate matrix
#  A = sp.coo_matrix(A)
#  G = sp.coo_matrix(G)
# or explicitly in coordinate form
#  iAfun = row indices of A
#  jAvar = col indices of A
#  A     = matrix values of A
#
#  iGfun = row indices of G
#  jGvar = col indices of G

snopt.snopta(name='sntoyaFG',usrfun=sntoya_objFG,x0=x0,xlow=xlow,xupp=xupp,
             Flow=Flow,Fupp=Fupp,ObjRow=ObjRow,A=A,G=G,xnames=xNames,Fnames=FNames)
