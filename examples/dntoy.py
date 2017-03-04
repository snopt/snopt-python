"""
This is a toy example problem.

    min (x1 + x2 + x3)**2 + 3x3 + 5x4
    st.
         x3 >= 0, x4 >= 0
         x1**2 + x2**2 + x3      = 2
                 x2**4      + x4 = 4
         2x1   + 4x2            >= 0


 The linear part of the objective 3x3 + 5x4 is given
 as part of the constraints:

    min (x1 + x2 + x3)**2 + e_4
    st.
           x3 >= 0, x4 >= 0
           x1**2 + x2**2 + x3        =  2
                   x2**4      +  x4  =  4
           2x1   + 4x2               >= 0
     inf >                3x3 + 5x4  > -inf

 The Jacobian matrix is
    [2*x1   2*x2    1.0    0]
    [0      4*x3      0  1.0]
    [2.0     4.0      0    0]
    [0.0     0.0    3.0  5.0]
 with
    iObj = 4 (indicating the linear objective term)
"""

import numpy          as np
import scipy.sparse   as sp
from   optimize import dnopt, DNOPT_options


def toycon(mode,x,fCon,jCon,nState):
    # Nonlinear terms of the gradient only
    if mode == 0 or mode == 2:
        fCon[0] = x[0]**2 + x[1]**2
        fCon[1] = x[1]**4

    if mode >= 1:
        jCon[0][0] = 2.0*x[0]
        jCon[0][1] = 2.0*x[1]
        jCon[1][1] = 4.0*x[1]**3

    return mode, fCon, jCon


def toyobj(mode,x,fObj,gObj,nState):
    sumx = x[0] + x[1] + x[2]

    # Nonlinear objective term only
    if mode == 0 or mode == 2:
        fObj = sumx**2

    if mode == 1 or mode == 2:
        gObj[0] = 2.0*sumx
        gObj[1] = 2.0*sumx
        gObj[2] = 2.0*sumx

    return mode, fObj, gObj


options = DNOPT_options()
inf      = 1.0e+20

options.setOption('Infinite bound',inf)
options.setOption('Verify level',3)
options.setOption('Print filename','dntoy.out')

mLCon = 2
mNCon = 2
n     = 4
m     = mLCon + mNCon

nnJac = 2
nnObj = 3

x0    = np.zeros(n+m,float)


# H need not be initialized for Cold starts.
# For a Warm start, H provides the initial approximation of the
#  Hessian of the Lagrangian (usually the output H from a previous run)
H = np.zeros((n,n))

# A contains the Jacobian of the linear constraints.
# J contains the sparsity pattern of the nonlinear Jacobian matrix.
#  Only assign the constant elements.  Others will be assigned in toycon/funcon.

A = np.array([ [2.0, 4.0,   0,   0],
               [0.0, 0.0, 3.0, 5.0] ])

J = np.array([ [0.0, 0.0, 1.0,   0],
               [0.0, 0.0, 0.0, 1.0] ])

bl    = -inf*np.ones(n+m)
bu    =  inf*np.ones(n+m)

bl[2] = 0.0
bl[3] = 0.0

bl[4] = 2.0
bu[4] = 2.0

bl[5] = 4.0
bu[5] = 4.0

bl[6] = 0.0

iObj  = 2

options.setOption('Verbose',True)
options.setOption('Derivative level',0)

result = dnopt(toyobj,toycon,nnObj=nnObj,nnJac=nnJac,
               x0=x0,H=H,A=A,J=J,name='  sntoyb',
               iObj=iObj,bl=bl,bu=bu,options=options)
