"""
  Python interface for SNOPT
  http://ccom.ucsd.edu/~optimizers

  Elizabeth Wong, University of California, San Diego
  Philip Gill,    University of California, San Diego

  Feb 2015
"""

# Python modules
import os, sys
import numpy as np
import scipy.sparse as sp

import snopt7_python        as     snopt
from   optimize.sn_options  import SNOPT_options

date = '''Feb 2015'''

#-------------------------------------------------------------------------------#

class SNOPT_solver(object):
    '''
    SNOPT_solver class contains:

    Functions:
    '''
    def __init__(self,name=''):
        self.name    = name
        self.options = SNOPT_options()

        # Set up workspace
        self.lencw   = 500
        self.leniw   = 500
        self.lenrw   = 500
        self.cw      = np.array(['        ']*self.lencw,dtype='c')
        self.iw      = np.zeros(self.leniw,'i',order='F')
        self.rw      = np.zeros(self.lenrw,float,order='F')

        # Initialize SNOPT workspace
        # parameters are set to undefined values
        # title printing is disabled here
        snopt.sninit_wrap(self.cw,self.iw,self.rw)


    def getOption(self,key):
        ''' Return the value of a given option

        Keyword arguments:
        key     is the string containing the name of the option

        '''
        return self.options.getOption(key)


    def setOption(self,key,value):
        ''' Set the value of a given option

        Keyword arguments:
        key     is the string containing the name of the option
        value   is the value of the option

        '''
        self.options.setOption(key,value)


    def printOptions(self):
        ''' Print all option keywords and their values
        '''
        self.options.printOptions()


    def resetOption(self,key):
        ''' Reset a given option to its default value

        Keyword arguments:
        key     is the string containing the name of the option
        '''
        self.options.resetOption(key)


    def snopta(self,usrfun,**kwargs):
        """
        snopta() calls the SNOPTA solver to solve the optimization
        problem:
              min        F_objrow(x)
              s.t.   xlow <=   x  <= xupp
                     Flow <= F(x) <= Fupp

        where F(x) is a vector of smooth linear and nonlinear constraint functions
        and F_objrow(x) is one of the components of F to be minimized, as specified by input
        argument ObjRow.

        The user should provide a function that returns the function F(x) and ideally
        their gradients.


        Keyword arguments:
        usrfun     is the user-defined function that returns the objective, constraints
                   and gradients
        name       is the name of the problem (default '')
        n          is the number of constraints (default determined by size of x-arrays)
        nF         is the number of variables (default determined by size of F-arrays)

        A          is the linear part of the constraint Jacobian.  It can be provided as a
                   dense nF x n matrix via A or in its coordinate form via keyword arguments:
                   iAfun   containing the row indices
                   jAvar   containing the column indices
                   A       containing the elements of A such that A(j) is in location
                           (iAfun(j),jAvar(j))

        G          is the nonlinear part of the constraint Jacobian.  It can be provided as a
                   dense nF x n matrix via G or in its coordinate form via keyword arguments:
                   iGfun   containing the row indices
                   jGvar   containing the column indices
                   (Note the values of G are provided in the user-defined subroutine.)

        ObjRow     is an integer indicating the location of the objective in F(x) (default 0)
        ObjAdd     is a constant term of the objective (default 0.0)

        xnames     is an array of names for the x-variables (default '')
        xstate     is an array of names for the x-variables (default 0)
        x          are the intial values of x (default 0.0)
        xlow       are the lower bounds of x (default -infinity)
        xupp       are the upper bounds of x (default +infinity)

        Fnames     is an array of names for the F-variables (default '')
        Fstate     is an array of names for the F-variables (default 0)
        F          are the intial value of F (default 0.0)
        Flow       are the lower bounds of F (default -infinity)
        Fupp       are the upper bounds of F (default +infinity)


        """
        verbose  = self.options.getOption('Verbose')
        maxTries = self.options.getOption('Max memory attempts')

        name    = kwargs.get('name',self.name)
        nF      = kwargs.get('nF',0)
        n       = kwargs.get('n',0)

        ObjRow  = kwargs.get('ObjRow',0)
        ObjAdd  = kwargs.get('ObjAdd',0.0)

        xnamesU = kwargs.get('xnames',None)
        FnamesU = kwargs.get('Fnames',None)

        xstate  = kwargs.get('xstate',None)
        x       = kwargs.get('x0',None)
        xlow    = kwargs.get('xlow',None)
        xupp    = kwargs.get('xupp',None)

        Fstate  = kwargs.get('Fstate',None)
        F       = kwargs.get('F',None)
        Flow    = kwargs.get('Flow',None)
        Fupp    = kwargs.get('Fupp',None)

        if n is 0:
            nx      = 0 if x      is None else x.size
            nxlow   = 0 if xlow   is None else xlow.size
            nxupp   = 0 if xupp   is None else xupp.size
            nxstate = 0 if xstate is None else xstate.size
            n       = max(nx,nxlow,nxupp,nxstate)

        if nF is 0:
            nfF     = 0 if F      is None else F.size
            nFlow   = 0 if Flow   is None else Flow.size
            nFupp   = 0 if Fupp   is None else Fupp.size
            nFstate = 0 if Fstate is None else Fstate.size
            nF      = max(nfF,nFlow,nFupp,nFstate)

        if nF == 0 or n == 0:
            return 99

        usrA   = kwargs.get('A',None)
        neA    = 0
        if   type(usrA) is tuple: # A = (A,iAfun,jAvar)
            (A,iAfun,jAvar) = usrA
            neA             = A.size
        elif type(usrA) is np.ndarray and usrA.ndim == 2: # A = matrix
            A      = sp.coo_matrix(usrA)
            iAfun  = A.row
            jAvar  = A.col
            neA    = A.nnz
            A      = A.data
        elif type(usrA) is sp.coo_matrix: # A = scipy coo mtx
            iAfun  = usrA.row
            jAvar  = usrA.col
            A      = usrA.data
            neA    = usrA.nnz
        else:
            if verbose:
                print ' --> Linear component of Jacobian not provided'


        usrG   = kwargs.get('G',None)
        neG    = 0
        if   type(usrG) is tuple:  # G = (iGfun,jGvar)
            (iGfun,jGvar) = usrG
            neG           = iGfun.size
        elif type(usrG) is np.ndarray and usrG.ndim == 2: # G = matrix
            G      = sp.coo_matrix(usrG)
            iGfun  = G.row
            jGvar  = G.col
            neG    = G.nnz
        elif type(usrG) is sp.coo_matrix: # G = scipy coo mtx
            iGfun  = usrG.row
            jGvar  = usrG.col
            neG    = usrG.nnz
        else:
            if verbose:
                print ' --> Nonlinear component of Jacobian not provided'


        if xnamesU is None:
            nxname = 1
            xnames = None
        else:
            try:
                nxname = xnamesU.shape[0]
                if nxname != n:
                    nxname = 1
                    xnames = None
                else:
                    if xnamesU.dtype.char == 'c':
                        xnames  = xnamesU
                    else:
                        xnames = convertC(xnamesU,nxname)
            except:
                nxname = 1
                xnames = None

        if FnamesU is None or xnames is None:
            nFname = 1
            Fnames = None
        else:
            try:
                nFname = FnamesU.shape[0]
                if nFname != nF:
                    nFname = 1
                    Fnames = None
                else:
                    if FnamesU.dtype.char == 'c':
                        Fnames = FnamesU
                    else:
                        Fnames = convertC(FnamesU,nFname)
            except:
                nFname = 1
                Fnames = None

        if (nxname != n and nxname != 1) or (nFname != nF and nFname != 1 ):
            nxname = 1
            xnames = None
            nFname = 1
            Fnames = None


        lencu = 1
        leniu = 1
        lenru = 1

        cu    = np.array(['        ']*lencu,'c')
        iu    = np.zeros(leniu,'i',order='F')
        ru    = np.zeros(lenru,float,order='F')


        # Copy python-defined options to workspace
        #    This should also set iPrint,iSumm to
        #    defaults or user-def values.
        self.options.copyOptions(self.iw,self.rw)


        # Open files for printing
        iPrint  = self.options.getOption('Print unit')
        iSumm   = self.options.getOption('Summary unit')
        prtfile = self.options.getOption('Print file')
        sumfile = self.options.getOption('Summary file')
        snopt.snopen_wrap(iPrint,prtfile,iSumm,sumfile)

        # Print the title
        snopt.sntitle_wrap(self.iw)


        # Read specs file if one was given
        info = 0
        spcfile = self.options.getOption('Specs file')
        if spcfile is not None:
            iSpecs  = self.options.getOption('Specs unit')
            if iSpecs > 0:
                info = snopt.snspec_wrap(iSpecs,spcfile,self.cw,self.iw,self.rw)
            else:
                info = 131
            if info != 101 and info != 104:
                print 'Specs read failed', info
                return info


        # Get Jacobian structure if necessary:
        if neG == 0 and neA == 0:
            print '  Could not determine Jacobian structure from user input'
            print '  Calling snJac...'

            len   = nF*n
            count = 1
            while True:
                snjac = snopt.snjac_wrap(nF,usrfun,len,len,\
                                         x,xlow,xupp,cu,iu,ru,\
                                         self.cw,self.iw,self.rw)
                info  = snjac[0]
                iAfun = snjac[1]
                jAvar = snjac[2]
                neA   = snjac[3]
                A     = snjac[4]
                iGfun = snjac[5]
                jGvar = snjac[6]
                neG   = snjac[7]
                x     = snjac[8]
                mincw = snjac[9]
                miniw = snjac[10]
                minrw = snjac[11]

                if info == 102:
                    self.setOption('Derivative option',0)
                    self.iw[103] = 0
                    break

                if info == 82 or info == 83 or info ==84:
                    count+= 1
                    if count > maxTries:
                        print ' Could not allocate memory for SNOPT'
                        return info
                    self.resizeW(mincw,miniw,minrw)

        else:
            iAfun = np.array([1]) if neA == 0 else iAfun + 1
            jAvar = np.array([1]) if neA == 0 else jAvar + 1

            iGfun = np.array([1]) if neG == 0 else iGfun + 1
            jGvar = np.array([1]) if neG == 0 else jGvar + 1


        # Check memory
        [info,mincw,miniw,minrw] = snopt.snmema_wrap(nF,n,nxname,nFname,neA,neG,\
                                                     self.cw,self.iw,self.rw)
        self.resizeW(mincw,miniw,minrw)


        # Check problem info
        inf = self.options.getOption('Infinite bound')

        if xstate is None:
            if verbose:
                print ' --> Initial xstate not provided; setting to 0'
            xstate =  np.zeros(n,'i')
        if Fstate is None:
            if verbose:
                print ' --> Initial Fstate not provided; setting to 0'
            Fstate =  np.zeros(nF,'i')

        if x is None:
            if verbose:
                print ' --> Initial x not provided; setting to 0'
            x = np.zeros(n,float)
        if F is None:
            F = np.zeros(nF,float)

        if xlow is None:
            if verbose:
                print ' --> Lower bounds on x not provided; setting to -inf'
            xlow = -inf*np.ones(n,float)
        if xupp is None:
            if verbose:
                print ' --> Upper bounds on x not provided; setting to +inf'
            xupp =  inf*np.ones(n,float)
        if Flow is None:
            if verbose:
                print ' --> Lower bounds on F not provided; setting to -inf'
            Flow = -inf*np.ones(nF,float)
        if Fupp is None:
            if verbose:
                print ' --> Upper bounds on F not provided; setting to +inf'
            Fupp =  inf*np.ones(nF,float)

        xmul = np.zeros(n,float)
        Fmul = np.zeros(nF,float)


        # Solve problem
        Start = self.options.getOption('Start type')
        iStart = 0
        if   Start is 'Warm':
            iStart = 1
        elif Start is 'Hot':
            iStart = 2

        if verbose:
            print ''.join('-' for j in range(82))
            print ' SNOPT python interface   (' + date + ')'
            print '   Problem: ' + name
            print '   # variables = %d; # constraints = %d \n' %(n,nF)
            print '     Name state(j)      low(j)           x(j)           upp(j)        mul(j)'
            form  = ' {0:>8s}{1:8d}{2:16.6e}{3:16.6e}{4:16.6e}{5:16.6e}'

            pxnames = convertS(xnames,nxname,n)
            pFnames = convertS(Fnames,nFname,nF)
            print 'x:'
            for j in range(n):
                print form.format(pxnames[j],xstate[j],xlow[j],x[j],xupp[j],xmul[j])
            print 'F:'
            for j in range(nF):
                print form.format(pFnames[j],Fstate[j],Flow[j],F[j],Fupp[j],Fmul[j])
            print ''.join('-' for j in range(82))

        if xnames is None:
            nxname = 1
            xnames = np.array(['        '],dtype='c')
        if Fnames is None:
            nFname = 1
            Fnames = np.array(['        '],dtype='c')

        count = 1
        while True:
            res    = snopt.snopta_wrap(iStart, ObjAdd, ObjRow, name, usrfun,\
                                       iAfun, jAvar, neA, A, iGfun, jGvar, neG,\
                                       xlow, xupp, xnames, Flow, Fupp, Fnames,
                                       x, xstate, xmul, F, Fstate, Fmul,\
                                       cu, iu, ru, self.cw, self.iw, self.rw)
            if  res[6]/10 == 8:
                count += 1
                if count > maxTries:
                    print ' Could not allocate memory for SNOPT'
                    return info
                self.resizeW(res[9],res[10],res[11])
            else:
                break

        # Results
        # res[0]  = x
        # res[1]  = xstate
        # res[2]  = xmul
        # res[3]  = f
        # res[4]  = fstate
        # res[5]  = fmul
        # res[6]  = info
        # res[7]  = itn
        # res[8]  = mjritn
        # res[9]  = mincw
        # res[10] = miniw
        # res[11] = minrw
        # res[12] = nS
        # res[13] = nInf
        # res[14] = sInf
        # res[15] = Obj

        # Finish up
        snopt.close_wrap(iPrint,iSumm)


        # Return solution
        self.name      = name
        self.x         = res[0]
        self.xstate    = res[1]
        self.xmul      = res[2]
        self.F         = res[3]
        self.Fstate    = res[4]
        self.Fmul      = res[5]

        self.exit      = (res[6]/10) * 10
        self.info      = res[6]

        self.iterations= res[7]
        self.major_itns= res[8]

        self.nS        = res[12]
        self.nInf      = res[13]
        self.sInf      = res[14]
        self.objective = res[15]


        # Print solution?
        if verbose:
            print ''
            print ''.join('-' for j in range(82))
            print ' Results for problem '  + self.name
            print '   EXIT code          = ' + repr(self.exit)
            print '   INFO code          = ' + repr(self.info)
            print '   Final objective    = ' + repr(self.objective)
            print '   Total # major itns = ' + repr(self.major_itns)
            print '   Total # iterations = ' + repr(self.iterations)

            print '\n     Name state(j)      low(j)           x(j)           upp(j)        mul(j)'
            form  = ' {0:>8s}{1:8d}{2:16.6e}{3:16.6e}{4:16.6e}{5:16.6e}'
            print 'x:'
            for j in range(n):
                print form.format(pxnames[j],self.xstate[j],xlow[j],self.x[j],xupp[j],self.xmul[j])
            print 'F:'
            for j in range(nF):
                print form.format(pFnames[j],self.Fstate[j],Flow[j],self.F[j],Fupp[j],self.Fmul[j])
            print ''.join('-' for j in range(82))

        return


    def snoptb(self,funcon,funobj,**kwargs):
        """
        snoptb() calls the solver SNOPTB to solve the optimization
        problem:
              min       f_0(x)
              s.t.       [  x ]
                   bl <= [f(x)] <= bu
                         [ Ax ]
        where f(x) is a vector of smooth nonlinear constraint functions, f_0(x)
        is a smooth scalar objective function, A is a sparse matrix for the
        linear constraints.

        Keyword arguments:
        funcon   is a user-defined function that computes the constraints and the
                 Jacobian
        funobj   is a user-defined function that computes the objective and its
                 gradient

        name     is the name of the problem (default '')
        m        is the number of constraints
        n        is the number of variables
        ne       is the number of elements in J
        nnCon    is the number of nonlinear constraints
        nnJac    is the number of nonlinear Jacobian variables
        nnObj    is the number of nonlinear objective variables

        iObj     indicates which row of J is the free row containing the linear
                 objective vector
        ObjAdd   is the constant term of the objective

        hs       are the initial states of the variables (default 0)
        x        are the initial values of x (default 0.0)
        bl       are the lower bounds of the problem (default -infinity)
        bu       are the upper bounds of the problem (default +infinity)

        valJ
        locJ
        indJ

        """
        verbose  = self.options.getOption('Verbose')
        maxTries = self.options.getOption('Max memory attempts')

        name   = kwargs.get('name','')
        m      = kwargs.get('m',0)
        n      = kwargs.get('n',0)

        NamesU = kwargs.get('Names',None)

        nnCon  = kwargs.get('nnCon',0)
        nnJac  = kwargs.get('nnJac',0)
        nnObj  = kwargs.get('nnObj',0)

        iObj   = kwargs.get('iObj',0)
        ObjAdd = kwargs.get('ObjAdd',0.0)

        hs     = kwargs.get('hs',None)
        x      = kwargs.get('x0',None)
        bl     = kwargs.get('bl',None)
        bu     = kwargs.get('bu',None)

        usrJ   = kwargs.get('J',None)
        if   type(usrJ) is np.ndarray:
            if   usrJ.ndim == 1:
                locJ   = kwargs.get('locJ',None)
                indJ   = kwargs.get('indJ',None)
                J      = usrJ
                ne     = J.size
            elif usrJ.ndim == 2:
                J      = sp.csc_matrix(usrJ)
                J.sort_indices()

                indJ   = J.indices
                locJ   = J.indptr
                ne     = J.nnz
                J      = J.data
            else:
                return 99

        elif type(usrJ) is sp.csc_matrix:
            usrJ.sort_indices()
            indJ   = usrJ.indices
            locJ   = usrJ.indptr
            ne     = usrJ.nnz
            J      = usrJ.data
        else:
            return 99

        if indJ is None or locJ is None or J is None:
            return 99


        indJ = indJ + 1
        locJ = locJ + 1


        if NamesU is None:
            nNames = 1
            Names  = None
        else:
            try:
                nNames = NamesU.shape[0]
                if nNames != n+m:
                    nNames = 1
                    Names  = None
                else:
                    if NamesU.dtype.char == 'c':
                        Names  = NamesU
                    else:
                        Names = convertC(NamesU,nNames)
            except:
                nNames = 1
                Names  = None


        # Copy python-defined options to workspace
        #    This should also set iPrint,iSumm to
        #    defaults or user-def values.
        self.options.copyOptions(self.iw,self.rw)


        # Open files for printing
        iPrint  = self.options.getOption('Print unit')
        iSumm   = self.options.getOption('Summary unit')
        prtfile = self.options.getOption('Print file')
        sumfile = self.options.getOption('Summary file')
        snopt.snopen_wrap(iPrint,prtfile,iSumm,sumfile)

        # Print the title
        snopt.sntitle_wrap(self.iw)


        # Read specs file if one was given
        info = 0
        spcfile = self.options.getOption('Specs file')
        if spcfile is not None:
            iSpecs  = self.options.getOption('Specs unit')
            if iSpecs > 0:
                info = snopt.snspec_wrap(iSpecs,spcfile,self.cw,self.iw,self.rw)
            else:
                info = 131
            if info != 101 and info != 104:
                print ' Specs read failed', info
                return info


        # Check memory
        neG = nnCon*nnJac
        [info,mincw,miniw,minrw] = snopt.snmemb_wrap(m,n,ne,neG,\
                                                     nnCon,nnJac,nnObj,\
                                                     self.cw,self.iw,self.rw)
        self.resizeW(mincw,miniw,minrw)


        # Check problem info
        inf = self.options.getOption('Infinite bound')

        if hs is None:
            if verbose:
                print ' --> Initial states not provided; setting to 0'
            hs =  np.zeros(n+m,'i')
        if x is None:
            if verbose:
                print ' --> Initial x not provided; setting to 0'
            x =  np.zeros(n+m,float)
        if bl is None:
            if verbose:
                print ' --> Lower bounds not provided; setting to -inf'
            bl = -inf*np.ones(n+m,float)
        if bu is None:
            if verbose:
                print ' --> Upper bounds not provided; setting to +inf'
            bu =  inf*np.ones(n+m,float)

        pi = np.zeros(m,float)
        rc = np.zeros(n+m,float)

        # Solve problem
        Start = self.options.getOption('Start type')

        lencu = 1
        leniu = 1
        lenru = 1

        cu    = np.array(['        ']*lencu,'c')
        iu    = np.zeros(leniu,'i',order='F')
        ru    = np.zeros(lenru,float,order='F')


        if verbose:
            print ''.join('-' for j in range(82))
            print ' SNOPT python interface   (' + date + ')'
            print '   Problem: ' + name
            print '   # variables = %d; # constraints = %d \n' %(n,m)
            print '     Name state(j)      low(j)           x(j)           upp(j)        mul(j)'
            form  = ' {0:>8s}{1:8d}{2:16.6e}{3:16.6e}{4:16.6e}{5:16.6e}'

            pNames = convertS(Names,nNames,n+m)

            for j in range(n+m):
                print form.format(pNames[j],hs[j],bl[j],x[j],bu[j],rc[j])
            print ''.join('-' for j in range(82))

        if Names is None:
            nNames = 1
            Names  = np.array(['        '],dtype='c')

        count = 1
        while True:
            res = snopt.snoptb_wrap(Start,nnCon,nnObj,nnJac, \
                                    iObj, ObjAdd, name, \
                                    funcon, funobj, \
                                    J, indJ, locJ, \
                                    bl, bu, Names, \
                                    hs, x, pi, rc, \
                                    cu, iu, ru, self.cw, self.iw, self.rw)
            if  res[4]/10 == 8:
                count += 1
                if count > maxTries:
                    print ' Could not allocate memory for SNOPT'
                    return info
                self.resizeW(res[7],res[8],res[9])
            else:
                break

        # Results
        # res[0]  = hs
        # res[1]  = x
        # res[2]  = pi
        # res[3]  = rc
        # res[4]  = info
        # res[5]  = itn
        # res[6]  = mjritn
        # res[7]  = mincw
        # res[8]  = miniw
        # res[9]  = minrw
        # res[10] = nS
        # res[11] = nInf
        # res[12] = sInf
        # res[13] = objective


        # Finish up
        snopt.close_wrap(iPrint,iSumm)

        # Return solution
        self.name      = name
        self.states    = res[0]
        self.x         = res[1]
        self.z         = res[3]
        self.exit      = (res[4]/10) * 10
        self.info      = res[4]
        self.iterations= res[5]
        self.major_itns= res[6]
        self.nS        = res[10]
        self.nInf      = res[11]
        self.sInf      = res[12]
        self.objective = res[13]

        # Print solution?
        if verbose:
            print ''.join('-' for j in range(82))
            print ' Results for problem '  + self.name
            print '   EXIT code          = ' + repr(self.exit)
            print '   INFO code          = ' + repr(self.info)
            print '   Final objective    = ' + repr(self.objective)
            print '   Total # major itns = ' + repr(self.major_itns)
            print '   Total # iterations = ' + repr(self.iterations)

            print '\n     Name state(j)      low(j)           x(j)           upp(j)        mul(j)'
            form  = ' {0:>8s}{1:8d}{2:16.6e}{3:16.6e}{4:16.6e}{5:16.6e}'

            for j in range(n+m):
                print form.format(pNames[j],self.states[j],bl[j],self.x[j],bu[j],self.z[j])
            print ''.join('-' for j in range(82))

        return


    def resizeW(self,mincw,miniw,minrw):
        if mincw > self.lencw:
            tmp = np.copy(self.cw)
            self.cw = np.array(['        ']*mincw,dtype='c')
            self.cw[0:self.lencw] = tmp[0:self.lencw]
            del tmp
            self.lencw = mincw
            self.iw[6] = mincw

        if miniw > self.leniw:
            self.iw.resize(miniw)
            self.leniw = miniw
            self.iw[4] = miniw

        if minrw > self.lenrw:
            self.rw.resize(minrw)
            self.lenrw = minrw
            self.iw[2] = minrw



def convertC(array,n):
    try:
        c_array = np.array(['        ']*n,dtype='c')
        for j in range(n):
            s = ''
            for sj in array[j]:
                s += sj
            c_array[j] = s
        return c_array
    except:
        n = 1
        return None


def convertS(array,narray,n):
    try:
        if array is None:
            return [ repr(j) for j in range(n) ]

        s_array = [ repr(j) for j in range(n) ]
        if array.dtype.char == 'S':
            for j in range(narray):
                s_array[j] = array[j]
            return s_array

        elif array.dtype.char == 'c':
            for j in range(narray):
                s = ''
                for k in range(array[j].shape[0]):
                    s += array[j][k]
                s_array[j] = s
            return s_array

        return [ repr(j) for j in range(n) ]
    except:
        return [ repr(j) for j in range(n) ]


if __name__ == '__main__':
    myProb = SNOPT_solver()
    print myProb
