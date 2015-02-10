"""
  Python interface for SNOPT
  http://ccom.ucsd.edu/~optimizers

  Problem data:
   min    f(x)
   s.t.   l <= [  x ] <= u
               [ Ax ]
               [c(x)]

  Elizabeth Wong, University of California, San Diego
  Philip Gill,    University of California, San Diego

  Feb 2015
"""

import os, sys

#-------------------------------------------------------------------------------#
class SNOPT_options(object):
    '''
    SNOPT_options class:
    '''
    def __init__ (self):
        self.setup()


    def setup(self):
        self.options = {
            # [Current value, default value, type, snopt location]
            # "snopt location" is zero-index modified
            'Start type'            : ['Cold','Cold',str,-1],

            'Specs unit'            : [4,4,int,10],
            'Specs file'            : [None,None,str,-1],

            'Print unit'            : [9,9,int,11],
            'Print file'            : ['SNOPT.out','SNOPT.out',str,-1],
            'Print frequency'       : [None,None,int,60],
            'Print level'           : [None,None,int,92], # minor print level

            'Summary unit'          : [6,6,int,12],
            'Summary file'          : ['','',str,-1],
            'Summary frequency'     : [None,None,int,61],

            'Solution file'         : [None,None,int,130],
            'Solution print'        : [None,None,bool,83],
            'Major print level'     : [None,None,int,91],
            'Minor print level'     : [None,None,int,92],

            'Sticky parameters'     : [None,None,int,115],
            'Suppress'              : [None,None,int,80],
            'Time limit'            : [None,None,float,78],
            'Timing level'          : [None,None,int,181],
            'System information'    : [None,None,int,70],
            'Verify level'          : [None,None,int,77],

            'Max memory attempts'   : [10,10,int,-1],

            'Total character workspace' : [None,None,int,6],
            'Total integer workspace'   : [None,None,int,4],
            'Total real workspace'      : [None,None,int,2],

            #'Problem minmax'        : ['Minimize','Minimize',str,86],
            'Proximal point'        : [None,None,int,78],
            #'QP solver'             : [None,None,str,55],    # Cholesky/CG/QN

            'Major feasibility'     : [None,None,float,56],  #tolCon
            'Major optimality'      : [None,None,float,51],  #tolOptNP
            'Minor feasibility'     : [None,None,float,55],  #tolx
            'Minor optimality'      : [None,None,float,51],  #tolOptQP
            'Minor phase1'          : [None,None,float,50],  #tolOptFP
            'Feasibility tolerance' : [None,None,float,55],  #tolx
            'Optimality tolerance'  : [None,None,float,51],  #tolOptQP

            'Iteration limit'       : [None,None,int,88],  #itnlim
            'Major iterations'      : [None,None,int,89],  #mMajor
            'Minor iterations'      : [None,None,int,90],  #mMinor

            'CG tolerance'          : [None,None,float,53],
            'CG preconditioning'    : [None,None,int,76],
            'CG iterations'         : [None,None,int,96],

            'Crash option'          : [None,None,int,87],
            'Crash tolerance'       : [None,None,float,61],

            'Debug level'           : [None,None,int,84],

            'Derivative level'      : [None,None,int,69],
            'Derivative linesearch' : [None,None,int,75],
            'Derivative option'     : [None,None,int,103],

            'Elastic objective'     : [None,None,int,72],
            'Elastic mode'          : [None,None,int,55],
            'Elastic weight'        : [None,None,float,87],
            'Elastic weightmax'     : [None,None,float,89],

            'Hessian frequency'     : [None,None,int,63],
            'Hessian flush'         : [None,None,int,65],
            'Hessian type'          : [None,None,int,71],
            'Hessian updates'       : [None,None,int,53],

            'Infinite bound'        : [1.0e+20,1.0e+20,float,69],
            'Major step limit'      : [None,None,float,79],
            'Unbounded objective'   : [None,None,float,70],
            'Unbounded step'        : [None,None,float,71],

            'Linesearch tolerance'  : [None,None,float,83],
            'Linesearch debug'      : [None,None,int,81],

            #'LU type'               : [None,None,str,79],   #partial/complete/rook
            'LU swap'               : [None,None,float,64],
            'LU factor tolerance'   : [None,None,float,65],
            'LU update tolerance'   : [None,None,float,66],
            'LU density'            : [None,None,float,157],
            'LU singularity'        : [None,None,float,[154,155]],

            'New superbasics'       : [None,None,int,94],
            'Partial pricing'       : [None,None,int,100],
            'Penalty parameter'     : [None,None,float,88],
            'Pivot tolerance'       : [None,None,float,59],
            'Reduced Hessian limit' : [None,None,int,51],
            'Superbasics limit'     : [None,None,int,52],

            'Scale option'          : [None,None,int,74],
            'Scale tolerance'       : [None,None,float,91],
            'Scale print'           : [None,None,int,82],

            'Verbose'               : [False,False,bool,-1]  # python verbose
        }


    def __str__(self):
        text = '\nSNOPT options \n'
        keywords = self.options.keys()
        keywords.sort()
        for key in keywords:
            value = self.options[key][0]
            if value is None:
                text += ' ' + str(key) + ': undefined\n'
            else:
                text += ' ' + str(key) + ': ' + str(self.options[key][0]) + '\n'
        text += '\n' + 'SNOPT will set any parameters that are undefined to defaults.\n'
        text += 'Please refer to SNOPT documentation for details.\n'
        return text


    def printOptions(self):
        print self.__repr__()


    def setOption(self,name,value):
        try:
            if type(value) is self.options[name][2]:
                self.options[name][0] = value
        except:
            raise RuntimeError, 'Incorrect option keyword or type'
            return None

    def getOption(self,name):
        try:
            value = self.options[name][0]
            return value
        except:
            raise RuntimeError, 'Incorrect option keyword'
            return None


    def copyOptions(self,iw,rw):
        for key in self.options:
            if key == 'Solution print':
                k = self.options[key][3]
                if self.options[key][0] is True:
                    iw[k] = 2
                else:
                    iw[k] = 0
            else:
                k = self.options[key][3]
                if self.options[key][0] is not None and k >= 0:
                    if   self.options[key][2] is int:
                        iw[k] = self.options[key][0]
                    elif self.options[key][2] is float:
                        rw[k] = self.options[key][0]
                    elif self.options[key][2] is bool:
                        if self.options[key][0] is True:
                            iw[k] = 1
                        else:
                            iw[k] = 0


    def resetOption(self,name):
        try:
            self.options[name] = self.options[name][1]
        except:
            raise RuntimeError, 'Incorrect option keyword'
            return None



if __name__ == '__main__':
   print 'SNOPT options class...\n'
   myOptions = SNOPT_options()
   print myOptions
