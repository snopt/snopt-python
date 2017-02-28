#!/usr/bin/env python

from distutils.core import setup
import os,sys

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

#    sqic_lib_path = os.getenv('SQICLIB')
#    if ( sqic_lib_path == None ):
#        sqic_lib_path=os.getcwd()+'/..'

#    sqic_mod_path = os.getenv('SQICMOD')
#    if ( sqic_mod_path == None ):
#        sqic_mod_path=os.getcwd()+'/..'

    snopt7_lib_path = os.getenv('SNOPT7LIB')
    if ( snopt7_lib_path == None ):
        snopt7_lib_path=os.getcwd()+'/..'

    config = Configuration('solvers', parent_package, top_path )

#    config.add_extension(name='sqic_python',
#                         sources=['sqic/sqic_python.pyf','sqic/sqic_python.f90'],
#                         library_dirs=sqic_lib_path,
#                         libraries=['sqic'],
#                         extra_f90_compile_args=[sqic_mod_path])

    config.add_extension(name='snopt7_python',
                         sources=['snopt/snopt7_python.pyf','snopt/snopt7_python.f90'],
                         library_dirs=snopt7_lib_path,
                         libraries=['snopt7'])

    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
