#!/usr/bin/env python

from distutils.core import setup
import os,sys

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    snopt7_path = os.getenv('SNOPT7')
    if ( snopt7_path == None ):
        snopt7_path=os.getcwd()+'/..'
    py_path   = snopt7_path+'/lib'

    config = Configuration('snopt7', parent_package, top_path )
    config.add_extension(name='snopt7_python',
                         sources=['snopt7_python.pyf','snopt7_python.f90'],
                         library_dirs=py_path,
                         libraries=['snopt7'])
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())

