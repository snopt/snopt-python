#!/usr/bin/env python

try:
    from snopt_solver import SNOPT_solver
    __all__ = ['SNOPT_solver']
except:
    __all__ = []
