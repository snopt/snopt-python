snopt-python
================

This package provides Python support for the nonlinear optimization code SNOPT.  This package does not include the SNOPT libraries, only the python interface.

To use:
Set the environment variable SNOPT7LIB to the location of your SNOPT libraries, e.g.,
```
export SNOPT7LIB=$HOME/snopt7/lib
```

Build, install,
```
python setup.py build
python setup.py install --user
```

Run an example:
```
python dieta.py
```

For info on SNOPT, http://ccom.ucsd.edu/~optimizers
