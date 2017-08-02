snopt-python
================

This package provides Python support for the nonlinear optimization code SNOPT.  This package does not include the SNOPT libraries, only the python interface.

Note: Use this interface at your own risk.  I haven't had the time to do a lot of testing on it yet.

To use:
Set the environment variable SNOPT7LIB to the location of your SNOPT libraries, e.g.,
```
export SNOPT7LIB=$HOME/snopt7/lib
```
On Linux systems, you should also make sure that this location is on LD_LIBRARY_PATH.

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
