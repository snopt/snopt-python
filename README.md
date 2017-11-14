snopt-python
================

This package provides Python support for the nonlinear optimization code SNOPT.  This package does not include the SNOPT libraries, only the python interface.

Note: Use this interface at your own risk.  I haven't had the time to do a lot of testing on it yet.

To use:
Set the environment variable DNOPTLIB to the location of your DNOPT libraries, for example, the following commands assume your DNOPT Fortran librarie is in a directory called "lib" in your home directory:
```
export DNOPTLIB=$HOME/lib
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
