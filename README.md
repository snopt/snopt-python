snopt-python
================

This package provides Python support for the nonlinear optimization code SNOPT.  This package does not include the SNOPT libraries, only the python interface.

To use:
Export the location of the SNOPT to variable $SNOPT7, e.g.,
```
export SNOPT7=$HOME/snopt7
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
