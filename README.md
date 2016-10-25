snopt-python
================

This package provides Python support for the nonlinear optimization code SNOPT.  This package does not include the SNOPT libraries, only the python interface.

To use:
This package assumes the SNOPT libraries are located in $SNOPT7/lib, where $SNOPT7 is user-defined (The 'lib' part is fixed -- if you know what you're doing, you could modify snopt-python/snopt7/setup.py to change it).

For example if your libraries are in $HOME/snopt7/lib, then eeport the location of the SNOPT to variable $SNOPT7, e.g.,
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
