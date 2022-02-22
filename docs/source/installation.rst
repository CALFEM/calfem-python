Installation instructions
^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest way to install CALFEM for Python is via pip. 
This procedure will ensure that all required dependencies are fulfilled.

This can be achieved by executing the following command::

    pip install calfem-python

or::

    sudo pip install calfem-python

to install system-wide (not recommended if your system used a lot of python dependencies)::

    pip install -u calfem-python

to install just for your own user. You can use the argument `--user` which is 
same as `-u`. If you want to specify your Python version, use the command like 
the following::

    python3.9 -m pip install --user calfem-python

where python3.9 is the Python version you want to install CALFEM for
Python. Change this command with your preferable version. This last command is
the preferred one according to the Python community.

Installing in a Anaconda environment
------------------------------------

If you don't want to install all packages directly in the base Anaconda environment you could create a new enviroment for CALFEM::

    conda create -n calfem-dev python=3.9
    conda activate calfem-dev 

Now we can install CALFEM for Python and its dependencies directly into this environment using Pip.::

    pip install calfem-python


