Developing CALFEM for Python
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Creating an environment for development
---------------------------------------

It could be benficial to create an dedicated environment for development. This can be easily created in the Anaconda distribution.::

    conda create -n calfem-dev python=3.9
    conda activate calfem-dev 

Installing dependencies
-----------------------

To develop for CALFEM for Python you will need to install the required dependencies first. The easiest way to do this is to install the pip version of CALFEM::

    pip install calfem-python

Creating a fork in GitHub
-------------------------

To be able to submit code changes and addition it is a good idea to create a fork on GitHub. A fork is your own version of CALFEM for Python in which you can track changes. From the fork you can also easily create a pull request, that is a suggested change that you can submit to CALFEM for Python.

Please see

https://reflectoring.io/github-fork-and-pull/

for more information on how to create and manage forks.

Checking out code from github
-----------------------------

It is also possible to check out a local version of the code from the command line. The master or develop branches can be checkout with the following procedures:

Cloning the master branch on your local computer.::

    git clone https://github.com/CALFEM/calfem-python.git

Cloning the develop branch on your local computer.::

    git clone https://github.com/CALFEM/calfem-python.git
    git checkout develop

It is also possible to use the command line tools to clone your github fork:

    git clone https://github.com/USERNAME/calfem-python.git


Using GitHub Desktop to manage your code
----------------------------------------

GitHub Desktop is a graphical client for GitHub that makes the procedure working with git-repos much more easy. More information on how to use GitHub Desktop can be found here:

https://desktop.github.com/

https://docs.github.com/en/desktop


Modifying the Python search PATH
--------------------------------

To work on the checked out source directory, set the PYTHONPATH environment variable to the root of the source directory.

On Windows::

    set PYTHONPATH=C:\[Path to source directory]

On other platforms::

    export PYTHONPATH=[Path to source directory]

Guidelines
----------

If you want to develop standalone element routines for CALFEM it is probarbly a good idea to develop them standalone and then submit a change request for inclusion in the calfem.extension module.

Element routines should be named according to the following rule:

    [name][1|2|3][e|s]

Where 1,2 and 3 indicates element dimensions. e - denotes that the routine returns element stiffness matrix and or force vector. s - denotes that the routines returns element forces given element nodal values.

Please look at existing routines to familiar yourself with how they are constructued. 

All routines should have documentation strings describing input and output parameters. Code comments describing the general flow should be added to the function at relevant positions.

Submit a Pull request
---------------------

If you want your changes to be included in future releases of CALFEM please create a pull request against GitHub. More information on this can be found here:

https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request






