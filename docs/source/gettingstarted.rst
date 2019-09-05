Getting started
===============

Trying ``pyscal``
----------------
You can try some of the examples provided with ``pyscal`` without installing the package using binder. Please use `this link <https://mybinder.org/v2/gh/srmnitc/pyscal/master?filepath=examples%2F>`_ to try the package. 

Installation
------------

**From the repository**

``pyscal`` needs ``cmake`` to compile the C++ part of the code. Please check `here <https://cmake.org/install/>`_ for documentation on how to install cmake.

Once cmake is available, clone the **pyscal** repository by `git clone https://github.com/srmnitc/pyscal.git`.
After cloning the repository, **pyscal** can be installed by running ``python setup.py install`` from main code directory. It can be uninstalled by ``pip uninstall pathsampling``. All the dependencies of **pyscal** are installed automatically.

Tests
-----
**pyscal** also contains automated tests which use the `pytest <https://docs.pytest.org/en/latest/>`_ python library, which can be installed by ``pip install pytest``. The tests can be run by executing the command ``pytest tests/`` from the main code directory.

It is good idea to run the tests to check if everything is installed properly.

