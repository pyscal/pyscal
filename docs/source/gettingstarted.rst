Installation
============

Trying pyscal
----------------
You can try some examples provided with pyscal using `Binder <https://mybinder.org/>`_  without installing the package. Please use `this link <https://mybinder.org/v2/gh/srmnitc/pyscal/master?filepath=examples%2F>`_ to try the package.

Supported operating systems
---------------------------
pyscal can be installed on both Linux and Mac OS based systems. For Windows systems, we recommend using Bash on Windows. Please check `this <https://lammps.sandia.gov/doc/Howto_bash.html>`_ tutorial
on how to set up Bash on Windows.

Installation using `conda <https://anaconda.org>`_
--------------------------------------------------

pyscal can be installed directly using `Conda <https://docs.conda.io/en/latest/>`_ from the `conda-forge channel <https://conda-forge.org/>`_ by the following statement-

.. code:: console

    conda install -c conda-forge pyscal


This is the recommended way to install if you have an `Anaconda <https://www.anaconda.com/>`_ distribution.

The above command installs the `latest release version <https://github.com/srmnitc/pyscal/releases>`_ of pyscal.
An always updated version can be installed using-

.. code:: console

    conda install -c pyscal pyscal


pyscal is no longer maintained for Python 2.

Quick installation
------------------

pyscal can be installed using the following steps-

* Download an archive of the pyscal library from `here <https://pyscal.readthedocs.io/en/latest/download.html>`_.

* Extract the downloaded version. From the extracted folder, run, ``python setup.py install --user``

.. tip::

    Pyscal can be installed system-wide using ``python setup.py install``.


Installation from the repository
--------------------------------

pyscal can be built from the repository by-

.. code:: console

    git clone https://github.com/srmnitc/pyscal.git
    cd pyscal
    python setup.py install --user

Using a conda environment
-------------------------

pyscal can also be installed in a conda environment, making it easier to manage dependencies.
A python3 Conda environment can be created by,

.. code:: console

    conda create -n myenv python=3

Once created, the environment can be activated using,

.. code:: console

    conda activate myenv

In case Cmake and C++11 are not available, these can be installed using,

.. code:: console

    (myenv) conda install -c anaconda gcc
    (myenv) conda install -c anaconda cmake

Now the pyscal repository can be cloned and the module can be installed. Python dependencies are installed automatically.

.. code:: console

    (myenv) git clone https://github.com/srmnitc/pyscal.git
    (myenv) cd pyscal
    (myenv) python setup.py install

.. tip::

    A good guide on managing Conda environments is available `here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.


Dependencies
------------

Dependencies for the C++ part

* `Cmake <https://cmake.org/>`_
* C++ 11

Dependencies for the python part

* `numpy <https://numpy.org/>`_

Optional dependencies

* `pytest <https://docs.pytest.org/en/latest/>`_
* `matplotlib <https://matplotlib.org/>`_

Tests
-----

In order to see if the installation worked, the following commands can be tried-

.. code:: python

    import pyscal.core as pc
    pc.test()

The above code does some minimal tests and gives a value of ``True`` if pyscal was installed successfully. However, pyscal also contains automated tests which
use the `pytest <https://docs.pytest.org/en/latest/>`_ python library, which can be installed by ``pip install pytest``.
The tests can be run by executing the command ``pytest tests/`` from the main code directory.

It is good idea to run the tests to check if everything is installed properly.
