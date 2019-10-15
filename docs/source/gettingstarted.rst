Getting started
===============

Trying pyscal
----------------
You can try some examples provided with pyscal using `Binder <https://mybinder.org/>`_  without installing the package. Please use `this link <https://mybinder.org/v2/gh/srmnitc/pyscal/master?filepath=examples%2F>`_ to try the package.

Installation using `conda <https://anaconda.org>`_
--------------------------------------------------

pyscal can be installed directly using `Conda <https://docs.conda.io/en/latest/>`_ by the following statement-

.. code:: console

    ``conda install -c pyscal pyscal``


This is the recommended way to install if you have an `Anaconda <https://www.anaconda.com/`_ distribution.

Quick Installation from repository
----------------------------------

Quick installation can be useful when you have administrator privileges, for example, on your personal computer, laptop etc.
It is recommended to install using Conda as described above. If building from repository is preferred,
directly clone the repository and set up an environment as discussed in the `Installation from the repository section <https://pyscal.readthedocs.io/en/latest/gettingstarted.html#installation-from-the-repository>`_

* Download an archive of the pyscal library from `here <https://pyscal.readthedocs.io/en/latest/download.html>`_.

* Extract the downloaded version. From the extracted folder, run, ``python setup.py install``.


Installation from the repository
--------------------------------

Dependencies for the C++ part

* `Cmake <https://cmake.org/>`_
* C++ 11

Dependencies for the python part

* `numpy <https://numpy.org/>`_

Optional dependencies

* `pytest <https://docs.pytest.org/en/latest/>`_
* `matplotlib <https://matplotlib.org/>`_

We recommend using a Conda environment to install the program. A good guide on managing Conda environments is available `here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

pyscal works with both python3 and python2. A python3 Conda environment can be created by,

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


Tests
-----
pyscal also contains automated tests which use the `pytest <https://docs.pytest.org/en/latest/>`_ python library, which can be installed by ``pip install pytest``. The tests can be run by executing the command ``pytest tests/`` from the main code directory.

It is good idea to run the tests to check if everything is installed properly.
