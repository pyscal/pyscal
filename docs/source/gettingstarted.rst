Getting started
===============

Trying pyscal
----------------
You can try some of the examples provided with pyscal without installing the package using binder. Please use `this link <https://mybinder.org/v2/gh/srmnitc/pyscal/master?filepath=examples%2F>`_ to try the package. 

Installation
------------

**From the repository**

Dependencies for the C++ part  

* `cmake <https://cmake.org/>`_  
* c++ 11  

Dependencies for the python part

* `numpy <https://numpy.org/>`_  

We recommend using a conda environment to install the program. A good guide on managing conda environments is available `here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_.

pyscal works with both python3 and python2. A python3 conda environment can be created by,  

.. code:: console
    
    conda create -n myenv python=3

Once created, the environment can be activated using,  

.. code:: console
    
    conda activate myenv

In case there is no cmake and c++11, these can be installed using,  

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
