# Getting started

## Trying pyscal

You can try some examples provided with pyscal using
[Binder](https://mybinder.org/) without installing the package. Please
use [this
link](https://mybinder.org/v2/gh/srmnitc/pyscal/master?urlpath=lab) to
try the package.

## Installation

### Supported operating systems

pyscal can be installed on Linux, Mac OS and Windows based systems.

### Installation using [conda](https://anaconda.org)

pyscal can be installed directly using
[Conda](https://docs.conda.io/en/latest/) from the [conda-forge
channel](https://conda-forge.org/) by the following statement-

``` console
conda install -c conda-forge pyscal
```

This is the recommended way to install if you have an
[Anaconda](https://www.anaconda.com/) distribution.

The above command installs the [latest release
version](https://github.com/pyscal/pyscal/releases) of pyscal and works
on all three operating systems.

```{warning}
pyscal is no longer maintained for Python 2. Although quick installation method might work for Python 2, all features may not work as expected.
```


### Installation using pip

pyscal is not available on pip directly. However pyscal can be installed
using pip by

``` console
pip install pybind11
pip install git+https://github.com/pyscal/pyscal
```

### Installation from the repository

pyscal can be built from the repository by-

``` console
git clone https://github.com/pyscal/pyscal.git
pip install pybind11
cd pyscal
python setup.py install --user
```

### Using a conda environment

pyscal can also be installed in a conda environment, making it easier to
manage dependencies. A python3 Conda environment can be created by,

``` console
conda create -n myenv python=3
```

Once created, the environment can be activated using,

``` console
conda activate myenv
```

In case C++11 is not available, these can be installed using,

``` console
(myenv) conda install -c anaconda gcc
```

Now the pyscal repository can be cloned and the module can be installed.
Python dependencies are installed automatically.

``` console
(myenv) git clone https://github.com/pyscal/pyscal.git
(myenv) conda install -c conda-forge pybind11
(myenv) cd pyscal
(myenv) python setup.py install
```

``` {tip}
A good guide on managing Conda environments is available
[here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).
```

### Dependencies

Dependencies for the C++ part

-   [pybind11](https://github.com/pybind/pybind11)
-   C++ 11

Dependencies for the python part

-   [numpy](https://numpy.org/)
-   [ase](https://wiki.fysik.dtu.dk/ase/)
-   [plotly](https://plotly.com/)
-   [ipywidgets](https://ipywidgets.readthedocs.io/en/latest/)

Optional dependencies

-   [pytest](https://docs.pytest.org/en/latest/)
-   [matplotlib](https://matplotlib.org/)
-   [LAMMPS](https://www.lammps.org/)

### Tests

In order to see if the installation worked, the following commands can
be tried-

``` python
import pyscal.core as pc
pc.test()
```

The above code does some minimal tests and gives a value of `True` if
pyscal was installed successfully. However, pyscal also contains
automated tests which use the
[pytest](https://docs.pytest.org/en/latest/) python library, which can
be installed by `pip install pytest`. The tests can be run by executing
the command `pytest tests/` from the main code directory.

It is good idea to run the tests to check if everything is installed
properly.
