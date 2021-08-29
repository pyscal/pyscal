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
version](https://github.com/srmnitc/pyscal/releases) of pyscal and works
on all three operating systems.

```{warning}
pyscal is no longer maintained for Python 2. Although quick installation method might work for Python 2, all features may not work as expected.
```

### Quick installation

pyscal can be installed using the following steps-

-   Download an archive of the pyscal library from
    [here](https://pyscal.readthedocs.io/en/latest/download.html).
-   Extract the downloaded version. From the extracted folder, run,
    `python setup.py install --user`

```{note}
Pyscal can be installed system-wide using `python setup.py install`.
```


### Installation using pip

pyscal is not available on pip directly. However pyscal can be installed
using pip by

``` console
pip install git+https://github.com/srmnitc/pyscal
```

### Installation from the repository

pyscal can be built from the repository by-

``` console
git clone https://github.com/srmnitc/pyscal.git
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

In case Cmake and C++11 are not available, these can be installed using,

``` console
(myenv) conda install -c anaconda gcc
(myenv) conda install -c anaconda cmake
```

Now the pyscal repository can be cloned and the module can be installed.
Python dependencies are installed automatically.

``` console
(myenv) git clone https://github.com/srmnitc/pyscal.git
(myenv) cd pyscal
(myenv) python setup.py install
```

``` {tip}
A good guide on managing Conda environments is available
[here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).
```

### Dependencies

Dependencies for the C++ part

-   [Cmake](https://cmake.org/)
-   C++ 11

Dependencies for the python part

-   [numpy](https://numpy.org/)
-   [ase](https://wiki.fysik.dtu.dk/ase/)

Optional dependencies

-   [pytest](https://docs.pytest.org/en/latest/)
-   [matplotlib](https://matplotlib.org/)

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
