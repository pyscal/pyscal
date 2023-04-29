
# pyscal - python Structural Environment Calculator

Complete documentation with examples available [here](https://pyscal.org/).

**pyscal** is a python module for the calculation of local atomic structural environments including [Steinhardt's bond orientational order parameters](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784) during post-processing of atomistic simulation data. The core functionality of pyscal is written in C++ with python wrappers using [pybind11](https://pybind11.readthedocs.io/en/stable/intro.html) which allows for fast calculations with possibilities for easy expansion in python.

Steinhardt's order parameters are widely used for [identification of crystal structures](https://aip.scitation.org/doi/full/10.1063/1.4774084). They are also used to identify if an atom is [solid or liquid](https://link.springer.com/chapter/10.1007/b99429).


## Installation

pyscal can be installed directly using conda by the following statement-

```
    conda install -c conda-forge pyscal
```

**From repository**

pyscal can be built from the repository by-

```
    git clone https://github.com/pyscal/pyscal.git
    cd pyscal
    python setup.py install --user
```

## Citing the work

If you use pyscal in your work, the citation of the [following article](https://joss.theoj.org/papers/10.21105/joss.01824) will be greatly appreciated:

Sarath Menon, Grisell Díaz Leines and Jutta Rogal (2019). pyscal: A python module for structural analysis of atomic environments. Journal of Open Source Software, 4(43), 1824, https://doi.org/10.21105/joss.01824
