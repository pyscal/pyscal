
# pyscal - python Structural Environment Calculator

> [!WARNING]
> ## 📢 This repository is no longer maintained
>
> This repository contains **pyscal v2**. It is no longer actively
> developed and will not receive bug fixes or new features.
>
> The current version of pyscal is **v4**, developed at
> [**pyscal/pyscal3**](https://github.com/pyscal/pyscal3).
>
> - Install: `pip install pyscal3`  *or*  `conda install -c conda-forge pyscal3`
> - Then: `import pyscal`  *or*  `import pyscal3`  (both work)
> - Documentation: <https://pyscal.org>
>
> New users should start with v4. This repository is kept online so that
> existing references and citations continue to resolve.

---

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

## Works using pyscal

For a complete list of publications which used pyscal, see [here](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=315020929885190486&as_sdt=5).
