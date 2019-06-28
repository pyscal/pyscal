---
title: 'pybop : A hybrid python/C++ package for calculation of bond orientational order parameters'
tags:
  - Python
  - materials science
  - bond order parameters
authors:
  - name: Sarath Menon
    orcid: 0000-0002-6776-1213
    affiliation: 1
  - name: Grisell D$\'{i}$az Leines
    orcid: 0000-0000-0000-0000
    affiliation: 1
  - name: Jutta Rogal
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: Interdisciplinary centre for advanced materials simulation, Ruhr-Universit$\"{a}t$ Bochum, 44780     Bochum, Germany
   index: 1
date: 26 June 2019
bibliography: paper.bib
---


# Summary

Bond-orientational order parameters introduced by [@Steinhardt:1983,@Auer:2005] have been widely used in distinction of different crystal structures like body-centred cubic, face-centred cubic and hexagonal closed packed structures in computational studies in materials science[@Lechner:2011,@DiazLeines:2017]. Furthermore, the bond order parameters parameters are used to distinguish solid and liquid particles in computer simulations of crystallisation during solidification[@DiazLeines:2017]. 

Although bond-orientational order parameters can be calculated using atomistic simulation program like LAMMPS [@Plimpton:1995] or other codes like, they lack in ease of use and flexibility in a post-processing environment. ``pybop`` is designed to address the above issues, by providing a python library for calculation of bond order parameters intended for post-processing of simulation data. While python provides easy extension and flexibility, ``pybop`` ensures the speed and efficiency of calculations as the core code for ``pybop`` is written in C++, which is then ported to python using ``pybind11`` [@Jakob:2016]. 

``pybop`` can perform a variety of calculations including the bond order parameters $q_{i}$ where $i = \{2,3 \to 12\}$, and the averaged versions [@Lechner:2008] which has been to improve the resolution in identification of crystal structures. Furthermore, ``pybop`` can also be used to calculate weighted $q_{i}$ [@Mickel:2013] where the contributions are weighted by the voronoi face area shared with adjacent atoms, the calculation of which is enabled with the help of Voro++ code [@Rycroft:2009] integrated into ``pybop``. Additionally, distinction of liquid and solid atoms based on these parameters [@DiazLeines:2017] is also incorporated into the package. Various other quantities like radial distribution function, coordination number and voronoi volume of individual particles, for example, can also be calculated.

``pybop`` can read in output files containing atomistic simulation data in the LAMMPS [@Plimpton:1995] [dump format](https://lammps.sandia.gov/doc/dump.html) and POSCAR format used by ab initio simulation package VASP. In order to improve the functionality, an easy interface is provided to extended the type of input file formats. In addition, both cubic and triclinic simulation cells are also supported in this package.  

Various features of ``pybop`` including documentation and example usage is available in the [pybop website](https://srmnitc.github.io/pybop/html/index.html). Examples       

## Statement of need

# Acknowledgements

# References