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
date: 05 August 2019
bibliography: paper.bib
---


# Summary

Atomistic simulation methods provide valuable insight into materials and their properties from an atomistic perspective. These methods often produce large amounts of data in terms of atomic positions over time, the analysis of which is crucial to study underlying mechanisms. Many computational analysis methods have been developed to determine the local environment of a particle, which is used for various purposes such as the distinction of a particle as solid or liquid, to identify the crystal structure they belong to, to identify defects in a crystalline state, to name a few [@Stukowski:2012]. Bond orientational order parameters introduced by Steinhardt and coworkers [@Steinhardt:1983], have been used extensively in this regard [@Mickel:2013 and the references within]. These parameters, which are rotationally invariant are defined by,


$$ \bar{q}_l (i) =  \Big(  \frac{2\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\pmb{r}_{ij}) \Big|^2 \Big )^{\frac{1}{2}} $$

where $Y_{lm}$ are the spherical harmonics and $N(i)$ is the number of nearest neighbours of particle $i$. $\pmb{r}_{ij}$ is the displacement vector connecting particles $i$ and $j$, and $l$ and $m$ are both intergers with $m \in [-l,+l]$.  

However, at high temperatures, atomic positions are subject to thermal noise which leads to difficulty in the resolution of crystal structures. The averaged version of the order parameters ($\bar{q}_l$), was introduced [@Lechner:2008] in order to counter this problem. ($\bar{q}_l$), calculated by averaging over the bond order parameters of nearest neighbour particles, was shown to improve the resolution of crystal structures by few order of magnitude. Another problem at high temperatures is the use of a fixed cutoff radius which leads to ambiguity in definition of the local neighbourhood of a particle. A parameter free definition of the local neighbourhood using voronoi tessellation was proposed[@Mickel:2013] and it was shown that the distinction of crystal structures can be improved by weighting the contribution of each neighbour to $q_l$ given by,         

$$ {q}_l (i) =  \Big(  \frac{2\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{N(i)} \sum_{j=1}^{N(i)} \frac{A_{ij}}{A} Y_{lm}(\pmb{r}_{ij}) \Big|^2 \Big )^{\frac{1}{2}} $$

where $A_{ij}$ is area of the voronoi facet between atoms $i$ and $j$, and $A$ is the total surface area of the voronoi polyhedra. Futher approaches include the using higher powers, $A_{ij}^\alpha$, with $\alpha=2,3$ to weight the contribution to $q_l$[@Haeberle:2019] and using an adaptive cutoff to identify neighbours similar to the approach followed in adaptive common neighbour analysis[@Stukowski:2012].

Originally inspired by [``BondOrderAnalysis``](https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html), a C++ code for calculation of bond order parameters, ``pybop`` is a python module designed for the computation of bond orientational order parameters during postprocessing of simulation data. ``pybop`` brings together the various approaches described above in a single python module while ensuring the speed and efficiency of calculations as the core code for ``pybop`` is written in C++, which is then ported to python using ``pybind11`` [@Jakob:2016]. The module includes the following functionality-   

 * calculation of bond order parameters ($q_{l}$) and the averaged version ($\bar{q}_{l}$).
 * weighted bond order parameters using face area of voronoi polyhedra, the calculation of which is  enabled using Voro++ code [@Rycroft:2009] integrated into ``pybop``.
 * distinction of atoms as solid or liquid based on $q_6$ [@Auer:2005].
 * improvement of local neighbourhood description using an adaptive cutoff.
 * clustering algorithm to find clusters of particles based on a user defined property.
 * other quantities like radial distribution function, voronoi volume of individual particles and coordination number.  

``pybop`` uses a list of particle positions and simulation box vectors as input. It can also read in output files containing atomistic simulation data in the LAMMPS [@Plimpton:1995] [dump format](https://lammps.sandia.gov/doc/dump.html) and [POSCAR](https://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html) format used by _ab initio_ simulation package [VASP](https://www.vasp.at/). Furthermore, ``pybop`` can read in orthogonal as well as triclinic simulation cells. Additionally, ``pybop`` offers a set of support tools including those for reading, writing and splitting of atomic trajectory files and creation of perfect crystal structures. Various features of ``pybop`` including documentation and example usage is available in the [pybop website](https://pybop.readthedocs.io/en/latest/). Examples including installation of the package to calculation of bond orientational order parameters are provided. Alternatively, ``pybop`` provides a [binder environment](https://mybinder.org/v2/gh/srmnitc/pybop/master?filepath=examples%2F) to try the example cases before installation.



# Acknowledgements

# References