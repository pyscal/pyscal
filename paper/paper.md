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
    orcid: 0000-0002-3048-5315
    affiliation: 1
  - name: Jutta Rogal
    orcid: 0000-0002-6268-380X
    affiliation: 1
affiliations:
 - name: Interdisciplinary centre for advanced materials simulation, Ruhr-Universit$\"{a}t$ Bochum, 44780     Bochum, Germany
   index: 1
date: 05 August 2019
bibliography: paper.bib
---


# Summary

Atomistic simulation methods are widely used in the study of condensed matter systems, often producing large amounts of data in terms of atomic positions over time. Analysis of this data to identify the local atomic environment, which is then used to determine whether a particle is part of a crystal structure or not is a common challenge [@Stukowski:2012]. A widely used approach to identify solid particles was proposed by Auer and Frenkel [@Auer:2005], based on the bond orientational order parameters introduced by Steinhardt and coworkers [@Steinhardt:1983]. These parameters ($q_l$, where $l$ is an integer from 2-12), are rotationally and translationally invariant and are based on spherical harmonics. $q_l$ take into account the local neighbourhood of an particle through the interatomic distance between a particle and its neighbors. These parameters have been used extensively, especially for the identification of crystal structures [@Mickel:2013 and references within].    

However, at high temperatures, atomic positions are subject to thermal vibrations which leads to a broad distribution of these parameters, which may make the resolution of crystal structures harder. The averaged version of the order parameters ($\bar{q}_l$), was introduced [@Lechner:2008] in order to counter this problem. $\bar{q}_l$, calculated by averaging over the bond order parameters of nearest neighbour particles, was shown to improve the resolution of crystal structures. Another problem at high temperatures is the use of a fixed cutoff radius which leads to ambiguity in the definition of the local neighbourhood of a particle. A parameter free definition of the local environment can be achieved by using a distinct cutoff radius for each particle using methods such as the SANN algorithm [@VanMeel:2012] or an adaptive cutoff method similar to to the approach followed in adaptive common neighbour analysis [@Stukowski:2012]. Additionally, Voronoi tessellation can be employed for a parameter free definition of the local neighbourhood. It was shown [@Mickel:2013] that the distinction of crystal structures can be improved by weighting the contribution of each neighbour to $q_{l}$ by the area of the Voronoi facet shared between the atom and its neighbour. Futher approaches include the use of higher powers of area to weight the contribution to $q_{lm}$[@Haeberle:2019].

Originally inspired by [``BondOrderAnalysis``](https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html) [@Lechner:2010], a C++ code for calculation of bond order parameters, ``pybop`` is a python module designed for the computation of bond orientational order parameters during post-processing of simulation data. ``pybop`` brings together the various approaches described above in a single python module while ensuring the speed and efficiency of calculations as the core code for ``pybop`` is written in C++, which is then ported to python using ``pybind11`` [@Jakob:2016]. The module includes the following functionality-   

 * calculation of bond order parameters ($q_{l}$) and the averaged version ($\bar{q}_{l}$).
 * weighted bond order parameters using face area of voronoi polyhedra, the calculation of which is  enabled using Voro++ code [@Rycroft:2009] integrated into ``pybop``.
 * distinction of atoms as solid or liquid based on $q_6$ [@Auer:2005].
 * improvement of local neighbourhood description using an adaptive cutoff.
 * clustering algorithm to find clusters of particles based on a user defined property.
 * inbuilt functions for other structural features like radial distribution function, voronoi volume of individual particles and coordination number.
 * calculation of other Voronoi tessellation based parameters such as number of vertices and face area.

``pybop`` uses a list of particle positions and simulation box vectors as input. It can also read in output files containing atomistic simulation data in the LAMMPS [@Plimpton:1995] [dump format](https://lammps.sandia.gov/doc/dump.html) and [POSCAR](https://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html) format used by _ab initio_ simulation package [VASP](https://www.vasp.at/), including orthogonal as well as triclinic simulation cells. ``pybop`` also offers a set of support tools including those for reading, writing and splitting of atomic trajectory files and creation of perfect crystal structures. Overall, ``pybop`` also acts as a tool for calculation of quantities based on atomic positions, and can be easily extended using either C++ for computationally intensive operations or python for other tasks that wrap over the existing features. Various features of ``pybop`` including documentation and example usage is available in the [pybop website](https://pybop.readthedocs.io/en/latest/). Examples including installation of the package to calculation of bond orientational order parameters are provided. Alternatively, ``pybop`` provides a [binder environment](https://mybinder.org/v2/gh/srmnitc/pybop/master?filepath=examples%2F) to try the example cases before installation.



# Acknowledgements
S.M. acknowledges a scholarship from the International Max Planck Research School for Interface Controlled Materials for Energy Conversion. G.D.L acknowledges the support by the Mexican National Council for Science and Technology (CONACYT) through project 232090 and by the German Research Foundation (DFG) through project RO 3073/6-1.

# References