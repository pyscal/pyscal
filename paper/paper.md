---
title: 'pyscal : python Structural Environment Calculator'
tags:
  - Python
  - materials science
  - bond orientational order parameters
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

Atomistic simulation methods are widely used in the study of condensed matter systems, often producing large amounts of data in terms of atomic positions over time. The analysis of this data to identify the local atomic environment and determine whether a particle is part of a crystal structure is a common challenge [@Stukowski:2012].

A widely used approach to identify solid particles from liquid was proposed by Auer and Frenkel [@Auer:2005],  which measures the structural correlations between a particle and its neighbors and is based on the bond orientational order parameters introduced by Steinhardt et al. [@Steinhardt:1983]. These parameters are defined by spherical harmonics in the local neighborhood of a particle within a cutoff distance, and are rotationally and translationally invariant. Steinhardt parameters have therefore been extensively used for the identification of crystal structures [@Mickel:2013 and references within]. At high temperatures, however, the atomic positions are subject to thermal vibrations which leads to a broad distribution of these parameters, and may result in poor resolution of the crystal structures. An averaged version of the Steinhardt bond order parameters introduced by Lechner and Dellago [@Lechner:2008] was shown to improve accuracy in the identification of crystal structures at finite temperatures. Another common problem at high temperatures is the use of a fixed cutoff distance which can lead to an ambiguity in the definition of the local neighborhood of a particle. An adaptive definition of the local environment can be achieved by using a distinct cutoff radius for each particle using methods such as the SANN algorithm [@VanMeel:2012] or the adaptive common neighbor analysis [@Stukowski:2012]. Additionally, Voronoi tessellation can be employed for a parameter-free definition of the local neighborhood of a particle. Mickel et al. showed [@Mickel:2013] that the identification of crystal structures can be improved by weighting the contribution of each neighbor to the Steinhardt parameters by the area of the Voronoi facet shared between the atom and its neighbors. Further approaches to improve the resolution of the crystal structures make use of higher exponents of the area to weight the bond orientational order parameters[@Haeberle:2019].

Originally inspired by [``BondOrderAnalysis``](https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html) [@Lechner:2010], a C++ code for calculation of bond order parameters, ``pyscal`` is a python module designed for the computation of bond orientational order parameters during post-processing of simulation data. In addition to Steinhardt parameters and their averaged version calculated by ``BondOrderAnalysis``, ``pyscal`` brings together the various approaches described above in a single python module. While python offers the advantages of extensibility and flexibility, the module ensures the speed and efficiency of the calculations by using a core code written in C++, which is then ported to python using ``pybind11`` [@Jakob:2016]. The ``pyscal``module includes the following functionality-   

 * calculation of Steinhardt bond orientational order parameters and their averaged version.
 * weighted bond order parameters using face area of voronoi polyhedra, the calculation of which is  enabled using Voro++ code [@Rycroft:2009] integrated into ``pyscal``.
 * distinction of atoms as solid or liquid based on Frenkel and Auer approach [@Auer:2005].
 * improvement of local neighborhood description using an adaptive cutoff distance.
 * clustering algorithm of particles based on a user defined property.
 * inbuilt functions for other structural features like radial distribution function, voronoi volume of individual particles and coordination number.
 * calculation of other Voronoi tessellation based parameters such as number of vertices and face area.

``pyscal`` uses a list of particle positions and simulation box vectors as input. It can also read in output files containing atomistic simulation data in the LAMMPS [@Plimpton:1995] [dump format](https://lammps.sandia.gov/doc/dump.html) and [POSCAR](https://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html) format used by _ab initio_ simulation package [VASP](https://www.vasp.at/). Orthogonal and triclinic simulation cells can be employed in the module. ``pyscal`` also offers a set of supporting tools including those for reading, writing and splitting of atomic trajectory files and the creation of perfect crystal structures. Overall, ``pyscal`` also acts as a tool for calculation of quantities based on atomic positions, and can be easily extended using either C++ for computationally intensive operations or python for other tasks that wrap over the existing features. ``pyscal`` includes a documentation and various usage examples, available in the [pyscal website](https://pyscal.readthedocs.io/en/latest/). Other examples, including installation of the package and the calculation of bond orientational order parameters, are also available. Alternatively, ``pyscal`` also provides a [binder environment](https://mybinder.org/v2/gh/srmnitc/pyscal/master?filepath=examples%2F) to use the example cases before installation.



# Acknowledgements
S.M. acknowledges a scholarship from the International Max Planck Research School for Interface Controlled Materials for Energy Conversion. G.D.L acknowledges the support by the Mexican National Council for Science and Technology (CONACYT) through project 232090 and by the German Research Foundation (DFG) through project RO 3073/6-1.

# References