---
title: 'pyscal: A python module for structural analysis of atomic environments'
tags:
  - Python
  - materials science
  - bond orientational order parameters
authors:
  - name: Sarath Menon
    orcid: 0000-0002-6776-1213
    affiliation: 1
  - name: Grisell Díaz Leines
    orcid: 0000-0002-3048-5315
    affiliation: 1
  - name: Jutta Rogal
    orcid: 0000-0002-6268-380X
    affiliation: 1
affiliations:
 - name: Interdisciplinary Centre for Advanced Materials Simulation, Ruhr-Universität Bochum, Universitätsstr. 150, 44801 Bochum, Germany.
   index: 1
date: 05 August 2019
bibliography: paper.bib
---


# Summary

The structural characterization of local atomic environments is essential in our understanding and design of materials with customized properties.  Atomistic simulations have become a powerful tool to provide insight into the structural patterns and mechanisms that signify transformations between different crystalline and liquid phases, or the formation and dynamics of point and extended defects.  The development and availability of methods that provide a reliable and accurate structural analysis of atomistic simulation data constitute an indispensable task that continues to be an active field of research in materials science.

Several techniques have been proposed over the years to identify local structural motifs in crystalline materials [@Stukowski:2012].
To distinguish solid- and liquid-like environments, Auer and Frenkel [@Auer:2005] introduced a criterion that measures the structural correlations between a particle and its neighbors using the bond orientational order parameters suggested by Steinhardt et al. [@Steinhardt:1983].
The Steinhardt parameters are based on the spherical harmonics in the local neighborhood of a particle within a cutoff distance, and are rotationally and translationally invariant. They have been used extensively for the identification of crystal structures [@Mickel:2013 and references within].
At finite temperatures the Steinhardt parameters for a given crystal structure exhibit a certain distribution of values instead of a single one due to thermal vibrations of the atoms.  This can lead to an overlap in the distributions for different phases which makes an accurate identification difficult or even impossible.
An averaged version of the Steinhardt bond order parameters introduced by Lechner and Dellago [@Lechner:2008] significantly reduces this overlap and improves the characterization of common crystal structures.
An additional complication at finite temperatures is the use of a fixed cutoff distance which might cause ambiguity in the definition of the local neighborhood of a particle. An adaptive definition of the local environment can be achieved by  methods such as the solid angle based nearest neighbor (SANN) algorithm [@VanMeel:2012] or the adaptive common neighbor analysis (a-CNA) [@Stukowski:2012].
Furthermore, Voronoi tessellation can be employed for a parameter-free definition of the local neighborhood of a particle. Mickel et al. showed [@Mickel:2013] that the identification of crystal structures can be improved by weighting the contribution of each neighbor to the Steinhardt parameters by the area of the Voronoi facet shared between the central atom and its corresponding neighbor.
Higher exponents of the facet area to weight the bond orientational order parameters can further improve the recognition of common crystal phases [@Haeberle:2019].

Originally inspired by the C++ code [``BondOrderAnalysis``](https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html) [@Lechner:2010] for the calculation of bond order parameters, ``pyscal`` is a Python module designed for the computation of local structural order parameters  during post-processing of atomistic simulation data.
In addition to the calculation of Steinhardt parameters which is offered by tools like [``BondOrderAnalysis``](https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html) and ``pyboo`` [@Leocmach2017], ``pyscal`` brings together the various approaches for structure identification discussed above in a single Python module. While Python offers the advantages of extensibility and flexibility, the module ensures the speed and efficiency of the calculations by using a core code written in C++, which is  ported to Python using ``pybind11`` [@Jakob:2016]. The ``pyscal`` module includes the following functionality:   

 * calculation of Steinhardt bond orientational order parameters and their averaged version.
 * weighted bond order parameters using face areas of Voronoi polyhedra; the calculation of the Voronoi polyhedra is  enabled using the Voro++ code [@Rycroft:2009] integrated into ``pyscal``.
 * distinction of atoms in solid or liquid environments based on the criterion by Auer and Frenkel [@Auer:2005].
 * flexible local neighborhood definition using an adaptive cutoff distance.
 * clustering algorithm for particles based on user defined properties.
 * inbuilt functions for additional structural features including radial distribution functions, Voronoi volumes of individual particles, and coordination numbers.
 * calculation of further parameters based on Voronoi tesselation such as the number of vertices and face areas.

``pyscal`` uses a list of particle positions and simulation cell vectors as input, which can also be read in from a file.  Currently supported file formats include the [dump format](https://lammps.sandia.gov/doc/dump.html) of the molecular dynamics code  LAMMPS [@Plimpton:1995] and the [POSCAR](https://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html) format used by the _ab initio_ simulation package  [VASP](https://www.vasp.at/) [@Kresse:1993; @Kresse:1996; @Furth:1996].  Extensions to other file formats are easily implemented.
The module can handle arbitrary simulation cells.
 ``pyscal`` also offers a set of supporting tools including those for reading, writing, and splitting atomic trajectory files and the setup of simulation cells for perfect crystal structures. In general, the ``pyscal`` module acts as a tool for the calculation of quantities based on atomic positions, and can easily be extended using either C++ for computationally intensive operations or Python for other tasks that wrap over the existing features. ``pyscal`` includes a documentation and various usage examples, available on the [pyscal website](https://pyscal.readthedocs.io/en/latest/). Further examples, including installation of the package and the calculation of bond orientational order parameters, are also available. Alternatively, ``pyscal`` provides a [binder environment](https://mybinder.org/v2/gh/srmnitc/pyscal/master?filepath=examples%2F) to test the example cases before installation.



# Acknowledgements
S.M. acknowledges a scholarship from the International Max Planck Research School for Interface Controlled Materials for Energy Conversion. The authors acknowledge support by the Mexican National Council for Science and Technology (CONACYT) through project 232090 and by the German Research Foundation (DFG) through project 262052203.

# References
