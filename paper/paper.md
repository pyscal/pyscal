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

Atomistic simulation methods are widely used in the study of condensed matter systems, often producing large amounts of data in terms of atomic positions over time. This data is then analysed to determine the local environment of a particle, which is used for purposes such as the distinction of a particle as solid or liquid, to identify the crystal structure, and to identify defects in a crystalline state. Many analysis methods have been developed for this purpose[@Stukowski:2012], of which the bond orientational order parameters introduced by Steinhardt and coworkers [@Steinhardt:1983], have been used extensively[@Mickel:2013 and references within]. These parameters, which are rotationally and translationally invariant are defined by,


$$ q_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l | q_{lm}(i) |^2 \Big )^{\frac{1}{2}} $$ 

where, 

$$ q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\pmb{r}_{ij}) $$ 

in which $Y_{lm}$ are the spherical harmonics and $N(i)$ is the number of neighbours of particle $i$, $\pmb{r}_{ij}$ is the vector connecting particles $i$ and $j$, and $l$ and $m$ are both intergers with $m \in [-l,+l]$.  

However, at high temperatures, atomic positions are subject to thermal vibrations which leads to a broad distribution of $q_l$ values, and the overlap between these distributions may make the resolution of crystal structures harder. The averaged version of the order parameters ($\bar{q}_l$), was introduced [@Lechner:2008] in order to counter this problem. $\bar{q}_l$ is given by,

$$ \bar{q}_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{\tilde{N}(i)} \sum_{k=0}^{\tilde{N}(i)} q_{lm}(k) \Big|^2 \Big )^{\frac{1}{2}} $$ 

where the sum from $k=0$ to $\tilde{N}(i)$ is over all the neighbors and the particle itself. ($\bar{q}_l$) was shown to decrease the overlap between the $q_l$ distributions, leading to an improvement in the resolution of crystal structures. Another problem at high temperatures is the use of a fixed cutoff radius which leads to ambiguity in the definition of the local neighbourhood of a particle. A parameter free definition of the local environment can be achieved by using a distinct cutoff radius for each particle using methods such as the SANN algorithm [@VanMeel:2012] or an adaptive cutoff method similar to to the approach followed in adaptive common neighbour analysis [@Stukowski:2012]. Additionally, Voronoi tessellation can be employed for a parameter free definition of the local neighbourhood. It was shown [@Mickel:2013] that the distinction of crystal structures can be improved by weighting the contribution of each neighbour to $q_{lm}$ given by,         

$$ q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} \frac{A_{ij}}{A} Y_{lm}(\pmb{r}_{ij}) $$

where $A_{ij}$ is area of the Voronoi facet between atoms $i$ and $j$, and $A$ is the total surface area of the Voronoi polyhedra. Futher approaches include the use of higher powers, $A_{ij}^\alpha$, with $\alpha=2,3$ to weight the contribution to $q_{lm}$[@Haeberle:2019].

Originally inspired by [``BondOrderAnalysis``](https://homepage.univie.ac.at/wolfgang.lechner/bondorderparameter.html) [@Lechner:2010], a C++ code for calculation of bond order parameters, ``pybop`` is a python module designed for the computation of bond orientational order parameters during post-processing of simulation data. ``pybop`` brings together the various approaches described above in a single python module while ensuring the speed and efficiency of calculations as the core code for ``pybop`` is written in C++, which is then ported to python using ``pybind11`` [@Jakob:2016]. The module includes the following functionality-   

 * calculation of bond order parameters ($q_{l}$) and the averaged version ($\bar{q}_{l}$).
 * weighted bond order parameters using face area of voronoi polyhedra, the calculation of which is  enabled using Voro++ code [@Rycroft:2009] integrated into ``pybop``.
 * distinction of atoms as solid or liquid based on $q_6$ [@Auer:2005].
 * improvement of local neighbourhood description using an adaptive cutoff.
 * clustering algorithm to find clusters of particles based on a user defined property.
 * inbuilt functions for other structural features like radial distribution function, voronoi volume of individual particles and coordination number.  

``pybop`` uses a list of particle positions and simulation box vectors as input. It can also read in output files containing atomistic simulation data in the LAMMPS [@Plimpton:1995] [dump format](https://lammps.sandia.gov/doc/dump.html) and [POSCAR](https://cms.mpi.univie.ac.at/vasp/vasp/POSCAR_file.html) format used by _ab initio_ simulation package [VASP](https://www.vasp.at/), including orthogonal as well as triclinic simulation cells.``pybop`` also offers a set of support tools including those for reading, writing and splitting of atomic trajectory files and creation of perfect crystal structures. Various features of ``pybop`` including documentation and example usage is available in the [pybop website](https://pybop.readthedocs.io/en/latest/). Examples including installation of the package to calculation of bond orientational order parameters are provided. Alternatively, ``pybop`` provides a [binder environment](https://mybinder.org/v2/gh/srmnitc/pybop/master?filepath=examples%2F) to try the example cases before installation.



# Acknowledgements

# References