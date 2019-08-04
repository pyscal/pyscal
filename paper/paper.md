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

Atomistic simulation methods provide an insight into materials and their properties from an atomistic perspective. These methods often produce large amounts of data in terms of atomic positions over time, the analysis of which is crucial to study underlying mechanisms. Many computational analysis methods have been developed to determine the local environment of a particle, which are used for various purposes such as the distinction of a particle as solid or liquid, to identify the crystal structure they long to, to identify defects in a crystalline state etc. to name a few. (Stukowski). 

Bond orientational order parameters introduced by Steinhardt(steinhardt), have been used extensively(citations). These parameters, which are rotationally invariant are defined by,


$$ \bar{q}_l (i) =  \Big(  \frac{2\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\pmb{r}_{ij}) \Big|^2 \Big )^{\frac{1}{2}} $$

where $Y_{lm}$ are the spherical harmonics and $N(i)$ is the number of nearest neighbours of particle $i$. $\pmb{r}_{ij}$ is the displacement vector connecting particles $i$ and $j$ and $l$ and $m$ are both intergers with $m \in [-l,+l]$. $q_4$ and $q_6$ are often used for distinction of crystal structures. Furthermore, the averaged version of the order parameters was introduced to improve the resolution by averaging over the bond order parameters of nearest neighbour particles.

At high temperatures, atomic positions   





# Acknowledgements

# References