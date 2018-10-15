
# Steinhardt Tools

Steinhardt tools is a python library written in C++ to calculate and analyse Steinhardts parameters[1]. It provides the necessary tools to easily calculate and extend upon the Steinhardt parameters calculations.

### Installation

The module can be installed by using `pip install ./SteinhardtTools` after cloning the repository.

The module can be imported as-


```python
import steinhardt as st
```

### Documentation

This module has two main member classes - `Atom` and `System`. `Atom` class is used for storing the properties of a single atom whereas `System` has all the properties of the system. The functions in each modules can be accessed by
`dir(st.Atom)` or `dir(st.System)`. The documentation for a member function can be accessed by `help(function_name)`. 

##### Example


```python
help(st.System.calculate_nucsize)
```

    Help on method calculate_nucsize in module steinhardt:
    
    calculate_nucsize(...) unbound steinhardt.System method
        calculate_nucsize(self: steinhardt.System) -> int
        
        
        Calculate the size of the largest cluster in the given system. Calculation
        the size of the largest cluster needs various prerequisites that can be set
        by the functions set_neighbordistance and set_nucsize_parameters. 
        For a detailed description of how the calculation works see-
        Diaz Leines et al, JCP 146(2017). http://doi.org/10.1063/1.4980082
        
        Parameters
        ----------
        None
        
        Returns
        -------
        cluster size : int
            size of the largest cluster in number of atoms
    


### Examples

The first step involved is creating a system.


```python
sys = st.System()
```

The system created now is empty and we need to fill it with data about atomic positions and simulation box. A sample [LAMMPS](https://lammps.sandia.gov/) dump file is given in `tests/conf.dump`. This file can be specified as the input file by 


```python
sys.set_inputfile('tests/conf.dump')
```

After setting the name of the input file, it needs to be read in, which is done using-


```python
sys.read_particle_file()
```

*Accessing atom information*  
All access functions start with the letter 'g' for get and all the functions used to set values start with the letter 's' for set. Atom class has access functions for the position and id. However, we need to get an atom from the system. This can be done as follows- 



```python
atom1 = sys.gatom(2)
type(atom1)
```




    steinhardt.Atom



The `gatom` function return the atom from the list index that is given as the argument. Now we have the atom, and we can check its coordinates and id. There are more features of the atom that can be checked, which will be discussed further in the examples.


```python
coordinates = atom1.gx()
coordinates
```




    [-0.947172, 10.3193, 11.3033]



Now the coordinates of the atom can be modified, for example to check how the value of $\bar{q}$ parameters changes and it can be put back in its original location.


```python
coordinates[0] += 2.00
atom1.sx(coordinates)
sys.satom(atom1)
```

Now it is possible to extract another atom and calculate the distance between two atoms.


```python
atom2 = sys.gatom(3)
sys.get_abs_distance(atom1,atom2)
```




    2.7748271598685195



*Calculating the size of the largest cluster*  
Next we will take a look at the calculation of the largest cluster size. A detailed description of each parameter and the method is available in [2], hence, we will skip the explanation here. There are certain parameters we need to set - namely cutoff distance to find neighbors, minimum number of solid neighbors, threshold and average threshold which can be set by-



```python
sys.set_neighbordistance(3.63)
sys.set_nucsize_parameters(7,0.5,0.5)
```

Now the largest cluster size can be calculated.


```python
nuc = sys.calculate_nucsize()
nuc
```




    63L



There were a lot of properties that can be accessed. We can, for example, find the neighbors of an atom.


```python
atom1 = sys.gatom(2)
atom1.gneighbors()
```




    [3L, 62L, 64L, 421L, 422L, 424L, 473L, 487L, 488L, 490L, 492L, 493L]



It is possible to check if an atom is solid or not by the `gissolid` function which return a value of 1 if it is solid or 0 if it is not. 


```python
atom1.gissolid()
```




    0L



This makes it very easy to calculate further properties, for example, the centre of mass of the solid cluster(assuming mass as unity). Start by getting all the atoms-


```python
atoms = sys.gallatoms()
```

Now we can loop over the solid ones


```python
#define sum variables
count=0
com = [0,0,0]

#loop over all atoms
for atom in atoms:
    #check if it is a solid
    if atom.gissolid():
        #then add the coordinates
        coord = atom.gx()
        com[0]+=coord[0]
        com[1]+=coord[1]
        com[2]+=coord[2]
        count+=1
com[0]/=float(count)
com[1]/=float(count)
com[2]/=float(count)
com
```




    [-3.918379215384615, 1.2204852923076925, 1.225490959076923]





### References

[1] Steinhardt, P. J., Nelson, D. R., & Ronchetti, M. (1983). Bond-orientational order in liquids and glasses. Physical Review B, 28(2), 784–805. http://doi.org/10.1103/PhysRevB.28.784  
[2] Díaz Leines, G., Drautz, R., & Rogal, J. (2017). Atomistic insight into the non-classical nucleation mechanism during solidification in Ni. Journal of Chemical Physics, 146(15). http://doi.org/10.1063/1.4980082
