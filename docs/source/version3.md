# Why version 3?

pyscal v3 is a new version with mostly updated codebase and breaking changes. Anybody who has working pyscal code will need to update it to get it working with this new version. Therefore, it is necessary to discuss why this new version was needed and the benefits of updating.

## Version 3 is much faster

In the plot below, the time needed to calculate neighbors with the 'cutoff' method for systems with varying number of atoms with versions 2.10.15 and 3.0 is shown.

<img src="img_time_neighbor.png"  width="60%">

v3 is faster for all system sizes. At a system size of about 50,000 atoms, v3 is about 4x faster.

## Version 3 uses less memory

A major issue with pyscal v2.x series was that it not useful for large system sizes due to the large amount of memory needed. In the plot below, the memory usage of both versions for the same calculation above is shown.

<img src="img_time_memory.png"  width="60%">

v3 uses less memory, for a system size of 50,000 atoms, v3 uses 14x less memory. A more interesting feature is the slope of the data, or how much the memory scales with the system size. For v3 it is only 0.008, while for v2 it is .12! For a system of 1 million atoms, v2 would use 117 GB of memory while v3 would need only 8 GB, making larger calculations accessible (these numbers will be updated after real use-case tests).

## What are reasons for these benefits?

- The older C++ atoms class is deprecated. Instead, it is store as python dictionary. Therefore the copying between python and C++ sides is avoided.

- The atoms python dictionary is directly exposed to the C++ side. The dictionary is passed by reference, which allows in-place modification directly.

## What are the other feature updates?

The new version includes a number of new features and quality of life improvements. Please check the examples for details.