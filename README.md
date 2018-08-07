

```python
import steinhardt as st
```


```python
import numpy as np
```


```python
%load_ext memory_profiler
```


```python
help(np.arange)
```

    Help on built-in function arange in module numpy.core.multiarray:
    
    arange(...)
        arange([start,] stop[, step,], dtype=None)
        
        Return evenly spaced values within a given interval.
        
        Values are generated within the half-open interval ``[start, stop)``
        (in other words, the interval including `start` but excluding `stop`).
        For integer arguments the function is equivalent to the Python built-in
        `range <http://docs.python.org/lib/built-in-funcs.html>`_ function,
        but returns an ndarray rather than a list.
        
        When using a non-integer step, such as 0.1, the results will often not
        be consistent.  It is better to use ``linspace`` for these cases.
        
        Parameters
        ----------
        start : number, optional
            Start of interval.  The interval includes this value.  The default
            start value is 0.
        stop : number
            End of interval.  The interval does not include this value, except
            in some cases where `step` is not an integer and floating point
            round-off affects the length of `out`.
        step : number, optional
            Spacing between values.  For any output `out`, this is the distance
            between two adjacent values, ``out[i+1] - out[i]``.  The default
            step size is 1.  If `step` is specified, `start` must also be given.
        dtype : dtype
            The type of the output array.  If `dtype` is not given, infer the data
            type from the other input arguments.
        
        Returns
        -------
        arange : ndarray
            Array of evenly spaced values.
        
            For floating point arguments, the length of the result is
            ``ceil((stop - start)/step)``.  Because of floating point overflow,
            this rule may result in the last element of `out` being greater
            than `stop`.
        
        See Also
        --------
        linspace : Evenly spaced numbers with careful handling of endpoints.
        ogrid: Arrays of evenly spaced numbers in N-dimensions.
        mgrid: Grid-shaped arrays of evenly spaced numbers in N-dimensions.
        
        Examples
        --------
        >>> np.arange(3)
        array([0, 1, 2])
        >>> np.arange(3.0)
        array([ 0.,  1.,  2.])
        >>> np.arange(3,7)
        array([3, 4, 5, 6])
        >>> np.arange(3,7,2)
        array([3, 5])
    



```python
%%timeit -o -n 5
%memit
! python calculate_opmofast.py traj.light
```

    peak memory: 31.71 MiB, increment: 0.00 MiB
    peak memory: 31.71 MiB, increment: 0.00 MiB
    peak memory: 31.71 MiB, increment: 0.00 MiB
    peak memory: 31.71 MiB, increment: 0.00 MiB
    peak memory: 31.71 MiB, increment: 0.00 MiB
    peak memory: 31.71 MiB, increment: 0.00 MiB
    peak memory: 31.71 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    peak memory: 31.72 MiB, increment: 0.00 MiB
    5 loops, best of 3: 691 ms per loop





    <TimeitResult : 5 loops, best of 3: 691 ms per loop>




```python
%%timeit -o -n 5
%memit
! python calculate_opmofast2.py traj.light
```

    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    peak memory: 31.74 MiB, increment: 0.00 MiB
    5 loops, best of 3: 798 ms per loop





    <TimeitResult : 5 loops, best of 3: 798 ms per loop>




```python
%%timeit -o -n 5
%memit
! python calculate_opmolib.py traj.light
```

    peak memory: 31.75 MiB, increment: 0.01 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    peak memory: 31.75 MiB, increment: 0.00 MiB
    5 loops, best of 3: 589 ms per loop





    <TimeitResult : 5 loops, best of 3: 589 ms per loop>




```python
%%timeit -o -n 5
%memit
! python calculate_opmolib2.py traj.light
```

    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    peak memory: 31.77 MiB, increment: 0.00 MiB
    5 loops, best of 3: 569 ms per loop





    <TimeitResult : 5 loops, best of 3: 569 ms per loop>




```python
%%timeit -o -n 5
%memit
! python calculate_opmolib3.py traj.light
```

    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    peak memory: 31.78 MiB, increment: 0.00 MiB
    5 loops, best of 3: 601 ms per loop





    <TimeitResult : 5 loops, best of 3: 601 ms per loop>




```python
sys =st.System()
```


```python
sys.set_inputfile('tests/conf.dump')
```


```python
sys.set_neighbordistance(3.63)
sys.set_nucsize_parameters(7,0.5,0.5)
```


```python
nuc = sys.calculate_nucsize()
nuc
```




    63L




```python
a = st.Atom() 
b = st.Atom()
```


```python
a = sys.gatom(1)
b = sys.gatom(3)
print a.gneighbors()
print b.gx()
```

    [0L, 11L, 14L, 15L, 60L, 62L, 71L, 73L, 74L, 75L, 84L, 421L, 422L]
    [-0.10301, -6.35752, -6.44787]



```python
sys.get_abs_distance(a,b)

```




    3.71396015167099




```python

```
