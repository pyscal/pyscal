
Angular criteria for identification of diamond structure
--------------------------------------------------------

Angular parameter introduced by Uttormark et al [1]_ is used to measure the tetrahedrality
of local atomic structure. An atom
belonging to diamond structure has four nearest neighbors which gives
rise to six three body angles around the atom. The angular parameter
:math:`A` is then defined as,

:math:`A = \sum_{i=1}^6 (\cos(\theta_i)+\frac{1}{3})^2`

An atom belonging to diamond structure would show the value of angular
params close to 0. Angular parameter can be calculated in pyscal using the following
method -

.. code:: python

    import pyscal.core as pc
    sys = pc.System()
    sys.read_inputfile('conf.dump')
    sys.find_neighbors(method='cutoff', cutoff='adaptive')
    sys.calculate_angularcriteria()


The calculated angular criteria value can be accessed for each atom using :attr:`~pyscal.catom.Atom.angular`.

.. [1] Uttormark, MJ, Thompson, MO, Clancy, P, Phys. Rev. B 47, 1993
