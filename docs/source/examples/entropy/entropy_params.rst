Entropy parameters
------------------

In this example, the entropy parameters are calculated and used for
distinction of solid and liquid. For a description of entropy
parameters, see `here <add%20link>`__.

.. code:: python

    import pyscal.core as pc
    import matplotlib.pyplot as plt
    import numpy as np

We have two test configurations for Al at 900 K, one is fcc structured
and the other one is in liquid state. We calculate the entropy
parameters for each of these configurations. First we start by reading
in the fcc configuration. For entropy parameters, the values of the
integration limit :math:`r_m` and the value of cutoff :math:`r_a` are
chosen as 1.4 and 0.9, based on the `original
publication <https://aip.scitation.org/doi/10.1063/1.4998408>`__.

.. code:: python

    sys = pc.System()
    sys.read_inputfile("../tests/conf.fcc.Al.dump")
    sys.find_neighbors(method="cutoff", cutoff=0)

The values of :math:`r_m` and :math:`r_a` are in units of lattice
constant, so we need to calculate the lattice constant first. Since is a
cubic box, we can do this by,

.. code:: python

    lat = (sys.box[0][1]-sys.box[0][0])/5

Now we calculate the entropy parameter and its averaged version.

.. code:: python

    sys.calculate_entropy(1.4*lat, ra=0.9*lat, averaged=True)

The calculated values are stored for each atom. This can be accessed as
follows,

.. code:: python

    atoms = sys.atoms
    solid_entropy = [atom.entropy for atom in atoms]
    solid_avg_entropy = [atom.avg_entropy for atom in atoms]

Now we can quickly repeat the calculation for the liquid structure.

.. code:: python

    sys = pc.System()
    sys.read_inputfile("../tests/conf.lqd.Al.dump")
    sys.find_neighbors(method="cutoff", cutoff=0)
    lat = (sys.box[0][1]-sys.box[0][0])/5
    sys.calculate_entropy(1.4*lat, ra=0.9*lat, averaged=True)
    atoms = sys.atoms
    liquid_entropy = [atom.entropy for atom in atoms]
    liquid_avg_entropy = [atom.avg_entropy for atom in atoms]

Finally we can plot the results

.. code:: python

    xmin = -3.55
    xmax = -2.9
    bins = np.arange(xmin, xmax, 0.01)
    x = plt.hist(solid_entropy, bins=bins, density=True, alpha=0.5, color="#EF9A9A")
    x = plt.hist(solid_avg_entropy, bins=bins, density=True, alpha=0.5, color="#B71C1C")
    x = plt.hist(liquid_entropy, bins=bins, density=True, alpha=0.5, color="#90CAF9")
    x = plt.hist(liquid_avg_entropy, bins=bins, density=True, alpha=0.5, color="#0D47A1")
    plt.xlabel(r"$s_s^i$")



.. image:: fig_1.png


The distributions of :math:`s_s^i` given in light red and light blue are
fairly distinct but show some overlap. The averaged entropy parameter,
:math:`\bar{s}_s^i` show distinct peaks which can distinguish solid and
liquid very well.
