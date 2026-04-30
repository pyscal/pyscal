<!--
Draft text for a pinned issue on https://github.com/pyscal/pyscal
Title suggestion:
    📢 pyscal v4 is out — this repository (v2) is no longer maintained
Recommended labels: `announcement`, `migration`
After opening: Issue → ⋯ menu → "Pin issue" so it sits at the top.
-->

# pyscal v4 is out — this repository (v2) is no longer maintained

Hi everyone — a quick note for anyone landing here.

This repository (`pyscal/pyscal`) is **pyscal v2**, the original release.
It is **no longer actively developed**: no new features, no bug fixes,
no python-version updates. The repository is kept online so that links
from papers, tutorials, and the JOSS citation continue to resolve.

The current version of pyscal is **v4**, developed in a separate repo:
👉 **<https://github.com/pyscal/pyscal3>**

## Installing v4

```bash
pip install pyscal3
# or
conda install -c conda-forge pyscal3
```

After install, both work and refer to the same library:

```python
import pyscal     # preferred going forward
# or
import pyscal3
```

(The PyPI name is `pyscal3` because the name `pyscal` on PyPI refers to
this v2 package, which we are leaving in place for backward compatibility.)

## What changed

v4 is a ground-up rewrite:

- **ASE `Atoms` is the data structure.** `pyscal.System` is gone. Every
  function takes an `ase.Atoms` object and writes results back into
  `atoms.arrays["pyscal_*"]` and `atoms.info["pyscal_*"]`.
- **Functional API.** Every descriptor is a top-level function — no
  classes to subclass, no method-chaining order to remember.
- **More descriptors:** Wigner $W_l$, Minkowski metrics, Ackland–Jones
  classifier, three coordination-number variants, angular & bond-length
  distributions, atomic strain / $D^2_{\min}$ / slip vector, Wigner–Seitz
  defect analysis, ACE descriptors, and more.

## Migration sketch

```python
# v2 (this repo)
import pyscal.core as pc
sys = pc.System()
sys.read_inputfile('conf.dump')
sys.find_neighbors(method='cutoff', cutoff=3)
sys.calculate_q([4, 6])
atoms = sys.atoms
q6 = [a.get_q(6) for a in atoms]
```

```python
# v4
import pyscal
from ase.io import read
atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='cutoff', cutoff=3)
pyscal.steinhardt_parameter(atoms, l=[4, 6])
q6 = atoms.arrays['pyscal_q6']
```

Full migration notes and the new API reference are at <https://pyscal.org>.

## What about existing v2 installs?

`pip install pyscal` will continue to install this v2 package — nothing
breaks for existing users. We just won't be releasing new v2 versions.

## Reporting issues / asking questions

For anything related to v4, please open issues at
<https://github.com/pyscal/pyscal3/issues>. New issues filed here on the
v2 repo will not be triaged.

Thanks to everyone who has used and cited pyscal over the years 🙏
