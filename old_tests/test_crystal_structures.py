
import pyscal.crystal_structures as pcs
import os

def test_create_structure():
    structure = ['bcc', 'fcc', 'hcp', 'diamond', 'a15', 'l12', 'b2']
    natoms = [2, 4, 4, 8, 8, 4, 2]

    for count, struct in enumerate(structure):
        #dumpfile = os.path.join(os.getcwd(), "test.dat")
        atoms, boxdims = pcs.make_crystal(struct)
        assert len(atoms) == natoms[count]
