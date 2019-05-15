import pytest
import os,sys,inspect
import numpy as np
import pybop.core as pc
import pybop.crystal_structures as pcs
import pybop.traj_process as ptp

def test_create_multislice_dump():
    """
    Create a multitest dump file and test it
    """
    atoms, boxdims = pcs.make_crystal('bcc', repetitions=[6,6,6])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)
    ptp.write_structure(sys, "tests/bcc1.dump")

    atoms2, boxdims2 = pcs.make_crystal('bcc', repetitions=[6,6,6])
    #modify the coordinates of one atom
    x  = atoms2[0].get_x()
    x[0] += 0.01
    atoms2[0].set_x(x)
    #write it out
    sys2 = pc.System()
    sys2.assign_atoms(atoms2, boxdims2)
    ptp.write_structure(sys2, "tests/bcc2.dump")

    #now merge the two dump files
    os.system("cat tests/bcc1.dump tests/bcc2.dump > tests/bcc3.dat")
    os.remove("tests/bcc1.dump")
    os.remove("tests/bcc2.dump")

    #now this file should have info of both - read it in
    sys3 = pc.System()
    sys3.read_inputfile("tests/bcc3.dat", frame=1)
    atoms = sys3.get_allatoms()
    assert atoms[0].get_x() == [0.01,0,0]    

    #now this file should have info of both - read it in
    sys4 = pc.System()
    sys4.read_inputfile("tests/bcc3.dat", frame=0)
    atoms = sys4.get_allatoms()
    assert atoms[0].get_x() == [0.0,0,0]

    #now cleanup
    os.remove("tests/bcc3.dat")    