import pytest
import os,sys,inspect
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs
import pyscal.traj_process as ptp

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
    atoms = sys3.get_atoms()
    assert atoms[0].pos == [0.01,0,0]

    #now this file should have info of both - read it in
    sys4 = pc.System()
    sys4.read_inputfile("tests/bcc3.dat", frame=0)
    atoms = sys4.get_atoms()
    assert atoms[0].pos == [0.0,0,0]

    #now cleanup
    os.remove("tests/bcc3.dat")

def test_customvals_dump():
    """
    Test writing customvals
    """
    atoms, boxdims = pcs.make_crystal('bcc', repetitions=[1,1,1])
    sys = pc.System()
    sys.assign_atoms(atoms, boxdims)

    #test for multiple customvals
    customks = 'one'
    customvs = [1,1]
    ptp.write_structure(sys, "tests/bcc4.dump", customkey=customks, customvals=customvs)

    #now read this file
    lines = []
    for line in open("tests/bcc4.dump", 'r'):
        lines.append(line)

    #now check the atoms
    last1line = lines[-1].strip().split()
    last2line = lines[-2].strip().split()
    last3line = lines[-3].strip().split()

    #now verify
    assert last1line[-1] == '1'
    assert last2line[-1] == '1'
    assert last3line[-1] == 'one'

    #clean up
    if os.path.exists("tests/bcc4.dat"):
        os.remove("tests/bcc4.dat")

    #test for multiple customvals
    customks = ['one', 'two']
    customvs = [[1,1], [2,2]]
    ptp.write_structure(sys, "tests/bcc4.dump", customkey=customks, customvals=customvs)

    #now read this file
    lines = []
    for line in open("tests/bcc4.dump", 'r'):
        lines.append(line)

    #now check the atoms
    last1line = lines[-1].strip().split()
    last2line = lines[-2].strip().split()
    last3line = lines[-3].strip().split()

    #now verify
    assert last1line[-1] == '2'
    assert last1line[-2] == '1'
    assert last2line[-1] == '2'
    assert last2line[-2] == '1'
    assert last3line[-1] == 'two'
    assert last3line[-2] == 'one'

    #clean up
    if os.path.exists("tests/bcc4.dat"):
        os.remove("tests/bcc4.dat")
