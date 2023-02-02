import pytest
import os
import numpy as np
import pyscal.core as pc
from pyscal.crystal_structures import Structure
import pyscal.traj_process as ptp

def test_create_multislice_dump():
    """
    Create a multitest dump file and test it
    """
    sys = Structure().lattice.bcc(repetitions = [6, 6, 6])
    ptp.write_file(sys, "tests/bcc1.dump")

    atoms2, boxdims2 = pcs.make_crystal('bcc', repetitions=[6,6,6])
    #modify the coordinates of one atom
    x  = atoms2["positions"][0]
    x[0] += 0.01
    atoms2["positions"][0] = x
    assert len(atoms2["positions"]) == 432
    #write it out
    sys2 = pc.System()
    sys2.box = boxdims2
    sys2.atoms = atoms2

    ptp.write_file(sys2, "tests/bcc2.dump")

    #now cleanup
    if os.path.exists("tests/bcc1.dat"):
        os.remove("tests/bcc1.dat")
    if os.path.exists("tests/bcc2.dat"):
        os.remove("tests/bcc2.dat")

def test_customvals_dump():
    """
    Test writing customvals
    """
    sys = Structure().lattice.bcc(repetitions = [1, 1, 1])


    #test for multiple customvals
    customks = ['one']
    customvs = [[1],[1]]
    ptp.write_file(sys, "tests/bcc4.dump", customkeys=customks, customvals=customvs)

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
    ptp.write_file(sys, "tests/bcc4.dump", customkeys=customks, customvals=customvs)

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
    assert last1line[-2] == '2'
    assert last2line[-1] == '1'
    assert last2line[-2] == '1'
    assert last3line[-1] == 'two'
    assert last3line[-2] == 'one'

    #clean up
    if os.path.exists("tests/bcc4.dat"):
        os.remove("tests/bcc4.dat")
