import pytest
import os
import numpy as np
import pyscal.core as pc
#from pyscal.trajectory import Trajectory, hdf_to_dump
from pyscal.trajectory import Trajectory

def test_traj():
	traj = Trajectory("examples/traj.light")
	assert traj.nblocks == 10

	assert len(traj.get_block(0)) == 509

	traj.load(0)
	data = traj.data[0]
	assert data["box"][0][0] ==  -7.34762
	assert data["atoms"]["x"][0] == -4.72745
	traj.unload(0)

def test_timeslice():
	traj = Trajectory("examples/traj.light")
	assert traj.nblocks == 10
	sys = traj[0].to_system()
	assert sys[0].box[0][0] == 18.21922

	aseobj = traj[0].to_ase(species=["Au"])
	assert aseobj[0].positions[0][0] == -4.72745

	od = traj[0].to_dict()
	assert od[0]["box"][0][0] ==  -7.34762

	traj[0].to_file("test.out")
	assert os.path.exists("test.out") == True

	#traj[0].to_hdf("test.hdf")
	#assert os.path.exists("test.hdf") == True

	hdf_to_dump("test.hdf", "test.dat")
	assert os.path.exists("test.dat") == True	