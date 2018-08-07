import pytest
import os,sys,inspect
#we need import from previous folder - just to organise things
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
#at this point we should have the module
 
import steinhardt as st
#this is the testing file for the steinhardt module- it will create a test bed for all the functions
#this has to be cleared before each commit to make sure everything works
#the test module requires the file conf.dump which is available in this folder.

#just a test to see if everything works
def test_this_file():
        assert 1==1

#tests for atom class. create an atom class and check everything

def test_atom_class():

        atom = st.Atom()

        #check the setting of positions and id
        atom.sx([1,2,3])
        atom.sid(1)

        assert atom.gx()==[1,2,3]
        assert atom.gid()==1

def test_system_readfile():
        #now set up a system module
        sys = st.System()
        sys.set_inputfile('tests/conf.dump')
        sys.read_particle_file()

        #fetch an atom
        a = st.Atom()
        a = sys.gatom(0)
        assert a.gid()==3

def test_get_abs_distance():

        sys = st.System()
        sys.set_inputfile('tests/conf.dump')
        sys.read_particle_file()

        a =sys.gatom(1)
        b =sys.gatom(3)

        assert sys.get_abs_distance(a,b)==3.71396015167099


def test_get_all_neighbors():

        sys = st.System()
        sys.set_inputfile('tests/conf.dump')
        sys.read_particle_file()

        sys.set_neighbordistance(3.63)
        sys.get_all_neighbors()

        a =sys.gatom(1)

        assert a.gneighbors()==[0,11,14,15,60,62,71,73,74,75,84,421,422]


# now a compund test for nucsize
def test_calculate_nucsize():
        sys = st.System()
        sys.set_inputfile('tests/conf.dump')
        sys.read_particle_file()

        sys.set_neighbordistance(3.63)
        sys.set_nucsize_parameters(7,0.5,0.5)
        nuc = sys.calculate_nucsize()

        assert nuc==63
