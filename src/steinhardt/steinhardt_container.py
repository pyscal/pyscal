"""
Now a npy file format is used instead of pickled data.

--old info--
Even though super short on time, this module has to be reimagined a bit.
This stems from the fact that pybind11 objects are hard to pickle. So,
we will use alternative classes here - Atomc, and Systemc. These classes
will be used as container ones. Afterwards, we have to find a way to 
exchange information between the C++ classes and these python ones. 

Pickling pybind11 code is important. So is thinkiing about migrating to
Çython. This is a temporary solution. Steinhardt tools will evolve.

"""
import numpy as np
import dask.bag as db
from dask import delayed
import os
import pickle
import shutil

class Atomc(object):
    def __init__(self,idd,x,y,z):
        #basic vals
        self.id = idd
        self.x = x
        self.y = y
        self.z = z
        #other items
        #neighbor related items
        self.neighbors = None
        self.neighborweight = None
        self.neighbordist = None
        self.diff = None
        self.r = None
        self.phi = None
        self.theta = None
        self.n_neighbors = None

        #q related items
        self.q = None
        self.aq = None
        self.realq = None
        self.imgq = None
        self.arealq = None
        self.aimgq = None

        #structure related items
        self.frenkelnumber = None
        self.belongsto = None
        self.issolid = None
        self.structure = None


class Systemc(object):
    """
    System is barebones. the idea is that if we store all the info for
    atoms we dont need to recalculate.
    """
    def __init__(self):
        #main items
        self.atoms = None
        self.natoms = None
        self.inputfile = None
        self.boxdims = None
        
def write_pickle(outfile, systems):
    """
    Write a pickled file with systems

    Parameters
    ----------
    outfile : string 
        name of the output file
    systems : array like - System
        array of Systemc objects
    
    Returns
    -------
    None

    """
    fout = open(outfile,'wb')
    for sys in systems:
        pickle.dump(sys, fout)
    fout.close()

def pickle_systems(infile, natoms, **kwargs):
    """
    --updated function--
    Reads in trajectory and creates Atoms which would be assigned to a class.
    These classes are then stored as a pickled array. New output is a pickled file.

    --old documentation--
    Read in a MD trajectory file and create an atom
    structure from it. For each timeslice, the return array has a slice, which in turn
    contains two sub-slices. The first one is the box dimensions and the second one is
    a list of all atoms in that particular time slices.

    The output data is .npy file. The array is of the dimension [[Box,[Atom1,Atom2..AtomN]],
    [Box,[Atom1,Atom2..AtomN]], .....(number of time slices)]. The output is in delayed
    format (using Dask) to avoid actual reading through of the input data. The x coordinate
    of Atom 0 at time step 1 can be accessed by-
    nsystems[0][1][0].x : The first index is time step, the second to show its atomm and 
    the third finally the index of atom. However this would be delayed object. The actual
    value can be accessed as required by nsystems[0][1][0].x.compute()   

    Parameters
    ----------
    infile : string
        name of the input trajectory file

    natoms : array like
        the number of atoms in each time slice. If the number of atoms in each slice
        is same, a single integer value can be provided. Otherwise a list of entries
        equal to the number of time slices has to be provided.

    **kwargs : 
        delayed (bool) : False, True if a dask delayed system is to be saved. 
        save_file (bool) : True, True if an outfile is to be saved, False otherwise
            save_file is now deprecated and is on by default.
        return_array (bool) : False, True if array is to be returned, False otherwise
            return_array is now deprecated and is off by default
        nslices (int) : 1, number of time slices in the trajectory
        format (string) : lammps, the format of the trajectory file
        outfile (string) : name of the output file
        compressed (bool) : False, if True, compress system and use npz format

    Returns
    -------
    nsystems : array like
        An array of box and atom information for each time slices. Only returned if
        return_array is True.

    """
    #process kwargs
    delay = kwargs.get('delay', False)
    compressed = kwargs.get('compressed', False)
    save_file = kwargs.get('save_file', True)
    return_array = kwargs.get('return_array', False)
    nslices = kwargs.get('nslices', 1)
    format = kwargs.get('format', "lammps")
    outfile = kwargs.get('outfile', os.path.join(os.getcwd(),".".join([infile,"dump"])))

    #read in the dask bag and convert to delayed object
    b = (db.read_text) (infile)
    c = b.to_delayed()

    #initialise systems
    #nsystems = []

    #if natoms is a single value - make it into an array. This allows for providing
    #variable atom numbers in each time slice
    if not isinstance(natoms, list):
        natoms = np.ones(nslices)*natoms

    #move out from pickling and use npy arrays
    outfolder = os.path.join(os.getcwd(),"npydata")
    if os.path.exists(outfolder):
        shutil.rmtree(outfolder)
    os.mkdir(outfolder)

    #this part would be  md code specific
    if format ==  "lammps":
        #fout = open(outfile,'wb')
        for slice in range(nslices):
            fout = ".".join(["snap",str(slice)])
            fout = os.path.join(outfolder, fout)
            nblock = int(natoms[slice]) + 9
            raw = c[0][slice*nblock + 5].strip().split()
            dimxlow = (delayed)(float)(raw[0])
            dimxhigh = (delayed)(float)(raw[1])          
            raw = c[0][slice*nblock + 6].strip().split()
            dimylow = (delayed)(float)(raw[0])
            dimyhigh = (delayed)(float)(raw[1])    
            raw = c[0][slice*nblock + 7].strip().split()
            dimzlow = (delayed)(float)(raw[0])
            dimzhigh = (delayed)(float)(raw[1])
            if not delay:
                dimxlow = dimxlow.compute()
                dimxhigh = dimxhigh.compute()
                dimylow = dimylow.compute()
                dimyhigh = dimyhigh.compute()
                dimzlow = dimzlow.compute()
                dimzhigh = dimzhigh.compute()

            boxdims = [[dimxlow,dimxhigh], [dimylow,dimyhigh], [dimzlow,dimzhigh]]
            atoms = []

            for i in range(9,natoms+9):
                line = c[0][slice*nblock + i].strip().split()
                #print type(slice*nblock + i)
                #print line.compute()
                idd = (delayed)(int)(line[0]) 
                x = delayed (float)(line[3])
                #print line[3].compute()
                y = (delayed)(float)(line[4])
                z = (delayed)(float)(line[5])
                if not delay:
                    idd = idd.compute()
                    x = x.compute()
                    y = y.compute()
                    z = z.compute()

                a = Atomc(idd,x,y,z)
                atoms.append(a)

            #create system
            sys = Systemc()
            sys.atoms = atoms
            sys.boxdims = boxdims
            #nsystems.append(sys)
            #pickle.dump(sys, fout)
            if compressed:
                np.savez(fout, [sys])
            else:
                np.save(fout, [sys])
        #fout.close()
    
        #now save pickled file
        #create a function
        return outfolder

def fetch_system(slice, outfolder=""):
    """
    Fetch a npy style system from outfolder.

    Parameters
    ----------
    outfolder : string
        the outfolder of the files
        if not specified outfolder is npydata folder
        in the work directory.
    slice : int
        number of slice

    Returns
    -------
    sys : Systemc
        The read system object
    """
    if outfolder == "":
        outfolder = os.path.join(os.getcwd(),"npydata")

    filefound = False
    if os.path.exists(outfolder):
        fout = ".".join(["snap",str(slice),"npy"])
        fout = os.path.join(outfolder, fout)
        if os.path.exists(fout):
            sys = np.load(fout)
            filefound = True
    
    if not filefound:
        sys = False

    return sys

def transfer_steinhardt(sys):
    """
    Transfer a steinhardt sys to Systemc sys
    """
    nsys = Systemc()
    satoms = sys.get_allatoms()
    nsys.natoms = len(atoms)

    nsys.atoms = []
    for atom in satoms:
        natom = Atomc(atom.id, atom.x, atom.y, atom.z)
        nsys.atoms.append(natom)
    nsys.boxdims = sys.get_box()

    return nsys

def untransfer_steinhardt(sys):
    """
    Transfer a Systemc sys to steinhardt sys
    """
    nsys = st.System()
    satoms = sys.atoms
    nsys.natoms = len(atoms)

    atoms = []
    for atom in satoms:
        natom = st.Atom()
        natom.id = atom.id 
        natom.x = atom.x 
        natom.y = atom.y 
        natom.z = atom.z

        atoms.append(natom)
    
    nsys.assign_particles(atoms, sys.boxdims)

    return nsys





def pickle_steinhardt(systems):
    """
    Pickle a steinhardt system to Systemc and save to disk
    """
    #move out from pickling and use npy arrays
    outfolder = os.path.join(os.getcwd(),"npydata")
    if os.path.exists(outfolder):
        shutil.rmtree(outfolder)
    os.mkdir(outfolder)

    for slice,sys in enumerate(systems):
        fout = ".".join(["snap",str(slice)])
        fout = os.path.join(outfolder, fout)


def unpickle_steinhardt():
    """
    Unpickle a steinhardt system. Converts a Systemc to steinhardt
    system and return. 
    """