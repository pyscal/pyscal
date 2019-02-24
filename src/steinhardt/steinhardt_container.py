import numpy as np
import dask.bag as db
from dask import delayed
import os

class Atomc(object):
    def __init__(self,idd,x,y,z):
        self.id = idd
        self.x = x
        self.y = y
        self.z = z

class Boxc(object):
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

def create_delayed_system(infile, natoms, **kwargs):
    """
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
        save_file (bool) : True, True if an outfile is to be saved, False otherwise
        return_array (bool) : False, True if array is to be returned, False otherwise
        nslices (int) : 1, number of time slices in the trajectory
        format (string) : lammps, the format of the trajectory file
        outfile (string) : name of the output file

    Returns
    -------
    nsystems : array like
        An array of box and atom information for each time slices. Only returned if
        return_array is True.

    """
    #read in the dask bag and convert to delayed object
    b = db.read_text(infile)
    c = b.to_delayed()

    #initialise systems
    nsystems = []

    #process kwargs
    save_file = kwargs.get('save_file', True)
    return_array = kwargs.get('return_array', False)
    nslices = kwargs.get('nslices', 1)
    format = kwargs.get('format', "lammps")
    outfile = kwargs.get('outfile', os.path.join(os.getcwd(),".".join([infile,"npy"])))

    #if natoms is a single value - make it into an array. This allows for providing
    #variable atom numbers in each time slice
    if not isinstance(natoms, list):
        natoms = np.ones(nslices)*natoms

    #this part would be  md code specific
    if format ==  "lammps":
        for slice in range(nslices):
            nblock = natoms[slice] + 9
            raw = c[0][nslices*nblock + 5].strip().split()
            dimx = (delayed)(float)(raw[1])-(delayed)(float)(raw[0])          
            raw = c[0][nslices*nblock + 6].strip().split()
            dimy = (delayed)(float)(raw[1])-(delayed)(float)(raw[0])    
            raw = c[0][nslices*nblock + 7].strip().split()
            dimz = (delayed)(float)(raw[1])-(delayed)(float)(raw[0])
            boxd = [dimx,dimy,dimz]
            box = Boxc(dimx, dimy, dimz)
            atoms = []

            for i in range(9,natoms+9):
                line = c[0][nslices*nblock + i].strip().split()
                idd = (delayed)(int)(line[0]) 
                x = (delayed)(float)(line[3])
                y = (delayed)(float)(line[4])
                z = (delayed)(float)(line[5])
                a = Atomc(idd,x,y,z)
                atoms.append(a)
                nsystems.append([box,atoms])
    
    #additional codes can be added at this point
    
    #process the output
    if save_file:
        np.save(outfile, nsystems)

    if return_array:
        return nsystems