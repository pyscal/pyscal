
from steinhardt import *
import numpy as np
#import sys
import mmap


def sum_nos(a,b):
    return a+b

def traj_to_systems3(filename,natoms):
    
    #filename = sys.argv[1]
    #natoms=int(sys.argv[2])




def traj_to_systems(filename,map=False):
    """
    Function to convert a trajfile to systems and return them.

    Parameters
    ----------
    filename : string
        the name of the input file
    map : bool
        use memory map if True, False otherwise

    Returns
    -------
    systems : array of System
        array consisting of a system for each time slice
    """
    atoms = []
    systems=[]

    with open(filename,"r+") as infile:
        if map:
            map = mmap.mmap(infile.fileno(),0)
            iterval = iter(map.readline,"")
        else:
            iterval = infile

        count = 1
        timestep = 0
        for line in iterval:
            #print count
            if count==4:
                natoms = int(line.strip())
            elif count==6:
                box = []
                raw = line.strip().split()
                dimx = float(raw[1])-float(raw[0])          
            
            elif count==7: 
                raw = line.strip().split()
                dimy = float(raw[1])-float(raw[0])
                
            elif count==8:
                raw = line.strip().split()
                dimz = float(raw[1])-float(raw[0])
                boxd = [dimx,dimy,dimz]
                
            elif count>9:
                line = line.strip().split()
                idd = int(line[0]) 
                x = float(line[3])
                y = float(line[4])
                z = float(line[5])
                xx = [x,y,z]
                a = Atom()
                a.sx(xx)
                a.sid(idd)  
                atoms.append(a)
                if count==natoms+9:
                    count=0
                    #do stuff here
                    sys = System()
                    sys.assign_particles(atoms,boxd)
                    systems.append(sys)
                    del atoms[:]
                    del boxd[:]
                    atoms = []
                    boxd = []
            count+=1
    return systems


#def systems_to_traj(systems,titles=None,values=None):
    """
    Function to write down a list of systems as traj,along with other per atom
    quantities.
    The other per quantities are provided by titles and values.
    As of now, the type of atom is irrelevant and hence is written.
    Also, velocities are ignored and set to zero.
    This can be easily modifies to read in extra quantities.
    For easiness, everything is handled as a string

    Parameters
    ----------
    systems : list of System
        list which contains all the system classes.
    titles: string
        this is the header values to be written
        for example, if you want to add three columns with titles ax, ay and az, the
        string to be passed would be-
        " ax ay az"
    values: 3d array of values
        This contains the per atom values to be written for each system
        the length along the first dimension must be equal to number of systems
        Each element must have 

    """

