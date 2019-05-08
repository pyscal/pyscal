"""
This file contains the methods for processing of a trajectory

"""
import pybop.ccore as pc
import numpy as np
import gzip

#functions that are not wrapped from C++
def read_lammps_dump(infile, compressed = False):
    """
    Function to read a lammps dump file format - single time slice. This is essentially the same as 
    C++ read, but can have variable headers, reads in type and so on. So if verstility is required,
    this function is better. C++ version might be faster.
    Zipped files which end with a `.gz` can also be read automatically. However, if the file does not
    end with a `.gz` extension, keyword `compressed = True` can also be used.

    Parameters
    ----------
    infile : string
        name of the input file
    compressed : bool, Default False
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary.

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input
    box : list of list of floats
        list of the type [[xlow, xhigh], [ylow, yhigh], [zlow, zhigh]] where each of them are the lower
        and upper limits of the simulation box in x, y and z directions respectively.

    Examples
    --------
    >>> atoms, box = read_lammps_dump('conf.dump')
    >>> atoms, box = read_lammps_dump('conf.dump.gz')
    >>> atoms, box = read_lammps_dump('conf.d', compressed=True)    

    """
    #first depending on if the extension is .gz - use zipped read
    raw = infile.split('.')
    if raw[-1] == 'gz' or  compressed:
        f = gzip.open(infile,'rt')          
    else:
        f = open(infile,'r')

    #now go through the file line by line            
    paramsread = False
    atoms = []

    for count, line in enumerate(f):
        if not paramsread:
            #atom numer is at line 3
            if count == 3:
                natoms = int(line.strip())
            #box dims in lines 5,6,7
            elif count == 5:
                raw = line.strip().split()
                boxx = [float(raw[0]), float(raw[1])]
            elif count == 6:
                raw = line.strip().split()
                boxy = [float(raw[0]), float(raw[1])]
            elif count == 7:
                raw = line.strip().split()
                boxz = [float(raw[0]), float(raw[1])]
                boxdims = [boxx, boxy, boxz]
            #header is here
            elif count == 8:
                raw = line.strip().split()
                headerdict = { raw[x]:x-2 for x in range(0, len(raw)) }
                paramsread = True
        elif count > 8:
            raw = line.strip().split()
            idd = int(raw[headerdict["id"]])
            typ = int(raw[headerdict["type"]])
            x = float(raw[headerdict["x"]])
            y = float(raw[headerdict["y"]])
            z = float(raw[headerdict["z"]])
            atom = pc.Atom()
            atom.set_x([x, y, z])
            atom.set_id(idd)
            atom.set_type(typ)
            atoms.append(atom)

    #close files
    f.close()

    return atoms, boxdims


def read_poscar(infile, compressed = False):
    """
    Function to read a POSCAR format. Zipped files which end with a `.gz` can also be 
    read automatically. However, if the file does not end with a `.gz` extension, 
    keyword `compressed = True` can also be used.

    Parameters
    ----------
    infile : string
        name of the input file
    compressed : bool, Default False
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary.

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input
    box : list of list of floats
        list of the type [[xlow, xhigh], [ylow, yhigh], [zlow, zhigh]] where each of them are the lower
        and upper limits of the simulation box in x, y and z directions respectively.

    Examples
    --------
    >>> atoms, box = read_poscar('POSCAR')
    >>> atoms, box = read_poscar('POSCAR.gz')
    >>> atoms, box = read_poscar('POSCAR.dat', compressed=True)    

    """
    raw = infile.split('.')
    if raw[-1] == 'gz' or  compressed:
        f = gzip.open(infile,'rt')            
    else:
        f = open(infile,'r')
    
    data = []
    for line in f:
        data.append(line)

    no_atoms = data[5].split()
    no_atoms = np.array(no_atoms)
    no_atoms = no_atoms.astype(int)

    natoms = np.sum(no_atoms)
    atom_list = no_atoms        

    scaling_factor = np.array(data[1].strip()).astype(float)
    xvector = np.array(data[2].strip().split()).astype(float)
    yvector = np.array(data[3].strip().split()).astype(float)
    zvector = np.array(data[4].strip().split()).astype(float)

    xlow = 0
    xhigh = scaling_factor*xvector[0]    
    ylow = 0
    yhigh = scaling_factor*yvector[1]    
    zlow = 0
    zhigh = scaling_factor*zvector[2]
    boxdims = [[xlow, xhigh],[ylow, yhigh],[zlow, zhigh]]    

    if (data[6].strip().split()[0]=='s' or data[6].strip().split()[0]=='S'):
        selective_dynamics=True
        cord_system=data[7].strip()
        atom_start = 8
    else:
        cord_system=data[6].strip()
        atom_start = 7

    species = 1
    count = 0

    cum_list = np.cumsum(atom_list)
    i = atom_start
    atoms = []

    while i in range(atom_start,atom_start+natoms):
        if (count<cum_list[species-1]):
            raw = np.array(data[i].strip().split()).astype(float)
            typ = species
            x = float(raw[0])*xhigh
            y = float(raw[1])*yhigh
            z = float(raw[2])*zhigh
            #if x,y,z are out of the box, they need to be put in
            if (x < xlow):
                x = x + (xhigh - xlow)
            elif (x > xhigh):
                x = x - (xhigh - xlow)
            if (y < ylow):
                y = y + (yhigh - ylow)
            elif (y > yhigh):
                y = y - (yhigh - ylow)
            if (z < zlow):
                z = z + (zhigh - zlow)
            elif (z > zhigh):
                z = z - (zhigh - zlow)

            count+=1
            idd = count
            atom = pc.Atom()
            atom.set_x([x, y, z])
            atom.set_id(idd)
            atom.set_type(typ)
            #atom = pc.Atom(pos=, id=idd, type=typ)
            atoms.append(atom)
            i+=1
        else:
            species+=1

    return atoms, boxdims

def write_structure(sys, outfile, format = 'lammps-dump', compressed = False):
    """
    Write the state of the system to a trajectory file. 

    Parameters
    ----------
    sys : `System` object
        the system object to be written out
    outfile : string
        name of the output file
    format : string, default `lammps-dump`
        the format of the output file, as of now only `lammps-dump` format
        is supported.
    compressed : bool, default false
        write a `.gz` format

    Returns
    -------
    None

    """
    boxdims = sys.get_box()
    atoms = sys.get_allatoms()

    #open files for writing
    if compressed:
        gz = gzip.open(outfile,'w')
        dump = io.BufferedReader(gz)            
    else:
        gz = open(outfile,'w')
        dump = gz

    #now write
    dump.write("ITEM: TIMESTEP\n")
    dump.write("0\n")
    dump.write("ITEM: NUMBER OF ATOMS\n")
    dump.write("%d\n" % len(atoms))
    dump.write("ITEM: BOX BOUNDS\n")
    dump.write("%f %f\n" % (boxdims[0][0], boxdims[0][1]))
    dump.write("%f %f\n" % (boxdims[1][0], boxdims[1][1]))
    dump.write("%f %f\n" % (boxdims[2][0], boxdims[2][1]))
    dump.write("ITEM: ATOMS id type x y z\n")
    for atom in atoms:
        pos = atom.get_x()
        dump.write(("%d %d %f %f %f\n") %
                   (atom.get_id(), atom.get_type(), pos[0], pos[1], pos[2]))

    dump.close()


"""
def traj_to_systems(filename,natoms,snaps=False,npy=False,zipped=False,usecols=[]):
"""
"""
    DEPRECATED
    Reads in LAMMPS trajectory file and converts it to systems.

    Parameters
    ----------
    
    filename : string
        name of input file
    natoms : int
        number of input atoms
    snaps : bool
        splits the input file into smaller files of one time slice each and
        return a list of filenames if True.
        If False, returns an array of systems with the atom and box coordinates
        read in.
    npy : bool
        If true saves the systems as a npy file.
    usecols : array like
        If custom column numbers are to be used, they can be specified here.
        for example if usecols = [3,4,5,8,9] - the first three columns are
        used to read in x,y, and z values are the rest are read in as custom
        properties. They can be accessed by atom.gcustom(). The minimum length
        of usecols has to be three(for positions). The column number starts 
        from 0. Atom ID always has to be column 0.
        Restrictions - * ONLY used if snaps is set to False *

    Returns
    -------
    
    files : array like
        list of files names of one time slice each. Only returned if
        snaps=True

"""
"""
    if zipped:
        gz = gzip.open(filename,'rb')
        infile = io.BufferedReader(gz)
    else:
        gz = open(filename,'r')
        infile = gz        

    if snaps:
        nblock = natoms+9
        startblock = 0
        count=1
        files = []


        for line in infile:
            if(count==1):
                ff = ".".join([filename, 'snap', str(startblock), 'dat'])
                #ff = 'snap.'+str(startblock)+'.dat'
                files.append(ff)
                fout = open(ff,'w')
                fout.write(line)
                #count+=1
            elif(count<nblock):
                fout.write(line)
                #count+=1

            else:
                fout.write(line)
                fout.close()
                count=0
                startblock+=1
            count+=1
        infile.close()

        #qtraj = getQTrajectory(files)
        #systems=[]
        #for file in files:
            #sys = System()
            #sys.inputfile = file
            #systems.append(sys)
        return files

    else:
        atoms = []
        systems=[]

        custom = False
        if len(usecols) > 0:
            if len(usecols)>=3:
                xid = usecols[0]
                yid = usecols[1]
                zid = usecols[2]
                customcols = usecols[3:]
                custom = True
        else:
            xid = 3
            yid = 4
            zid = 5

        iterval = infile
        count = 1
        timestep = 0
        for line in iterval:
            if count==6:
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
                x = float(line[xid])
                y = float(line[yid])
                z = float(line[zid])
                cavals = []
                if custom:
                    for cindex in customcols:
                        cavals.append(float(line[cindex]))

               
                xx = [x,y,z]
                a = Atom()
                a.x = x
                a.y = y
                a.z = z
                a.id = idd
                a.scustom(cavals)
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
        infile.close()
        
        if npy:
            np.save(".".join([filename,"npy"]),systems)
        return systems        
"""
