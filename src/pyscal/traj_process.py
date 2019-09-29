"""
pyscal module containing methods for processing of a trajectory. Methods for
reading of input files formats, writing of output files etc are provided in
this module.

"""
import numpy as np
import gzip
import pyscal.catom as pca


#functions that are not wrapped from C++
def read_lammps_dump(infile, compressed = False, check_triclinic=False, box_vectors=False, customkeys=[]):
    """
    Function to read a lammps dump file format - single time slice.

    Parameters
    ----------
    infile : string
        name of the input file

    compressed : bool, optional
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary. Default True.

    check_triclinic : bool, optional
        If true check if the sim box is triclinic. Default False.

    box_vectors : bool, optional
        If true, return the full box vectors along with `boxdims` which gives upper and lower bounds.
        default False.

    customkeys : list of strings, optional
        A list of extra keywords to read from trajectory file.

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input

    boxdims : list of list of floats
        The dimensions of the box. This is of the form `[[xlo, xhi],[ylo, yhi],[zlo, zhi]]` where `lo` and `hi` are
        the upper and lower bounds of the simulation box along each axes. For triclinic boxes, this is scaled to
        `[0, scalar length of the vector]`.

    box : list of list of floats
        list of the type `[[x1, x2, x3], [y1, y2, y3], [zz1, z2, z3]]` which are the box vectors. Only returned if
        `box_vectors` is set to True.

    triclinic : bool
        True if the box is triclinic. Only returned if `check_triclinic` is set to True

    .. note::

        Values are always returned in the order `atoms, boxdims, box, triclinic` if all
        return keywords are selected. For example, ff `check_triclinic` is not selected, the return
        values would still preserve the order and fall back to  `atoms, boxdims, box`.

    Notes
    -----
    Read a lammps-dump style snapshot that can have variable headers, reads in type and so on.
    Zipped files which end with a `.gz` can also be read automatically. However, if the file does not
    end with a `.gz` extension, keyword `compressed = True` can also be used.

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
    triclinic = False
    volume_fraction = 1.00

    #now if custokeys are provided - read those in too
    customread = False
    customlength = len(customkeys)
    if customlength > 0:
        customread = True


    for count, line in enumerate(f):
        if not paramsread:
            #atom numer is at line 3
            if count == 3:
                natoms = int(line.strip())
            #box dims in lines 5,6,7
            elif count == 5:
                raw = line.strip().split()
                boxx = [float(raw[0]), float(raw[1])]
                if len(raw) == 3:
                    xy = float(raw[2])
            elif count == 6:
                raw = line.strip().split()
                boxy = [float(raw[0]), float(raw[1])]
                if len(raw) == 3:
                    xz = float(raw[2])
            elif count == 7:
                raw = line.strip().split()
                boxz = [float(raw[0]), float(raw[1])]
                if len(raw) == 3:
                    yz = float(raw[2])
                    triclinic = True
                    tilts = [xy, xz, yz]
                #boxdims = [boxx, boxy, boxz]

            #header is here
            elif count == 8:
                raw = line.strip().split()
                headerdict = { raw[x]:x-2 for x in range(0, len(raw)) }
                paramsread = True
                if customread:
                    if not np.prod([(x in headerdict) for x in customkeys]):
                        raise KeyError("some values in custokeys was not found in the file")


        elif count > 8:
            raw = line.strip().split()
            idd = int(raw[headerdict["id"]])
            typ = int(raw[headerdict["type"]])
            x = float(raw[headerdict["x"]])
            y = float(raw[headerdict["y"]])
            z = float(raw[headerdict["z"]])

            atom = pca.Atom()
            atom.pos = [x, y, z]
            atom.id = idd
            atom.type = typ
            atom.loc = count-8

            customdict = {}
            #if customkeys need to be read, do it
            if customread:
                for cc, kk in enumerate(customkeys):
                    customdict[kk] = raw[headerdict[kk]]

            atom.custom = customdict
            atoms.append(atom)

    #close files
    f.close()

    if triclinic:
        #process triclinic box
        amin = min([0.0, tilts[0], tilts[1] ,tilts[0]+tilts[1]])
        amax = max([0.0, tilts[0], tilts[1] ,tilts[0]+tilts[1]])
        bmin = min([0.0, tilts[2]])
        bmax = max([0.0, tilts[2]])
        xlo = boxx[0] - amin
        xhi = boxx[1] - amax
        ylo = boxy[0] - bmin
        yhi = boxy[1] - bmax
        zlo = boxz[0]
        zhi = boxz[1]

        #triclinic cell
        a = np.array([xhi-xlo, 0, 0])
        b = np.array([tilts[0], yhi-ylo, 0])
        c = np.array([tilts[1], tilts[2], zhi-zlo])

        rot = np.array([a, b, c]).T
        rotinv = np.linalg.inv(rot)
        ortho_origin = np.array([boxx[0], boxy[0], boxz[0]])

        for atom in atoms:
            #correct zero of the atomic positions (shift box to origin)
            dist = np.array(atom.pos) - ortho_origin
            atom.pos = dist

        #finally change boxdims - to triclinic box size
        box = np.array([a, b, c])
        boxdims = np.array([[0, np.sqrt(np.sum(a**2))],[0, np.sqrt(np.sum(b**2))],[0, np.sqrt(np.sum(c**2))]])
    else:
        box = np.array([[boxx[1]-boxx[0], 0, 0],[0, boxy[1]-boxy[0], 0],[0, 0, boxz[1]-boxz[0]]])
        boxdims = np.array([[boxx[0], boxx[1]],[boxy[0], boxy[1]],[boxz[0], boxz[1]]])

    if box_vectors and check_triclinic:
        return atoms, boxdims, box, triclinic

    elif box_vectors:
        return atoms, boxdims, box

    elif check_triclinic:
        return atoms, boxdims, triclinic
    else:
        return atoms, boxdims

def read_poscar(infile, compressed = False):
    """
    Function to read a POSCAR format.

    Parameters
    ----------
    infile : string
        name of the input file

    compressed : bool, optional
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary, Default False

    Returns
    -------
    atoms : list of `Atom` objects
        list of all atoms as created by user input

    box : list of list of floats
        list of the type `[[xlow, xhigh], [ylow, yhigh], [zlow, zhigh]]` where each of them are the lower
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
            atom = pca.Atom()
            atom.pos = [x, y, z]
            atom.id = idd
            atom.type = typ
            atom.loc = i-atom_start
            #atom = pc.Atom(pos=, id=idd, type=typ)
            atoms.append(atom)
            i+=1
        else:
            species+=1

    return atoms, boxdims

def write_structure(sys, outfile, format = 'lammps-dump', compressed = False, customkey=None, customvals=None):
    """
    Write the state of the system to a trajectory file.

    Parameters
    ----------
    sys : `System` object
        the system object to be written out

    outfile : string
        name of the output file

    format : string, optional
        the format of the output file, as of now only `lammps-dump` format
        is supported.

    compressed : bool, default false
        write a `.gz` format

    customkey : string or list of strings, optional
        If specified, it adds this custom column to the dump file. Default None.

    customvals : list or list of lists, optional
        If `customkey` is specified, `customvals` take an array of the same length
        as number of atoms, which contains the values to be written out.

    Returns
    -------
    None

    """
    boxdims = sys.box
    atoms = sys.atoms

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


    #check customkey and its values
    #if single value, make it into array
    writecustom = False

    if customkey != None:
        if not isinstance(customkey, list) or  isinstance(customkey, np.ndarray):
            customkey = [str(customkey)]
            customvals = np.array([customvals])
        else:
            ccdummy = [str(x) for x in customkey]
            customkey = ccdummy

        #verify lengths
        if not len(customkey) == len(customvals):
            #print(cust)
            #raise TypeError("length of customkey and customvals should be same. lengths %d and %d not compatible"%(len(customkey), len(customvals)))
            raise TypeError(customkey, customvals)

        #now verify the lengths of customkey
        for count, c in enumerate(customkey):
            if not len(customvals[count]) == len(atoms):
                raise TypeError("Length of customvals should be equal to number of atoms. lengths %d and %d not compatible"%(len(customvals[count]), len(atoms)))

        #if everything works - change writecustom
        writecustom = True

    #now write header
    if writecustom:
        ckey = " ".join(customkey)
        title_str = "ITEM: ATOMS id type x y z %s\n"% ckey
    else:
        title_str = "ITEM: ATOMS id type x y z\n"

    #write it out
    dump.write(title_str)

    for cc, atom in enumerate(atoms):
        pos = atom.pos

        if writecustom:
            cvals = " ".join(np.array(customvals)[:, cc].astype(str))
            atomline = ("%d %d %f %f %f %s\n")%(atom.id, atom.type, pos[0], pos[1], pos[2], cvals)
        else:
            atomline = ("%d %d %f %f %f\n")%(atom.id, atom.type, pos[0], pos[1], pos[2])

        dump.write(atomline)

    dump.close()


def split_trajectory(infile, format='lammps-dump', compressed=False):
    """
    Read in a trajectory file and convert it to individual time slices.

    Parameters
    ----------

    filename : string
        name of input file

    format : format of the input file
        only `lammps-dump` is supported now.

    compressed : bool, optional
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary.

    Returns
    -------
    snaps : list of strings
        a list of filenames which contain individual frames from the main trajectory.

    Notes
    -----
    This is a wrapper function around `split_traj_lammps_dump` function.

    """

    snaps = []

    if format=='lammps-dump':
        snaps = split_traj_lammps_dump(infile, compressed = compressed)

    return snaps


def split_traj_lammps_dump(infile, compressed = False):
    """
    Read in a LAMMPS dump trajectory file and convert it to individual time slices.

    Parameters
    ----------

    filename : string
        name of input file

    compressed : bool, optional
        force to read a `gz` zipped file. If the filename ends with `.gz`, use of this keyword is not
        necessary, Default False.

    Returns
    -------
    snaps : list of strings
        a list of filenames which contain individual frames from the main trajectory.

    """
    raw = infile.split('.')
    if raw[-1] == 'gz' or  compressed:
        f = gzip.open(infile,'rt')
    else:
        f = open(infile,'r')


    #pre-read to find the number of atoms
    for count, line in enumerate(f):
            if count == 3:
                natoms = int(line.strip())
                break
    f.close()

    #now restart f()
    if raw[-1] == 'gz' or  compressed:
        f = gzip.open(infile,'rt')
    else:
        f = open(infile,'r')


    nblock = natoms+9
    startblock = 0
    count=1
    snaps = []



    for line in f:
        if(count==1):
            ff = ".".join([infile, 'snap', str(startblock), 'dat'])
            snaps.append(ff)
            fout = open(ff,'w')
            fout.write(line)

        elif(count<nblock):
            fout.write(line)

        else:
            fout.write(line)
            fout.close()
            count=0
            startblock+=1
        count+=1

    f.close()

    return snaps
