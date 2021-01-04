import os
import numpy as np
from pyscal.formats.ase import convert_snap

class Timeslice:
    """
    Timeslice containing info about a single time slice
    """
    def __init__(self, trajectory, blocklist):
        """
        Initialize instance with data
        """
        self.trajectory = trajectory
        self.blocklist = blocklist
        self.trajectories = [trajectory]
        self.blocklists = [blocklist]

    def __repr__(self):
        """
        String of the class
        """
        blockstring = ["%d-%d"%(x[0], x[-1]) for x in self.blocklists]
        blockstring = "/".join(blockstring)

        data = "Trajectory slice\n %s\n natoms=%d\n"%(blockstring, self.trajectory.natoms)
        return data

    def __add__(self, ntraj):
        """
        Add a method for summing
        """
        for traj in ntraj.trajectories:
            self.trajectories.append(traj)

        for block in ntraj.blocklists:
            self.blocklists.append(block)

        return self


    def __radd__(self, ntraj):
        """
        Reverse add method
        """
        if ntraj == 0:
            return self
        else:
            return self.__add__(ntraj)


    def to_system(self, customkeys=None):
        """
        Get block as pyscal system

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        sys : System
            pyscal System
        """
        sys = []
        for count, traj in enumerate(self.trajectories):
            for x in self.blocklists[count]:
                s = self.trajectories[count].get_block_as_system(x, customkeys=customkeys)
                sys.append(s)
        return sys

    def to_ase(self, species=None):
        """
        Get block as Ase objects

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        sys : ASE object
            
        """
        sys = []
        for count, traj in enumerate(self.trajectories):
            for x in self.blocklists[count]:
                s = self.trajectories[count].get_block_as_ase(x, species=species)
                sys.append(s)
        return sys

    def to_file(self, outfile, mode="w"):
        """
        Get block as outputfile

        Parameters
        ----------
        outfile : string
            name of output file

        mode : string
            write mode to be used, optional
            default "w" write
            also can be "a" to append.
        
        Returns
        -------
        None

        """
        fout = open(outfile, mode)
        for count, traj in enumerate(self.trajectories):
            self.trajectories[count].get_blocks_to_file(fout, self.blocklists[count])
        fout.close()


class Trajectory:
    """
    A Trajectory class for LAMMPS
    """
    def __init__(self, filename):
        """
        Initiaze the class

        Parameters
        ----------
        filename : string
            name of the inputfile

        customkeys : list of string
            keys other than position, id that needs to be read
            in from the input file
        """
        if os.path.exists(filename):
            self.filename = filename
        else:
            raise FileNotFoundError("%s file not found"%filename)
        self.natoms = 0
        self.blocksize = 0
        self.nblocks = 0
        self.loadlist = None
        self.data = None

        self.get_natoms()
        self.get_nblocks()
        
    def __repr__(self):
        """
        String of the class
        """
        return "Trajectory of %d slices with %d atoms"%(self.nblocks, self.natoms)

    def __getitem__(self, blockno):
        """
        Allow for slice indexing
        """
        if isinstance(blockno, slice):
            blocklist = range(*blockno.indices(self.nblocks))
            timeslice = Timeslice(self, blocklist)
            return timeslice
        else:
            blocklist = [blockno]
            timeslice = Timeslice(self, blocklist)
            return timeslice

    def get_natoms(self):
        """
        Get number of atoms in the system

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        with open(self.filename, "rb") as fout:
            data = [next(fout) for x in range(0, 4)]
        self.natoms = (int(data[-1]))

    def get_nlines(self):
        """
        Get total number of lines in the file

        Parameters
        ----------
        None

        Returns
        -------
        nlines : int
            number of lines
        """
        line_offset = []
        offset = 0
        nlines = 0
        for line in open(self.filename, 'rb'):
            line_offset.append(offset)
            offset += len(line)
            nlines += 1
        
        self.nlines = nlines
        self.line_offset = line_offset
        return nlines
    
    def get_nblocks(self):
        """
        Get number of blocks in the trajectory file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.get_natoms()
        nlines = self.get_nlines()
        self.blocksize = self.natoms+9
        self.nblocks = nlines//self.blocksize
        self.straylines = nlines - self.nblocks*self.blocksize
        #set load list to False
        self.loadlist = [False for x in range(self.nblocks)]
        self.data = [None for x in range(self.nblocks)]

    def get_block(self, blockno):
        """
        Get a block from the file as raw data

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        data : list
            list of strings containing data
        """
        start = blockno*self.blocksize
        stop = (blockno+1)*self.blocksize

        fin = open(self.filename, "rb")
        fin.seek(0)
        fin.seek(self.line_offset[start])

        data = []
        for i in range(self.blocksize):
            line = fin.readline().decode("utf-8")
            data.append(line)
        return data

    def load(self, blockno):
        """
        Load the data of a block into memory as a dictionary
        of numpy arrays

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        None

        Notes
        -----
        When the data of a block is loaded, it is accessible
        through `Trajectory.data[x]`. This data can then be
        modified. When the block is written out, the modified
        data is written instead of existing one. But, loaded
        data is kept in memory until unloaded using `unload`
        method.
        """
        data = self.get_block(blockno)
        box =  np.loadtxt(data[5:8])
        columns = np.loadtxt(data[9:])
        header = np.loadtxt(data[8:9], dtype=str)[2:]
        outdict = {}
        outdict["box"] = box
        outdict["atoms"] = {}
        for count, h in enumerate(header):
            outdict["atoms"][h] = columns[:,count]        

        self.data[blockno] = outdict
        self.loadlist[blockno] = True

    def unload(self, blockno):
        """
        Unload the data that is loaded to memory using
        `load` method

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        None        
        """
        self.data[blockno] = None
        self.loadlist[blockno] = False        

    def convert_data_to_lines(self, blockno):
        """
        Create lines from loaded data
        
        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0
        
        Returns
        -------
        data : list of strs
            list of lines
        """
        dd = self.data[blockno]

        data = []
        data.append("ITEM: TIMESTEP\n")
        data.append("".join([str(0), os.linesep]))
        data.append("ITEM: NUMBER OF ATOMS\n")
        data.append("".join([str(self.natoms), os.linesep]))
        data.append("ITEM: BOX BOUNDS pp pp pp\n")
        for b in dd["box"]:
            dstr = " ".join(b.astype(str))
            data.append("".join([dstr, os.linesep]))

        xf = []
        xd = []
        xfkeys = []
        xdkeys = []
        for key, val in dd["atoms"].items():
            if key in ["id", "type"]:
                val.astype(int)
                xd.append(val)
                xdkeys.append(key)
            else:
                xf.append(val)
                xfkeys.append(key)

        xdstrs = []
        if len(xd)>0:
            for i in range(len(xd[0])):
                substr = []
                for j in range(len(xdkeys)):
                    substr.append("%d"%xd[j][i])
                xdstrs.append(" ".join(substr))

        xdheader = " ".join(xdkeys)
        xdheader = " ".join(["ITEM: ATOMS", xdheader])

        xfstrs = []
        xf = np.array(xf)
        if len(xf)>0:
            for i in range(len(xf[0])):
                dstr = " ".join((xf[:,i]).astype(str))
                xfstrs.append("".join([dstr, os.linesep]))

        xfheader = " ".join(xfkeys)
        mainheader = " ".join([xdheader, xfheader])
        mainheader = "".join([mainheader, os.linesep])

        data.append(mainheader)

        for i in range(len(xfstrs)):
            valstr = " ".join([xdstrs[i], xfstrs[i]])
            #valstr = "".join([valstr, os.linesep])
            data.append(valstr)

        return data        

    def get_block_as_system(self, blockno, customkeys=None):
        """
        Get block as pyscal system
        
        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0
        
        Returns
        -------
        sys : System
            pyscal System
        """
        #convert to system and return
        data = self.get_block(blockno)

        sys = pc.System()
        sys.read_inputfile(data, customkeys=customkeys)
        return sys

    def get_block_as_ase(self, blockno, species=None):
        """
        Get block as pyscal system
        
        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0
        
        Returns
        -------
        sys : System
            pyscal System
        """
        #convert to system and return
        data = self.get_block(blockno)

        sys = pc.System()
        sys.read_inputfile(data, customkeys=None)
        asesys = convert_snap(sys, species=species)
        return asesys

    def get_blocks_to_file(self, fout, blocklist):
        """
        Get a series of blocks from the file as raw data

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        data : list
            list of strings containing data
        """
        xl = [x for x in blocklist]
        xl = np.array(xl)
                
        #open file
        #convert lines to start from end
        for x in xl:
            if self.loadlist[x]:
                data = self.convert_data_to_lines(x)
            else:
                data = self.get_block(x)
            for line in data:
                fout.write(line)
