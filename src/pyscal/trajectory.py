import os
import pyscal.core as pc

class Timeslice:
    """
    Timeslice containing info about a single time slice
    """
    def __init__(self, data):
        """
        Initialize instance with data
        """
        pass

class Trajectory:
    """
    A Trajectory class for LAMMPS
    """
    def __init__(self, filename, customkeys=None):
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
        self.customkeys = customkeys
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
            blocklist = blockno.indices(self.nblocks)
            sys = [self.get_block_as_system(x) for x in blocklist]
            return sys
        else:
            sys = self.get_block_as_system(blockno)
            return sys

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
        nlines = sum(1 for i in open(self.filename, 'rb'))
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
        #now we can get the lines quickly
        with open(self.filename, "rb") as fout:
            for i in range(start):
                _ = next(fout)
            data = [next(fout).decode("utf-8") for x in range(start, stop)]
        return data

    def get_block_as_system(self, blockno):
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
        sys.read_inputfile(data, customkeys=self.customkeys)
        return sys

    def get_block_as_file(self, blockno, outfile, format="lammps-dump"):
        """
        Get block as outputfile

        Parameters
        ----------
        blockno : int
            number of the block to be read, starts from 0

        Returns
        -------
        None

        """
        #convert to system and return
        sys = self.get_block_as_system(blockno)
        sys.to_file(outfile, customkeys=self.customkeys, format=format)

