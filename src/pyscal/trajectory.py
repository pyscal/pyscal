import os
import pyscal.core as pc
import numpy as np

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

    def to_file(self, outfile):
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
        if os.path.exists(outfile):
            os.remove(outfile)

        fout = open(outfile, "a")
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
        self.straylines = lines - self.nblocks*self.blocksize

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
        
        #get start and stop lines
        start = xl*self.blocksize
        stop = (xl+1)*self.blocksize

        args = np.argsort(start)
        sortedstart = start[args]
        sortedstop = stop[args]
        
        #check if things need to be reversed
        reverse = False
        if (args[0]!=0):
            reverse = True

        if reverse:
            #open file
            firstline = 0
            with open(self.filename, "rb") as fin:
                for count, i in enumerate(start):
                    secondline = sortedstart[count]
                    thirdline = sortedstop[count]
                    #discard lines
                    for j in range(firstline, secondline):
                        _ = next(fin)
                    #now write the rest
                    for j in range(secondline, thirdline):
                        line = next(fin).decode("utf-8")
                        fout.write(line)
                    #we have to seek back to beginning of the file
                    fout.seek(0)
                    #now we have to reset the variables
                    firstline = sortedstop[count]
        else:
            firstline = 0
            with open(self.filename, "rb") as fin:
                for count, i in enumerate(start):
                    secondline = start[count]
                    thirdline = stop[count]
                    #discard lines
                    for j in range(firstline, secondline):
                        _ = next(fin)
                    #now write the rest
                    for j in range(secondline, thirdline):
                        line = next(fin).decode("utf-8")
                        fout.write(line)
                    #now we have to reset the variables
                    firstline = stop[count]