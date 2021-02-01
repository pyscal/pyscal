import os
import numpy as np
from pyscal.formats.ase import convert_snap
import pyscal.core as pc
#import h5py
import warnings


#def hdf_to_dump(infile, outfile, keys=None):
"""
    A support function that can convert hdf formatted
    trajectory to dump format

    Parameters
    ----------
    infile : string
        name of the input hdf file

    outfile : string
        name of the output dump file

    keys : list, optional
        output keys to be written.
        default keys are box, [id, type, x, y, z]
    """
    """
    if keys is None:
        outkeys = ['x', 'y', 'z']
        mainkey = ['id', 'type', 'x', 'y', 'z']
    else:
        outkeys = np.concatenate((['x', 'y', 'z'], keys))
        mainkey = np.concatenate((['id', 'type', 'x', 'y', 'z'], keys))

    keyheader = " ".join(mainkey)
    keyheader = " ".join(["ITEM: ATOMS", keyheader, "\n"])

    with open(outfile, "w") as dump:
        with h5py.File(infile, "r") as hf:
            for key in hf.keys():
                
                natoms = len(np.array(hf[key]["atoms"]["x"]))
                box = np.array(hf[key]["box"])

                dump.write("ITEM: TIMESTEP\n")
                dump.write("%s\n" % key)
                dump.write("ITEM: NUMBER OF ATOMS\n")
                dump.write("%d\n" % natoms)
                dump.write("ITEM: BOX BOUNDS\n")
                dump.write("%f %f\n" % (box[0][0], box[0][1]))
                dump.write("%f %f\n" % (box[1][0], box[1][1]))
                dump.write("%f %f\n" % (box[2][0], box[2][1]))

                dump.write(keyheader)

                for i in range(natoms):
                    outval1 = ["%d %d"%(int(hf[key]["atoms"]["id"][i]), int(hf[key]["atoms"]["type"][i])) ]
                    outval2 = [str(hf[key]["atoms"][x][i]) for x in outkeys]
                    outvals = [*outval1, *outval2]
                    outvals.append("\n")
                    outline = " ".join(outvals)
                    dump.write(outline)
"""


class Timeslice:
    """
    Timeslice containing info about a single time slice
    Timeslices can also be added to each
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
                s = self.trajectories[count]._get_block_as_system(x, customkeys=customkeys)
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
                s = self.trajectories[count]._get_block_as_ase(x, species=species)
                sys.append(s)
        return sys

    def to_dict(self,):
        """
        Get the required block as data
        """
        data = []
        for count, traj in enumerate(self.trajectories):
            for x in self.blocklists[count]:
                self.trajectories[count].load(x)
                data.append(self.trajectories[count].data[x])
                self.trajectories[count].unload(x)
        return data

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
            self.trajectories[count]._get_blocks_to_file(fout, self.blocklists[count])
        fout.close()


    #def to_hdf(self, outfile, keys=None, mode='w', compression="gzip"):
    """
        Get the block as hdf file

        Parameters
        ----------
        outfile : string
            name of the output file

        headers : list, optional
            The keys to be stored in hdf format.
            Default values stored are [id, type, x, y, z]
        
        mode : string, optional
            h5 write mode.
            see here - https://docs.h5py.org/en/stable/high/file.html
        
        compression : string, optional
            the compression algorithm to choose. Default gzip

        Returns
        -------
        None
        """
        """
        if keys is None:
            outkeys = ['id', 'type', 'x', 'y', 'z']
        else:
            outkeys = np.concatenate((['id', 'type', 'x', 'y', 'z'], keys))

        c = 0
        with h5py.File(outfile, 'w') as hf:
            for count, traj in enumerate(self.trajectories):
                for x in self.blocklists[count]:
                    self.trajectories[count].load(x)
                    data = self.trajectories[count].data[x]
                    self.trajectories[count].unload(x)
                    tk = str(c)
                    #warnings.warn(tk)
                    hf.create_group(tk)
                    hf[tk].create_dataset('box', data=data['box'], compression=compression)
                    hf[tk].create_group("atoms")
                    for key in outkeys:
                       hf[tk]["atoms"].create_dataset(key, data=data['atoms'][key], compression=compression)
                    c += 1
    """

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

        self._get_natoms()
        self._get_nblocks()
        
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

    def _get_natoms(self):
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

    def _get_nlines(self):
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
    
    def _get_nblocks(self):
        """
        Get number of blocks in the trajectory file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self._get_natoms()
        nlines = self._get_nlines()
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

    def _convert_data_to_lines(self, blockno):
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

    def _get_block_as_system(self, blockno, customkeys=None):
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

    def _get_block_as_ase(self, blockno, species=None):
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

    def _get_blocks_to_file(self, fout, blocklist):
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
                data = self._convert_data_to_lines(x)
            else:
                data = self.get_block(x)
            for line in data:
                fout.write(line)
