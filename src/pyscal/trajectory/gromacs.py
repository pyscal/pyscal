# -*- coding: utf-8 -*-
# Copyright (c) 2017, Anders Lervik.
# Used with permission from Anders Lervik.

import numpy as np
import gzip
import pyscal.catom as pca
from ase import Atom, Atoms
import gzip
import io
import struct


GROMACS_MAGIC = 1993
DIM = 3
TRR_VERSION = 'GMX_trn_file'
TRR_VERSION_B = b'GMX_trn_file'
SIZE_FLOAT = struct.calcsize('f')
SIZE_DOUBLE = struct.calcsize('d')
HEAD_FMT = '{}13i'
HEAD_ITEMS = ('ir_size', 'e_size', 'box_size', 'vir_size', 'pres_size',
              'top_size', 'sym_size', 'x_size', 'v_size', 'f_size',
              'natoms', 'step', 'nre', 'time', 'lambda')
DATA_ITEMS = ('box_size', 'vir_size', 'pres_size',
              'x_size', 'v_size', 'f_size')

def read_snap(filename):
	"""
	Read a trr file
	"""
    with GroTrrReader(filename) as trrfile:
        for count, frame in enumerate(trrfile):
            data = trrfile.get_data()
    box = data["box"]
    positions = np.array(data['x'])

    atoms = []
    for count, position in enumerate(positions):
        atom = pca.Atom()
        atom.pos = list(position)
        atom.id = count+1
        atom.loc = count
        atoms.append(atom)

    return atoms, box


def write_snap(filename):
	raise NotImplementedError("write method for mdtraj is not implemented")

def split_snaps(infile, compressed = False):
    snaps = []
    if os.path.exists(infile):
        with GroTrrReader(infile) as trrfile:
            for count, frame in enumerate(trrfile):
                outfile = ".".join([infile, "frame", str(count), "trr"])
                snaps.append(outfile)
                data = trrfile.get_data()
                natoms = len(data["x"])
                data['natoms'] = natoms
                data['step'] = 0
                data['time'] = 0
                data['lambda'] = 0
                header = write_trr_frame(
                outfile,
                data,
                double=frame['double'],
                append=False,
                endian=frame['endian']
            )


    return snaps 

def convert_snap(**kwargs):
	raise NotImplementedError("convert method for mdtraj is not implemented")


def swap_integer(integer):
    """Convert little/big endian."""
    return (((integer << 24) & 0xff000000) | ((integer << 8) & 0x00ff0000) |
            ((integer >> 8) & 0x0000ff00) | ((integer >> 24) & 0x000000ff))


def swap_endian(endian):
    """Just swap the string for selecting big/little."""
    if endian == '>':
        return '<'
    elif endian == '<':
        return '>'
    else:
        raise ValueError('Undefined swap!')


def read_struct_buff(fileh, fmt):
    """Unpack from a filehandle with a given format.

    Parameters
    ----------
    fileh : file object
        The file handle to unpack from.
    fmt : string
        The format to use for unpacking.

    Returns
    -------
    out : tuple
        The unpacked elements according to the given format.

    Raises
    ------
    EOFError
        We will raise an EOFError if `fileh.read()` attempts to read
        past the end of the file.
    """
    buff = fileh.read(struct.calcsize(fmt))
    if not buff:
        raise EOFError
    else:
        return struct.unpack(fmt, buff)


def read_matrix(fileh, endian, double):
    """Read a matrix from the TRR file.

    Here, we assume that the matrix will be of
    dimensions (DIM, DIM).

    Parameters
    ----------
    fileh : file object
        The file handle to read from.
    endian : string
        Determines the byte order.
    double : boolean
        If true, we will assume that the numbers
        were stored in double precision.

    Returns
    -------
    mat : numpy.array
        The matrix as an numpy array.
    """
    if double:
        fmt = '{}{}d'.format(endian, DIM*DIM)
    else:
        fmt = '{}{}f'.format(endian, DIM*DIM)
    read = read_struct_buff(fileh, fmt)
    mat = np.zeros((DIM, DIM))
    for i in range(DIM):
        for j in range(DIM):
            mat[i, j] = read[i*DIM + j]
    return mat


def read_coord(fileh, endian, double, natoms):
    """Read a coordinate section from the TRR file.

    This method will read the full coordinate section from a TRR 
    file. The coordinate section may be positions, velocities or
    forces.

    Parameters
    ----------
    fileh : file object
        The file handle to read from.
    endian : string
        Determines the byte order.
    double : boolean
        If true, we will assume that the numbers
        were stored in double precision.
    natoms : int
        The number of atoms we have stored coordinates for.

    Returns
    -------
    mat : numpy.array
        The coordinates as a numpy array. It will have
        ``natoms`` rows and ``DIM`` columns.
    """
    if double:
        fmt = '{}{}d'.format(endian, natoms * DIM)
    else:
        fmt = '{}{}f'.format(endian, natoms * DIM)
    read = read_struct_buff(fileh, fmt)
    mat = np.array(read)
    mat.shape = (natoms, DIM)
    return mat


def is_double(header):
    """Determines we we should use double precision.

    This method determined the precision to use when reading
    the TRR file. This is based on the header read for a given
    frame which defines the sizes of certain "fields" like the box
    or the positions. From this size, the precision can be obtained.

    Parameters
    ----------
    header : dict
        The header read from the TRR file.

    Returns
    -------
    out : boolean
        True if we should use double precision.
    """
    key_order = ('box_size', 'x_size', 'v_size', 'f_size')
    size = 0
    for key in key_order:
        if header[key] != 0:
            if key == 'box_size':
                size = int(header[key] / DIM**2)
                break
            else:
                size = int(header[key] / (header['natoms'] * DIM))
                break
    if (size != SIZE_FLOAT) and (size != SIZE_DOUBLE):
        raise ValueError('Could not determine size!')
    else:
        return size == SIZE_DOUBLE


def read_trr_header(fileh):
    """Read a header from a TRR file.

    Parameters
    ----------
    fileh : file object
        The file handle for the file we are reading.

    Returns
    -------
    header : dict
        The header read from the file.
    """
    endian = '>'

    magic = read_struct_buff(fileh, '{}1i'.format(endian))[0]

    if magic == GROMACS_MAGIC:
        pass
    else:
        magic = swap_integer(magic)
        endian = swap_endian(endian)

    slen = read_struct_buff(fileh, '{}2i'.format(endian))
    raw = read_struct_buff(fileh, '{}{}s'.format(endian, slen[0]-1))
    version = raw[0].split(b'\0', 1)[0].decode('utf-8')
    if not version == TRR_VERSION:
        raise ValueError('Unknown format')

    head_fmt = HEAD_FMT.format(endian)
    head_s = read_struct_buff(fileh, head_fmt)
    header = {}
    for i, val in enumerate(head_s):
        key = HEAD_ITEMS[i]
        header[key] = val
    # The next are either floats or double
    double = is_double(header)
    if double:
        fmt = '{}2d'.format(endian)
    else:
        fmt = '{}2f'.format(endian)
    header_r = read_struct_buff(fileh, fmt)
    header['time'] = header_r[0]
    header['lambda'] = header_r[1]
    header['endian'] = endian
    header['double'] = double
    return header


def skip_trr_data(fileh, header):
    """Skip coordinates/box data etc.

    This method is used when we want to skip a data section in
    the TRR file. Rather than reading the data it will use the
    sized read in the header to skip ahead to the next frame.

    Parameters
    ----------
    fileh : file object
        The file handle for the file we are reading.
    header : dict
        The header read from the TRR file.
    """
    offset = sum([header[key] for key in DATA_ITEMS])
    fileh.seek(offset, 1)
    return None


def read_trr_data(fileh, header):
    """Read box, coordinates etc. from a TRR file.

    Parameters
    ----------
    fileh : file object
        The file handle for the file we are reading.
    header : dict
        The header read from the file.

    Returns
    -------
    data : dict
        The data we read from the file. It may contain the following
        keys if the data was found in the frame:

        - ``box`` : the box matrix,
        - ``vir`` : the virial matrix,
        - ``pres`` : the pressure matrix,
        - ``x`` : the coordinates,
        - ``v`` : the velocities, and
        - ``f`` : the forces
    """
    data = {}
    endian = header['endian']
    double = header['double']
    for key in ('box', 'vir', 'pres'):
        header_key = '{}_size'.format(key)
        if header[header_key] != 0:
            data[key] = read_matrix(fileh, endian, double)
    for key in ('x', 'v', 'f'):
        header_key = '{}_size'.format(key)
        if header[header_key] != 0:
            data[key] = read_coord(fileh, endian, double,
                                   header['natoms'])
    return data


def _write_trr_header(outfile, header, floatfmt, endian=None):
    """Helper method for writing a header to a TRR file.

    Parameters
    ----------
    outfile : filehandle
        The file we can write to.
    header : dict
        The header data for the TRR file.
    floatfmt : string
        The string which gives the format for floats. It should indicate
        if we are writing for double or single precision.
    endian : string, optional
        Can be used to force endianess.
    """
    slen = (13, 12)
    fmt = ['1i', '2i', '{}s'.format(slen[0] - 1), '13i']
    if endian:
        fmt = [endian + i for i in fmt]
    outfile.write(struct.pack(fmt[0], GROMACS_MAGIC))
    outfile.write(struct.pack(fmt[1], *slen))
    outfile.write(struct.pack(fmt[2], TRR_VERSION_B))
    head = [header[key] for key in HEAD_ITEMS[:13]]
    outfile.write(struct.pack(fmt[3], *head))
    outfile.write(struct.pack(floatfmt.format(1), header['time']))
    outfile.write(struct.pack(floatfmt.format(1), header['lambda']))


def write_trr_frame(filename, data, endian=None, double=False, append=False):
    """Write data in TRR format to a file.

    Parameters
    ----------
    filename : string
        The file we will write to.
    data : dict
        The data we will write to the file.
    endian : string, optional
        Select the byte order; big-endian or little-endian. If not
        specified, the native byte order will be used.
    double : boolean, optional
        If True, we will write in double precision.
    append : boolean, optional
        If True, we will append to the given file.
    """
    if double:
        size = SIZE_DOUBLE
        floatfmt = '{}d'
    else:
        size = SIZE_FLOAT
        floatfmt = '{}f'
    if endian:
        floatfmt = endian + floatfmt

    header = {}
    for key in HEAD_ITEMS:
        header[key] = 0

    header['natoms'] = data['natoms']
    header['step'] = data['step']
    header['box_size'] = size * DIM * DIM
    for i in ('x', 'v', 'f'):
        if i in data:
            header['{}_size'.format(i)] = data['natoms'] * size * DIM
    header['endian'] = endian
    header['double'] = double
    header['time'] = data['time']
    header['lambda'] = data['lambda']

    if append:
        mode = 'ab'
    else:
        mode = 'wb'
    with open(filename, mode) as outfile:
        _write_trr_header(outfile, header, floatfmt, endian=endian)
        for key in DATA_ITEMS:
            if header[key] != 0:
                # Note: We assume that the data is a numpy array, and that
                # we can find it as data['x'], data['v'], ... and so on.
                matrix = data[key.split('_')[0]]
                fmt = floatfmt.format(matrix.size)
                outfile.write(struct.pack(fmt, *matrix.flatten()))
    return header


class GroTrrReader():
    """A simple class for reading frames from a GROMACS TRR file.

    Attributes
    ----------
    filename : string
        The file name to open
    _skip : boolean
        If True, the next frame will be skipped. This is used
        to control the reading when iterating so that we do not
        read a full frame unless explicitly told to.
    fileh : file object
        The open file handle.
    header : dict
        The previously read header from the TRR file.
    """

    def __init__(self, filename):
        """Initiate the reader.

        Parameters
        ----------
        filename : string
            The name of the file to open.
        """
        self.filename = filename
        self._skip = False
        self.fileh = None
        self.header = None

    def __enter__(self):
        """Just open the file."""
        self.fileh = open(self.filename, 'rb')
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Ensure that we close the file."""
        self.fileh.close()

    def read_frame(self, read_data=True):
        """Read a new frame from the file.

        Parameters
        ---------
        read_data : boolean
            If False, we will not read the data, just the header and
            skip forward to the next header position.

        Returns
        -------
        out[0] : dict
            The header read from the file.
        out[1] : dict
            The data section read from the file.
        """
        header = read_trr_header(self.fileh)
        if read_data:
            data = read_trr_data(self.fileh, header)
        else:
            skip_trr_data(self.fileh, header)
            data = {}
        return header, data

    def __iter__(self):
        """Just set up the iterator."""
        return self

    def __next__(self):
        """Returns the next frame.

        Note
        ----
        In order to get the data for the current frame,
        ``self.get_data()`` needs to be called. The default
        behaviour is thus to skip the data in frames.
        """
        try:
            if self._skip:
                self.skip_data()
            header = read_trr_header(self.fileh)
            self._skip = True
            self.header = header
            return header
        except EOFError:
            raise StopIteration

    def get_data(self):
        """Get data from a frame and return it."""
        self._skip = False
        return read_trr_data(self.fileh, self.header)

    def skip_data(self):
        """Just skip data."""
        self._skip = False
        #print(self.header)
        skip_trr_data(self.fileh, self.header)
