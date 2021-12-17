"""


"""


import os
import numpy as np
import warnings
import itertools
from ase.io import write
import uuid
import gzip
import io

import pyscal.csystem as pc
#import pyscal.traj_process as ptp
#from pyscal.formats.ase import convert_snap
#import pyscal.routines as routines
#import pyscal.visualization as pv


class System:
    """
    Python class for holding the properties of an atomic configuration 
    """
    def __init__(self):
        self.initialized = True
        self.neighbors_found = False
        self.neighbor_method = None
        self.ghosts_created = False
        self.actual_box = None

        #box parameters
        self.rot = [[0,0,0], [0,0,0], [0,0,0]]
        self.rotinv = [[0,0,0], [0,0,0], [0,0,0]]
        self.boxdims = [0,0,0]
        self.triclinic = 0
        self._atoms = {}

    @property
    def natoms(self):
        if 'positions' in self.atoms.keys():
            nop = np.sum([1 for x in range(len(self.atoms["positions"])) if self.atoms["ghost"][x]==False])
            return nop
        else:
            return 0
    
    @property
    def box(self):
        """
        Wrap for inbuilt box
        """
        if self.actual_box is not None:
            return self.actual_box
        else:
            return self._box

    @box.setter
    def box(self, userbox):
        """
        Box setter 
        """
        #we should automatically check for triclinic cells here
        summ = 0
        for i in range(3):
            box1 = np.array(userbox[i-1])
            box2 = np.array(userbox[i])
            summ += np.dot(box1, box2)/(np.linalg.norm(box1)*np.linalg.norm(box2))

        #check if the summ is zero
        if (np.abs(summ) > 1E-6):
            #this is a triclinic box
            rot = np.array(userbox).T
            rotinv = np.linalg.inv(rot)
            self.triclinic = 1
            self.rot = rot
            self.rotinv = rotinv

        self.boxdims[0] = np.sum(np.array(userbox[0])**2)**0.5
        self.boxdims[1] = np.sum(np.array(userbox[1])**2)**0.5
        self.boxdims[2] = np.sum(np.array(userbox[2])**2)**0.5
        self._box = userbox

    @property
    def atoms(self):
        """
        Atom access
        """
        return self._atoms

    @atoms.setter
    def atoms(self, atoms):
        """
        Set atoms
        """
        #we need to check atoms and add necessary keys
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        for key, val in atoms.items():
            if not (len(val)==nop):
                raise ValueError("All times in the atoms dict should have same length as positions")
        #now add necessary keys-ids, types, ghost
        if not 'ids' in atoms.keys():
            atoms['ids'] = [x+1 for x in range(nop)]
        if not 'types' in atoms.keys():
            atoms['types'] = [1 for x in range(nop)]
        if not 'ghost' in atoms.keys():
            atoms['ghost'] = [False for x in range(nop)]

        if(len(atoms['positions']) < 200):
            #we need to estimate a rough idea
            needed_atoms = 200 - len(atoms)
            #get a rough cell
            needed_cells = np.ceil(needed_atoms/len(atoms))
            nx = int(needed_cells**(1/3))
            nx = int(np.ceil(nx/2))

            if np.sum(self.box) == 0:
                raise ValueError("Simulation box should be initialized before atoms")
            atoms = self.repeat((nx, nx, nx), atoms=atoms, ghost=True, scale_box=True)

        self._atoms = atoms

    #iterator for atoms
    def iter_atoms(self):
        """
        Iter over atoms
        """
        for i in range(self.natoms):
            rdict = {}
            for key in self.atoms.keys():
                rdict[key] = self.atoms[key][i]
            yield rdict 

    def add_atoms(self, atoms):
        """
        Cleanly add a given list of atoms

        Parameters
        ----------
        atoms : dict

        Returns
        -------
        None
        """ 
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        for key, val in atoms.items():
            if not (len(val)==nop):
                raise ValueError("All times in the atoms dict should have same length as positions")
        
        #now add necessary keys-ids, types, ghost
        maxid = max(self.atoms["ids"])
        if not 'ids' in atoms.keys():
            atoms['ids'] = [maxid+x+1 for x in range(nop)]
        else:
            for i in atoms['ids']:
                if i in self.atoms['ids']:
                    raise ValueError("Atom id already exists, unique ID is required")

        if not 'types' in atoms.keys():
            atoms['types'] = [1 for x in range(nop)]
        if not 'ghost' in atoms.keys():
            atoms['ghost'] = [False for x in range(nop)]

        for key in self.atoms.keys():
            self.atoms[key] = [*self.atoms[key], *atoms[key]]

    def repeat(self, reps, atoms=None, ghost=False, scale_box=True):
        """
        Replicate simulation cell
        
        Parameters
        ----------
        reps : list of ints of size 3
            repetitions in each direction
        
        atoms : list of atoms, optional
            if not provided, use atoms that are assigned

        ghost : bool, optional
            If True, assign the new atoms as ghost instead of actual atoms
        """
        
        box = self.box        
        self.actual_box = box.copy()

        if atoms is None:
            atoms = self.atoms

        newatoms = []
        idstart = len(atoms) + 1

        x1 = -reps[0]
        x2 = reps[0]+1
        y1 = -reps[1]
        y2 = reps[1]+1
        z1 = -reps[2]
        z2 = reps[2]+1
        xs = 2*reps[0] + 1
        ys = 2*reps[1] + 1
        zs = 2*reps[2] + 1
        tp = atoms['types'][0]

        positions = []
        ids = []
        types = []
        ghosts = []

        for i in range(x1, x2):
            for j in range(y1, y2):
                for k in range(z1, z2):
                    if (i==j==k==0):
                        continue
                    for pos in atoms['positions']:
                        pos = (pos + i*np.array(box[0]) + j*np.array(box[1]) + k*np.array(box[2]))
                        positions.append(pos)
                        ids.append(idstart)
                        idstart += 1
                        types.append(tp)
                        ghosts.append(ghost)

        if scale_box:
            box[0] = xs*np.array(box[0])
            box[1] = ys*np.array(box[1])
            box[2] = zs*np.array(box[2])
            self.box = box
        if ghost:
            self.ghosts_created = True

        atoms['positions'] = [*atoms['positions'], *positions]
        atoms['ids'] = [*atoms['ids'], *ids]
        atoms['types'] = [*atoms['types'], *types]
        atoms['ghost'] = [*atoms['ghost'], *ghosts]

        return atoms


    def embed_in_cubic_box(self,):
        """
        Embedded the triclinic box in a cubic box
        """
        #first task is to create a box representation
        
        box = self._box
        backupbox = box.copy()
        
        a = np.array(box[0])
        b = np.array(box[1])
        c = np.array(box[2])

        cosa = np.dot(b, c)/(np.linalg.norm(b)*np.linalg.norm(c))
        cosb = np.dot(c, a)/(np.linalg.norm(c)*np.linalg.norm(a))
        cosc = np.dot(a, b)/(np.linalg.norm(a)*np.linalg.norm(b))

        lx = np.linalg.norm(a)
        xy = np.linalg.norm(b)*cosc
        xz = np.linalg.norm(c)*cosb
        ly = np.sqrt(np.linalg.norm(b)*np.linalg.norm(b) - xy*xy)
        yz = (np.linalg.norm(b)*np.linalg.norm(c)*cosa - xy*xz)/ly
        lz = np.sqrt(np.linalg.norm(c)*np.linalg.norm(c) - xz*xz - yz*yz)

        xlo = ylo = zlo = 0
        xhi = lx
        yhi = ly
        zhi = lz

        xlo_bound = xlo + min(0.0,xy,xz,xy+xz)
        xhi_bound = xhi + max(0.0,xy,xz,xy+xz)
        ylo_bound = ylo + min(0.0,yz)
        yhi_bound = yhi + max(0.0,yz)
        zlo_bound = zlo
        zhi_bound = zhi

        newbox = np.array([[xhi_bound-xlo_bound, 0, 0], [0, yhi_bound-ylo_bound, 0], [0, 0, zhi_bound-zlo_bound]])
        
        self.newbox = newbox
        self.box = newbox
        self.actual_box = None

    """
    def extract_cubic_box(self, repeat=(3,3,3)):
        """
        #Extract a cubic representation of the box from triclinic cell

        #Parameters
        #----------
        #repeat : list of ints
        #    the number of times box should be repeat

        #Returns
        #-------
        #cubebox: list of list of floats
        #    cubic box

        #atoms: list of Atom objects
        #    atoms in the cubic box
        """
        ratoms = self.repeat(repeat)
        anchor = ratoms["positions"][0]

        nb = []

        for count, ind in enumerate([[1, 2], [0, 2], [0, 1]]): 
            mindist = 10000
            minpos = None
            for i, pos in enumerate(ratoms["positions"][1:]):
                if np.abs(pos[ind[0]]-anchor[ind[0]])< 1E-5:
                    if np.abs(pos[ind[1]]-anchor[ind[1]])< 1E-5:
                        if ratoms["types"][i] == ratoms["types"][0]:
                            dist = self.get_distance(anchor, pos)
                            print(dist)
                            if dist < mindist:
                                mindist = dist
                                minpos = pos
            if minpos is None:
                raise ValueError("Cubic box cound not be found, try increasing repeat")
            lo = min(anchor[count], minpos[count])
            hi = max(anchor[count], minpos[count])
            nb.append([lo, hi])

        print(mindist)
        #collect atoms
        collectedatoms = []
        for i, pos in enumerate(ratoms["positions"]):
            if (nb[0][0] <= pos[0] < nb[0][1]):
                if (nb[1][0] <= pos[1] < nb[1][1]):
                    if (nb[2][0] <= pos[2] < nb[2][1]):
                        collectedatoms.append(i)

        #remap atoms to box
        newdict = {}
        for key in ratoms.keys():
            newdict[key] = []

        for i in collectedatoms:
            pos = [ratoms["positions"][i][0]-nb[0][0], ratoms["positions"][i][1]-nb[1][0], ratoms["positions"][i][2]-nb[2][0]]
            newdict["positions"].append(pos)
            for key in ratoms.keys():
                newdict[key].append(ratoms[key][i])

        #create a new box
        cubebox = [[(nb[0][1]-nb[0][0]), 0, 0],
                       [0, (nb[1][1]-nb[1][0]), 0],
                       [0, 0, (nb[2][1]-nb[2][0])]]

        return cubebox, newdict
    """

    def get_distance(self, pos1, pos2):
        """
        Get the distance between two atoms.

        Parameters
        ----------
        pos1 : list
                first atom position
        pos2 : list
                second atom position

        Returns
        -------
        distance : double
                distance between the first and second atom.

        Notes
        -----
        Periodic boundary conditions are assumed by default.
        """
        return pc.get_abs_distance(pos1, pos2, self.triclinic, 
            self.rot, self.rotinv, self.boxdims, 0.0, 0.0, 0.0)