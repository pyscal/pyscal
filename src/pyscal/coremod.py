import pyscal.traj_process as ptp
from pyscal.formats.ase import convert_snap
import pyscal.routines as routines
import os
import numpy as np
import warnings
import itertools
from ase.io import write
import uuid
import gzip
import io

import pyscal.visualization as pv
import pyscal.cmodsystem as pc

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
        self._atoms = None
 	
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
        if(len(atoms) < 200):
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


        for i in range(x1, x2):
            for j in range(y1, y2):
                for k in range(z1, z2):
                    if (i==j==k==0):
                        continue
                    for atom in atoms:
                        pos = np.array(atom['pos'])
                        pos = (pos + i*np.array(box[0]) + j*np.array(box[1]) + k*np.array(box[2]))
                        a = {}
                        a['pos'] = pos
                        a['id'] = idstart
                        idstart += 1
                        a['type'] = atom['type']
                        a['ghost'] = ghost
                        newatoms.append(a)

        if scale_box:
            box[0] = xs*np.array(box[0])
            box[1] = ys*np.array(box[1])
            box[2] = zs*np.array(box[2])
            self.box = box
        if ghost:
            self.ghosts_created = True

        completeatoms = atoms + newatoms
        #print(len(completeatoms))
        return completeatoms