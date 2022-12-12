"""
This file includes the definition of an Atoms class which can be used with
System

TODO
----
- Iterators for masked/unmasked atoms
- Iterators for selected/unselected atoms
"""

import numpy as np
import warnings
from pyscal.attributes import generate_class

AttrSetter = generate_class()

class Atoms(dict, AttrSetter):
    def __init__(self, *args, **kwargs):
        self.update(*args, **kwargs)
        self._nreal = 0
        self._nghost = 0
        
    
    def __getitem__(self, key):
        if isinstance(key, slice):
            return self._get_atoms(key)
        elif isinstance(key, int):
            return self._get_atoms(key)
        else:
            val = dict.__getitem__(self, key)
            return val

    def __setitem__(self, key, val):
        dict.__setitem__(self, key, val)

    def __repr__(self):
        dictrepr = dict.__repr__(self)
        return '%s(%s)' % (type(self).__name__, dictrepr)

    def _repr_json(self):
        #convert to atom base dict
        disp_atoms = {f"atom {x}": self._get_atoms(x) for x in range(self.natoms)}
        return disp_atoms

        
    def update(self, *args, **kwargs):
        for k, v in dict(*args, **kwargs).items():
            self[k] = v

    def __add__(self, atoms):
        if not isinstance(atoms, Atoms):
            raise TypeError("Only Atoms can be added to Atoms")
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        val_length_check = np.prod([len(val)==nop for key, val in atoms.items()])
        if not val_length_check:
            raise ValueError("All times in the atoms dict should have same length as positions")
        
        #now add necessary keys-ids, types, ghost
        maxid = max(self["ids"])
        if not 'ids' in atoms.keys():
            atoms['ids'] = [maxid+x+1 for x in range(nop)]
        else:
            if len(set(atoms['ids']).intersection(set(self['ids']))):
                raise ValueError("Atom id already exists, unique ID is required")
        
        atoms['ghost'] = [False for x in range(nop)]
        if not 'types' in atoms.keys():
            atoms['types'] = [1 for x in range(nop)]
        if not 'mask_1' in atoms.keys():
            atoms['mask_1'] = [False for x in range(nop)]
        if not 'mask_2' in atoms.keys():
            atoms['mask_2'] = [False for x in range(nop)]
        if not 'condition' in atoms.keys():
            atoms['condition'] = [True for x in range(nop)]
        if not 'head' in atoms.keys():
            atoms['head'] = [self.natoms+x for x in range(nop)]
        
        common_keys = list(set(self.keys()).intersection(set(atoms.keys())))
        _ = [self[key].extend(atoms[key]) for key in common_keys]
        extra_keys_add = len(atoms.keys()) - len(common_keys)
        extra_keys_exist = len(self.keys()) - len(common_keys)
        
        if extra_keys_add > 0:
            warnings.warn("Some keys present in the atoms are add are not present in existing atoms, they were ignored")
        if extra_keys_exist > 0:
            warnings.warn("Some keys present in the existing Atoms were not present in the atoms to add, please recalculate")
        self._nreal += nop
        return self


    def __radd__(self, atoms):
        """
        Reverse add method
        """
        if ntraj == 0:
            return self
        else:
            return self.__add__(atoms)
    
    @property
    def natoms(self):
        return self._nreal

    @property
    def nreal(self):
        return self._nreal
    
    @property
    def nghost(self):
        return self._nghost
    
    @property
    def ntotal(self):
        return self._nreal + self._nghost

    def from_dict(self, atoms):
        if not 'positions' in atoms.keys():
            raise ValueError('positions is a necessary key in atoms')
        nop = len(atoms["positions"])
        
        if not 'ids' in atoms.keys():
            atoms['ids'] = [x+1 for x in range(nop)]
        
        atoms['ghost'] = [False for x in range(nop)]
        if not 'types' in atoms.keys():
            atoms['types'] = [1 for x in range(nop)]
        if not 'mask_1' in atoms.keys():
            atoms['mask_1'] = [False for x in range(nop)]
        if not 'mask_2' in atoms.keys():
            atoms['mask_2'] = [False for x in range(nop)]
        if not 'condition' in atoms.keys():
            atoms['condition'] = [True for x in range(nop)]
        if not 'head' in atoms.keys():
            atoms['head'] = [x for x in range(nop)]
        
        for key, val in atoms.items():
            self[key] = val
        self._nreal = len(val)

        #add attributes
        mapdict = {"positions": "positions",
        "ids": "ids",
        "types": "types",
        "species": "types",
        "mask": {"primary": "mask_1", "secondary": "mask_2"},
        "selection": "condition",
        "condition": "condition",
        "head": "head"}

        self._add_attribute(mapdict)

        
    def _get_atoms(self, index):
        atom_dict = {key: self[key][index] for key in self.keys()}
        return atom_dict

    def _delete_atoms(self, indices):
        del_real = np.sum([1 for x in indices if x < self._nreal])
        del_ghost = np.sum([1 for x in indices if x >= self._nreal])

        for key in self.keys():
            for index in sorted(indices, reverse=True):
                del self[key][index]

        td = len(indices)
        self._nreal = int(self.nreal - del_real)
        self._nghost = int(self.nghost - del_ghost)


    def iter_atoms(self):
        for index in range(self.nreal):
            atom_dict = {key: self[key][index] for key in self.keys()}
            yield atom_dict

    def _generate_bool_list(self, ids=None, indices=None, condition=None):
        #necessary checks
        non_nones = sum(x is not None for x in [ids, indices, condition])
        if non_nones > 1:
            raise ValueError("Only one of ids, indices or condition should be provided")
        elif non_nones == 0:
            warnings.warn("No conditions provided, all atoms will be included")
        #generate a list of indices
        if ids is not None:
            if not isinstance(ids, list):
                ids = [ids]
            indices = [x for x in range(len(self["ids"])) if self["ids"][x] in ids]
            if len(indices) == 0:
                raise ValueError("No ids found to delete")
            if len(indices) != len(ids):
                warnings.warn("Not all ids were found")
        elif condition is not None:
            indices = [x for x in range(self.nreal) if condition(self._get_natoms(x))]
        elif indices is None:
            indices = [x for x in range(self.nreal)]
        
        if not isinstance(indices, list):
            indices = [indices]

        bool_list = [ True if x in indices else False for x in range(self.nreal)]
        return bool_list

    def _apply_mask(self, masks, mask_type):
        if (mask_type == 'primary') or (mask_type == 'all'):
            for i in range(self.ntotal):
                self["mask_1"][i] = masks[self["head"][i]]
        if (mask_type == 'secondary') or (mask_type == 'all'):
            for i in range(self.ntotal):
                self["mask_2"][i] = masks[self["head"][i]]

    def mask(self, mask_type="primary", ids=None, indices=None, condition=None):
        masks = self._generate_bool_list(ids=ids, indices=indices, condition=condition)
        self._apply_mask(masks, mask_type)


    def unmask(self, mask="primary", ids=None, indices=None, condition=None):
        masks = self._generate_bool_list(ids=ids, indices=indices, condition=condition)
        masks = [not x for x in masks]
        self._apply_mask(masks, mask_type)

    def _apply_selection(self, condition):
        for i in range(self.ntotal):
            self["condition"][i] = condition[self["head"][i]]
        
    def select(self, ids=None, indices=None, condition=None):
        masks = self._generate_bool_list(ids=ids, indices=indices, condition=condition)
        self._apply_selection(masks)
    
    def unselect(self, ids=None, indices=None, condition=None):
        masks = self._generate_bool_list(ids=ids, indices=indices, condition=condition)
        masks = [not x for x in masks]
        self._apply_selection(masks)

    def delete(self, ids=None, indices=None, condition=None):
        #delete atoms
        #reassign ids
        #reassign indices
        #reassign heads
        masks = self._generate_bool_list(ids=ids, indices=indices, condition=condition)
        delete_list = [masks[self["head"][x]] for x in range(self.ntotal)]
        delete_ids = [x for x in range(self.ntotal) if masks[x]]
        self._delete_atoms(delete_ids)

    @property
    def composition(self):
        typelist = self["types"][:self.nreal]
        types, typecounts = np.unique(typelist, return_counts=True)
        concdict = {str(t): typecounts[c] for c, t in enumerate(types)}
        return concdict

    