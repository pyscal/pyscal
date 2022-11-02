"""
This file includes the definition of an Atoms class which can be used with
System
"""

class AttrSetter:
    def _add_attribute(self, mapdict, head=None):
        """
        Add mapping of attributes 
        """
        if head is None:
            head = self
        for key, val in mapdict.items():
            #here nesting is needed
            if isinstance(val, dict):
                #we have to create a nested class
                setattr(self, key, AttrSetter())
                getattr(self, key)._add_attribute(val, head=head)
            elif val in head.keys():
                setattr(self, key, self._filter_ghost(head, val))
    
    def _filter_ghost(self, head, val):
        return head[val][:head.natoms]
                
                
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
        
    def _get_atoms(self, index):
        atom_dict = {}
        for key in self.keys():
            atom_dict[key] = self[key][index]
        return atom_dict
