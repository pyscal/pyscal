import numpy as np

from pyscal.attributes import DocumentedKeywords

class Repeat(DocumentedKeywords):
    def __init__(self, system, repetitions):
        _input = {
            "system": {
                "doc": "pyscal System", 
                "value": None,
            },      
            "repetitions": {
                "doc": "Number of repetitions in each direction",
                "value": None, 
            },
            "ghost": {
                "doc": "Make the generated atoms ghost", 
                "value": None,
            },
            "scale_box": {
                "doc": "If True, scale the box", 
                "value": None,
            },                              
            "atoms": {
                "doc": "Atoms contained in the system", 
                "value": None,
            },
            "return_atoms": {
                "doc": "If True, return the atoms instead of the system", 
                "value": None,
            },                                                              

        }

        self._create_input_tree(_input)
        self.
        self.input.system = system
        self.input.repetitions = repetitions

    def calculate(self):
        box = self.input.system.box        
        self.input.system.actual_box = box.copy()

        if self.input.atoms is None:
            atoms = self.input.system.atoms
        else:
            atoms = self.input.atoms

        newatoms = []
        idstart = len(atoms) + 1

        x1 = -self.input.repetitions[0]
        x2 = self.input.repetitions[0]+1
        y1 = -self.input.repetitions[1]
        y2 = self.input.repetitions[1]+1
        z1 = -self.input.repetitions[2]
        z2 = self.input.repetitions[2]+1
        xs = 2*self.input.repetitions[0] + 1
        ys = 2*self.input.repetitions[1] + 1
        zs = 2*self.input.repetitions[2] + 1
        
        datadict = {key:[] for key in atoms.keys()}
        del datadict['positions']
        del datadict['ids']
        del datadict['head']
        del datadict['ghost']
        positions = []
        ids = []
        head = []
        ghosts = []

        for i in range(x1, x2):
            for j in range(y1, y2):
                for k in range(z1, z2):
                    if (i==j==k==0):
                        continue
                    for count, pos in enumerate(atoms['positions']):
                        #we should create ghost images for only real atoms
                        if not atoms["ghost"][count]:
                            pos = (pos + i*np.array(box[0]) + j*np.array(box[1]) + k*np.array(box[2]))
                            positions.append(list(pos))
                            ids.append(idstart)
                            head.append(count)
                            ghosts.append(self.input.ghost)
                            idstart += 1
                            for key in datadict.keys():
                                datadict[key].append(atoms[key][count])

        if self.input.scale_box:
            box[0] = xs*np.array(box[0])
            box[1] = ys*np.array(box[1])
            box[2] = zs*np.array(box[2])
            self.input.system.box = box
        if self.input.ghost:
            self.input.system.ghosts_created = True

        atoms['positions'].extend(positions)
        atoms['ids'].extend(ids)
        atoms['ghost'].extend(ghosts)
        atoms['head'].extend(head)
        for key in datadict.keys():
            atoms[key].extend(datadict[key])

        if self.input.return_atoms:
            return atoms
        else:
            self.input.system.atoms = atoms
            return self.input.system

