import numpy as np
from pyscal.attributes import DocumentedKeywords
import spglib

class Symmetry(DocumentedKeywords):
    def __init__(self, system, angle_tolerance=-1.0,
        symprec=1e-5, hall_number=0):
        _input = {
            "system": {
                "doc": "pyscal System", 
                "value": None,
            },      
            "angle_tolerance": {
                "doc": "An experimental argument that controls angle tolerance between basis vectors. Normally it is not recommended to use this argument.",
                "value": -1.0,
                "url": "https://spglib.github.io/spglib/variable.html#angle-tolerance", 
            },
            "symprec": {
                "doc": "Distance tolerance in Cartesian coordinates to find crystal symmetry.", 
                "value": 1e-5,
                "url": "https://spglib.github.io/spglib/variable.html#symprec"
            },
            "hall_number": {
                "doc": "The argument to constrain the space-group-type search only for the Hall symbol corresponding to it.", 
                "value": 0,
                "url": "https://spglib.github.io/spglib/dataset.html#hall-number"
            },                              
        }

        self._create_input_tree(_input)
        self.input.system = system
        self.input.angle_tolerance = angle_tolerance
        self.input.symprec = symprec
        self.input.hall_number = hall_number

    def calculate(self):
        results = spglib.get_symmetry_dataset((self.input.system.box,
            self.input.system.direct_coordinates, self.input.system.atoms.types))

        #now create outputs
        _output = {    
            "international_space_group_number": {
                "doc": "The space group type number defined in International Tables for Crystallography (ITA).",
                "value": None,
                "url": "https://it.iucr.org/Ac/contents/", 
            },
            "international_symbol": {
                "doc": "The (full) Hermann–Mauguin notation of space group type", 
                "value": None,
                "url": "https://www.chemeurope.com/en/encyclopedia/Hermann-Mauguin_notation.html"
            },
            "point_group": {
                "doc": "Symbol of the crystallographic point group in the Hermann–Mauguin notation", 
                "value": None,
                "url": "https://www.chemeurope.com/en/encyclopedia/Hermann-Mauguin_notation.html"
            },                              
        }

        self._create_input_tree(_output, name="output")
        self.output.international_space_group_number = results["number"]
        self.output.international_symbol = results["international"]
        self.output.point_group = results["pointgroup"]

        #this is temp, for checking
        self.results = results

