from pyscal.crystal_structures import Structure
from pyscal.operations.symmetry import Symmetry
from pyscal.core import System

def test_bcc():
	struct = Structure()
	sys = struct.element.Fe()
	sym = Symmetry(sys)
	sym.calculate()
	assert sym.output.international_symbol == 'Im-3m'	

def test_fcc():
	struct = Structure()
	sys = struct.element.Cu()
	sym = Symmetry(sys)
	sym.calculate()
	assert sym.output.international_symbol == 'Fm-3m'	

def test_hcp():
	struct = Structure()
	sys = struct.element.Ti()
	sym = Symmetry(sys)
	sym.calculate()
	assert sym.output.international_symbol == 'P6_3/mmc'	

def test_diamond():
	struct = Structure()
	sys = struct.element.Si()
	sym = Symmetry(sys)
	sym.calculate()
	assert sym.output.international_symbol == 'Fd-3m'	