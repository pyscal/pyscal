import pytest

try:
	import steinhardt as st
except ImportError:
	print "Steinhardt was not imported correctly"

def test_get_absDistance():
	atom1 = st.Atom([3,4,0],[0,0,0],1)
	atom2 = st.Atom([0,0,0],[0,0,0],2)
	box1 = st.Simbox([10,10,10])
	dist = st.get_absDistance(atom1,atom2,box1)
	assert dist == 5

def test_convert_SphericalCoordinates():
        atom1 = st.Atom([3,4,0],[0,0,0],1)
        atom2 = st.Atom([0,0,0],[0,0,0],2)
        box1 = st.Simbox([10,10,10])
        diff = st.get_diff(atom1,atom2,box1)
        res = st.convert_SphericalCoordinates(diff)
        assert res[0]==5
        assert res[1]== 1.5707963267948966
        assert res[2]== 0.9272952180016122

def test_get_AllNeighborsandDistances():
        
        a1 = st.Atom([0,0,0],[0,0,0],1)
        a2 = st.Atom([1,0,0],[0,0,0],2)
        a3 = st.Atom([0,1,0],[0,0,0],3)
        a4 = st.Atom([0,0,1],[0,0,0],4)
        a5 = st.Atom([1,1,0],[0,0,0],5)
        a6 = st.Atom([1,0,1],[0,0,0],6)
        a7 = st.Atom([0,1,1],[0,0,0],7)
        a8 = st.Atom([1,1,1],[0,0,0],8)
        a9 = st.Atom([2,1,1],[0,0,0],9)
        a10 = st.Atom([1,2,1],[0,0,0],10)
        a11 = st.Atom([1,1,2],[0,0,0],11)
        a12 = st.Atom([2,2,1],[0,0,0],12)
        a13 = st.Atom([2,1,2],[0,0,0],13)
        a14 = st.Atom([1,2,2],[0,0,0],14)
        box1 = st.Simbox([2,2,2])
        param = st.Params(1.2,[4])
        atoms = [a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14]
        atoms2 = st.get_AllNeighborsandDistances(atoms,box1,param)
        assert atoms2[0].gneighborcount() == 6


def unpack_snapshot(filename):
        atoms = []
        with open(filename) as infile:    
                count = 1
                timestep = 0
                for line in infile:
                        if count==6:
                                box = []
                                raw = line.strip().split()
                                dimx = float(raw[1])-float(raw[0])
                
                        elif count==7: 
                                raw = line.strip().split()
                                dimy = float(raw[1])-float(raw[0])
                
                        elif count==8:
                                raw = line.strip().split()
                                dimz = float(raw[1])-float(raw[0])
                
                        elif count>9:
                                line = line.strip().split()
                                idd = int(line[0]) 
                                x = float(line[3])
                                y = float(line[4])
                                z = float(line[5])
                                vx = float(line[6])
                                vy = float(line[7])
                                vz = float(line[8])
                                xx = [x,y,z]
                                vv = [vx,vy,vz]
                                a = st.Atom(xx,vv,idd)
                                atoms.append(a)
                        count+=1
                
        box = st.Simbox([dimx,dimy,dimz])
        return atoms,box


def  test_module_aq6():
        atoms,box = unpack_snapshot("conf.dump")
        param = st.Params(3.63,[6])
        atoms = st.get_AllNeighborsandDistances(atoms,box,param)
        atoms = st.calculate_complexQLM(atoms,box,param)
        atoms = st.calculate_aQ(atoms,box,param)

        assert atoms[0].gaq()[0] == 0.3127369144275015
        assert atoms[30].gaq()[0] == 0.4107729464466044


def  test_module_nucsize():
        atoms,box = unpack_snapshot("conf.dump")
        param = st.Params(3.63,[6])
        param.setnucsizeparams(0.5,0.5,7)
        atoms = st.get_AllNeighborsandDistances(atoms,box,param)
        atoms = st.calculate_complexQLM(atoms,box,param)
        atoms = st.calculate_aQ(atoms,box,param)
        atoms = st.calculate_aQ(atoms,box,param)
        atoms = st.find_solids(atoms,param)
        atoms = st.find_clusters(atoms)
        nuc = st.largest_cluster(atoms)
        
        assert nuc>=60 and nuc<=70