import numpy as np
from pyscal.attributes import DocumentedKeywords

def repeat(system, repetitions, ghost = False,
    scale_box = True, atoms = None, return_atoms = False):
    """
    Repeat the given system

    Parameters
    ----------
    system: pyscal System object
        the input system to be repeated

    repetitions: tuple of ints
        number of times the system is to be rotated in each direction

    ghost: bool, optional
        if True, make the new atoms ghost, default False

    scale_box: bool, optional
        if True, scale the simulation box, default True

    atoms: None, optional
        if provided use the given atoms, and not the atoms from the system

    return_atoms: bool, optional
        if True, return atoms instead of adding them to the system. Default False
    
    Returns
    -------
    system: pyscal System object
        the system with repetitions. Only returned if `return_atoms` is False.
    
    atoms: Atoms object
        only returned if `return_atoms` is True.


    """
    box = system.box        
    system.actual_box = box.copy()

    if atoms is None:
        atoms = system.atoms

    newatoms = []
    idstart = len(atoms) + 1

    x1 = -repetitions[0]
    x2 = repetitions[0]+1
    y1 = -repetitions[1]
    y2 = repetitions[1]+1
    z1 = -repetitions[2]
    z2 = repetitions[2]+1
    xs = 2*repetitions[0] + 1
    ys = 2*repetitions[1] + 1
    zs = 2*repetitions[2] + 1
    
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
                        ghosts.append(ghost)
                        idstart += 1
                        for key in datadict.keys():
                            datadict[key].append(atoms[key][count])

    if scale_box:
        box[0] = xs*np.array(box[0])
        box[1] = ys*np.array(box[1])
        box[2] = zs*np.array(box[2])
        system.box = box
    if ghost:
        system.ghosts_created = True

    atoms['positions'].extend(positions)
    atoms['ids'].extend(ids)
    atoms['ghost'].extend(ghosts)
    atoms['head'].extend(head)
    for key in datadict.keys():
        atoms[key].extend(datadict[key])

    if return_atoms:
        return atoms
    else:
        system.atoms = atoms
        return system


def embed_in_cubic_box(system, input_box=None, 
    return_box=False):
    """
    Embedded the triclinic box in a cubic box
    
    Parameters
    ----------
    system : pyscal System
        input system

    input_box: optional
        if specified, this box is used instead of the
        system simulation box

    return_box: bool, optional
        if specified, return the box instead of assigning
        to the system. Default False

    Returns
    -------
    system: 
        system object with the modified box
        only returned if `return_box` is False

    box:
        the cubic simulation box
        only returned if `return_box` is True
    """
    if input_box is None:
        box = system._box
        backupbox = box.copy()
    else:
        box = input_box
        backupbox = input_box.copy

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
    
    if not return_box:
        system.newbox = newbox
        system.box = newbox
        system.box_backup = backupbox
        system.actual_box = None
        return system
    else:
        return newbox
