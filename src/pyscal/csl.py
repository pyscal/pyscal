import sys
import random
from math import degrees, atan, sqrt, pi, ceil, cos, acos, sin, gcd, radians
import numpy as np
from numpy import dot, cross
from numpy.linalg import det, norm, inv

def get_cubic_sigma(uvw, m, n=1):
    """
    CSL analytical formula based on the book:
    'Interfaces in crystalline materials',
     Sutton and Balluffi, clarendon press, 1996.
     generates possible sigma values.
    arguments:
    uvw -- the axis
    m,n -- two integers (n by default 1)

    """
    u, v, w = uvw
    sqsum = u*u + v*v + w*w
    sigma = m*m + n*n * sqsum
    while sigma != 0 and sigma % 2 == 0:
        sigma /= 2
    return sigma if sigma > 1 else None


def get_cubic_theta(uvw, m, n=1):
    """
    generates possible theta values.
    arguments:
    uvw -- the axis
    m,n -- two integers (n by default 1)
    """
    u, v, w = uvw
    sqsum = u*u + v*v + w*w
    if m > 0:
        return 2 * atan(sqrt(sqsum) * n / m)
    else:
        return pi

def get_theta_m_n_list(uvw, sigma):
    """
    Finds integers m and n lists that match the input sigma.
    """
    if sigma == 1:
        return [(0., 0., 0.)]
    thetas = []
    max_m = int(ceil(sqrt(4*sigma)))

    for m in range(1, max_m):
        for n in range(1, max_m):
            if gcd(m, n) == 1:
                s = get_cubic_sigma(uvw, m, n)
            if s == sigma:
                theta = (get_cubic_theta(uvw, m, n))
                thetas.append((theta, m, n))
                thetas.sort(key=lambda x: x[0])
    return thetas

def get_sigma_list(uvw, limit):
    """
    prints a list of smallest sigmas/angles for a given axis(uvw).
    """
    sigmas = []
    thetas = []
    for i in range(limit):
        tt = get_theta_m_n_list(uvw, i)
        if len(tt) > 0:
            theta, _, _ = tt[0]
            sigmas.append(i)
            thetas.append(degrees(theta))
    return sigmas, thetas

def get_theta_m_n_list(uvw, sigma):
    """
    Finds integers m and n lists that match the input sigma.
    """
    if sigma == 1:
        return [(0., 0., 0.)]
    thetas = []
    max_m = int(ceil(sqrt(4*sigma)))

    for m in range(1, max_m):
        for n in range(1, max_m):
            if gcd(m, n) == 1:
                s = get_cubic_sigma(uvw, m, n)
            if s == sigma:
                theta = (get_cubic_theta(uvw, m, n))
                thetas.append((theta, m, n))
                thetas.sort(key=lambda x: x[0])
    return thetas

def rot(a, Theta):
    """
    produces a rotation matrix.
    arguments:
    a -- an axis
    Theta -- an angle
    """
    c = cos(Theta)
    s = sin(Theta)
    a = a / norm(a)
    ax, ay, az = a
    return np.array([[c + ax * ax * (1 - c), ax * ay * (1 - c) - az * s,
                      ax * az * (1 - c) + ay * s],
                    [ay * ax * (1 - c) + az * s, c + ay * ay * (1 - c),
                        ay * az * (1 - c) - ax * s],
                     [az * ax * (1 - c) - ay * s, az * ay * (1 - c) + ax * s,
                      c + az * az * (1 - c)]])


def create_minimal_cell_method_1(sigma, uvw, R):
    """
    finds Minimal cell by means of a numerical search.
    (An alternative analytical method can be used too).
    arguments:
    sigma -- gb sigma
    uvw -- rotation axis
    R -- rotation matrix
    """
    uvw = np.array(uvw)
    MiniCell_1 = np.zeros([3, 3])
    MiniCell_1[:, 2] = uvw

    lim = 20
    x = np.arange(-lim, lim + 1, 1)
    y = x
    z = x
    indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(len(x) ** 3, 3)

    # remove 0 vectors and uvw from the list
    indice_0 = indice[np.where(np.sum(abs(indice), axis=1) != 0)]
    condition1 = ((abs(dot(indice_0, uvw) / norm(indice_0, axis=1) /
                       norm(uvw))).round(7))
    indice_0 = indice_0[np.where(condition1 != 1)]

    if minicell_search(indice_0, MiniCell_1, R, sigma):

        M1, M2 = minicell_search(indice_0, MiniCell_1, R, sigma)
        return (M1, M2)
    else:
        return None


def minicell_search(indices, MiniCell_1, R, sigma):

    tol = 0.001
    # norm1 = norm(indices, axis=1)
    newindices = dot(R, indices.T).T
    nn = indices[np.all(abs(np.round(newindices) - newindices) < 1e-6, axis=1)]
    TestVecs = nn[np.argsort(norm(nn, axis=1))]
    # print(len(indices), len(TestVecs),TestVecs[:20])

    Found = False
    count = 0
    while (not Found) and count < len(TestVecs) - 1:
        if 1 - ang(TestVecs[count], MiniCell_1[:, 2]) > tol:
            # and  (ang(TestVecs[i],uvw) > tol):
            MiniCell_1[:, 1] = (TestVecs[count])
            count += 1
            for j in range(len(TestVecs)):
                if (1 - ang(TestVecs[j], MiniCell_1[:, 2]) > tol) and (
                        1 - ang(TestVecs[j], MiniCell_1[:, 1]) > tol):
                    if (ang(TestVecs[j],
                            cross(MiniCell_1[:, 2], MiniCell_1[:, 1])) > tol):
                        #  The condition that the third vector can not be
                        #  normal to any other two.
                        #  and (ang(TestVecs[i],uvw)> tol) and
                        # (ang(TestVecs[i],MiniCell[:,1])> tol)):
                        MiniCell_1[:, 0] = (TestVecs[j]).astype(int)
                        Det1 = abs(round(det(MiniCell_1), 5))
                        MiniCell_1 = (MiniCell_1).astype(int)
                        MiniCell_2 = ((np.round(dot(R, MiniCell_1), 7))
                                      .astype(int))
                        Det2 = abs(round(det(MiniCell_2), 5))

                        if ((abs(Det1 - sigma)) < tol and
                                (abs(Det2 - sigma)) < tol):
                            Found = True
                            break
    if Found:
        return MiniCell_1, MiniCell_2
    else:
        return Found

def create_possible_gb_plane_list(uvw, m=5, n=1, lim=5):
    """
    generates GB planes and specifies the character.

    arguments:
    uvw -- axis of rotation.
    m,n -- the two necessary integers
    lim -- upper limit for the plane indices

    """
    uvw = np.array(uvw)
    Theta = get_cubic_theta(uvw, m, n)
    Sigma = get_cubic_sigma(uvw, m, n)
    R1 = rot(uvw, Theta)

    # List and character of possible GB planes:
    x = np.arange(-lim, lim + 1, 1)
    y = x
    z = x
    indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(len(x) ** 3, 3)
    indice_0 = indice[np.where(np.sum(abs(indice), axis=1) != 0)]
    indice_0 = indice_0[np.argsort(norm(indice_0, axis=1))]

    # extract the minimal cell:
    Min_1, Min_2 = create_minimal_cell_method_1(Sigma, uvw, R1)
    V1 = np.zeros([len(indice_0), 3])
    V2 = np.zeros([len(indice_0), 3])
    GBtype = []
    tol = 0.001
    # Mirrorplanes cubic symmetry
    MP = np.array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1],
                   [1, 1, 0],
                   [0, 1, 1],
                   [1, 0, 1],
                   ], dtype='float')
    # Find GB plane coordinates:
    for i in range(len(indice_0)):
        if common_divisor(indice_0[i])[1] <= 1:
            V1[i, :] = (indice_0[i, 0] * Min_1[:, 0] +
                        indice_0[i, 1] * Min_1[:, 1] +
                        indice_0[i, 2] * Min_1[:, 2])
            V2[i, :] = (indice_0[i, 0] * Min_2[:, 0] +
                        indice_0[i, 1] * Min_2[:, 1] +
                        indice_0[i, 2] * Min_2[:, 2])

    V1 = (V1[~np.all(V1 == 0, axis=1)]).astype(int)
    V2 = (V2[~np.all(V2 == 0, axis=1)]).astype(int)
    MeanPlanes = (V1 + V2) / 2

    # Check the type of GB plane: Symmetric tilt, tilt, twist
    for i in range(len(V1)):
        if ang(V1[i], uvw) < tol:

            for j in range(len(symmetry_equivalent(MP))):
                if 1 - ang(MeanPlanes[i], symmetry_equivalent(MP)[j]) < tol:
                    GBtype.append('Symmetric Tilt')
                    break
            else:
                GBtype.append('Tilt')
        elif 1 - ang(V1[i], uvw) < tol:
            GBtype.append('Twist')
        else:
            GBtype.append('Mixed')

    return (V1, V2, MeanPlanes, GBtype)


def find_orthogonal_cell(uvw, sigma, theta, R, m, n, GB1, 
    Min_1, Min_2, tol = 0.001):
    """
    finds Orthogonal cells from the CSL minimal cells.
    arguments:
    basis -- lattice basis
    uvw -- rotation axis
    m,n -- two necessary integers
    GB1 -- input plane orientation
    """
    # Inputs
    lim = 15
    uvw = np.array(uvw)
    Theta = theta
    Sigma = sigma
    
    #create GB2 from GB1
    GB2 = np.round((dot(R, GB1)), 6)
    x = np.arange(-lim, lim + 1, 1)
    y = x
    z = x
    
    indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(len(x) ** 3, 3)
    indice_0 = indice[np.where(np.sum(abs(indice), axis=1) != 0)]
    indice_0 = indice_0[np.argsort(norm(indice_0, axis=1))]
    OrthoCell_1 = np.zeros([3, 3])
    OrthoCell_1[:, 0] = np.array(GB1)
    OrthoCell_2 = np.zeros([3, 3])
    OrthoCell_2[:, 0] = np.array(GB2)
    
    # Find Ortho vectors:
    if ang(OrthoCell_1[:, 0], uvw) < tol:
        OrthoCell_1[:, 1] = uvw
        OrthoCell_2[:, 1] = uvw
    else:
        for i in range(len(indice_0)):
            v1 = (indice_0[i, 0] * Min_1[:, 0] +
                  indice_0[i, 1] * Min_1[:, 1] +
                  indice_0[i, 2] * Min_1[:, 2])
            v2 = (indice_0[i, 0] * Min_2[:, 0] +
                  indice_0[i, 1] * Min_2[:, 1] +
                  indice_0[i, 2] * Min_2[:, 2])
            if ang(v1, OrthoCell_1[:, 0]) < tol:
                OrthoCell_1[:, 1] = v1
                OrthoCell_2[:, 1] = v2
                break
    OrthoCell_1[:, 2] = np.cross(OrthoCell_1[:, 0], OrthoCell_1[:, 1])
    OrthoCell_2[:, 2] = np.cross(OrthoCell_2[:, 0], OrthoCell_2[:, 1])

    if (common_divisor(OrthoCell_1[:, 2])[1] ==
            common_divisor(OrthoCell_2[:, 2])[1]):
        OrthoCell_1[:, 2] = common_divisor(OrthoCell_1[:, 2])[0]
        OrthoCell_2[:, 2] = common_divisor(OrthoCell_2[:, 2])[0]
    
    Volume_1 = (round(det(OrthoCell_1), 5))
    Volume_2 = (round(det(OrthoCell_2), 5))
    
    if Volume_1 == Volume_2:
        return True, OrthoCell_1.astype(float), OrthoCell_2.astype(float)

    return False, None, None


def get_tilt_twist_comp(v1, uvw, m, n, tol=0.001):
    """
    returns the tilt and twist components of a given GB plane.
    arguments:
    v1 -- given gb plane
    uvw -- axis of rotation
    m,n -- the two necessary integers
    """
    theta = get_cubic_theta(uvw, m, n)
    R = rot(uvw, theta)
    v2 = np.round(dot(R, v1), 6).astype(int)
    tilt = angv(v1, v2)
    twist = 0
    gb_type = ""
    if abs(tilt - degrees(theta)) < 10e-5:
        gb_type = "tilt"
    else:
        twist = 2 * acos(cos(theta / 2) / cos(radians(tilt / 2)))
    
    #Assign types
    MP = np.array([[1, 0, 0],
                   [0, 1, 0],
                   [0, 0, 1],
                   [1, 1, 0],
                   [0, 1, 1],
                   [1, 0, 1],
                   ], dtype='float')
        

    return tilt, twist

#Helpers

def common_divisor(a):
    """
    returns the common factor of vector a and the reduced vector.
    """
    CommFac = []
    a = np.array(a)
    for i in range(2, 100):
        while (a[0] % i == 0 and a[1] % i == 0 and a[2] % i == 0):
            a = a / i
            CommFac.append(i)
    return(a.astype(int), (abs(np.prod(CommFac))))

def integer_array(A, tol=1e-7):
    """
    returns True if an array is ineteger.
    """
    return np.all(abs(np.round(A) - A) < tol)

def integer_matrix(a):
    """
    returns an integer matrix from row vectors.
    """
    Found = True
    b = np.zeros((3, 3))
    a = np.array(a)
    for i in range(3):
        for j in range(1, 2000):
            testV = j * a[i]
            if integer_array(testV):
                b[i] = testV
                break
        if all(b[i] == 0):
            Found = False
            print("Can not make integer matrix!")
    return (b) if Found else None

def angv(a, b):
    """
    returns the angle between two vectors.
    """
    ang = acos(np.round(dot(a, b)/norm(a)/norm(b), 8))
    return round(degrees(ang), 7)


def ang(a, b):
    """
    returns the cos(angle) between two vectors.
    """
    ang = np.round(dot(a, b)/norm(a)/norm(b), 7)
    return abs(ang)

def symmetry_equivalent(arr):
    """
    returns cubic symmetric eqivalents of the given 2 dimensional vector.
    """
    Sym = np.zeros([24, 3, 3])
    Sym[0, :] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    Sym[1, :] = [[1, 0, 0], [0, -1, 0], [0, 0, -1]]
    Sym[2, :] = [[-1, 0, 0], [0, 1, 0], [0, 0, -1]]
    Sym[3, :] = [[-1, 0, 0], [0, -1, 0], [0, 0, 1]]
    Sym[4, :] = [[0, -1, 0], [-1, 0, 0], [0, 0, -1]]
    Sym[5, :] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    Sym[6, :] = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]
    Sym[7, :] = [[0, 1, 0], [1, 0, 0], [0, 0, -1]]
    Sym[8, :] = [[-1, 0, 0], [0, 0, -1], [0, -1, 0]]
    Sym[9, :] = [[-1, 0, 0], [0, 0, 1], [0, 1, 0]]
    Sym[10, :] = [[1, 0, 0], [0, 0, -1], [0, 1, 0]]
    Sym[11, :] = [[1, 0, 0], [0, 0, 1], [0, -1, 0]]
    Sym[12, :] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    Sym[13, :] = [[0, 1, 0], [0, 0, -1], [-1, 0, 0]]
    Sym[14, :] = [[0, -1, 0], [0, 0, 1], [-1, 0, 0]]
    Sym[15, :] = [[0, -1, 0], [0, 0, -1], [1, 0, 0]]
    Sym[16, :] = [[0, 0, 1], [1, 0, 0], [0, 1, 0]]
    Sym[17, :] = [[0, 0, 1], [-1, 0, 0], [0, -1, 0]]
    Sym[18, :] = [[0, 0, -1], [1, 0, 0], [0, -1, 0]]
    Sym[19, :] = [[0, 0, -1], [-1, 0, 0], [0, 1, 0]]
    Sym[20, :] = [[0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    Sym[21, :] = [[0, 0, -1], [0, 1, 0], [1, 0, 0]]
    Sym[22, :] = [[0, 0, 1], [0, -1, 0], [1, 0, 0]]
    Sym[23, :] = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

    arr = np.atleast_2d(arr)
    Result = []
    for i in range(len(Sym)):
        for j in range(len(arr)):
            Result.append(dot(Sym[i, :], arr[j]))
    Result = np.array(Result)
    return np.unique(Result, axis=0)


# GB building tools
def generate_ortho_unitcell_atoms(ortho, basis):

    """
    populates a unitcell from the orthogonal vectors.
    """
    Or = ortho.T
    Orint = integer_matrix(Or)
    LoopBound = np.zeros((3, 2), dtype=float)
    transformed = []
    CubeCoords = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0],
                          [0, 1, 1], [1, 0, 1], [1, 1, 1], [0, 0, 0]],
                          dtype=float)
    for i in range(len(CubeCoords)):
        transformed.append(np.dot(Orint.T, CubeCoords[i]))

    # Finding bounds for atoms in a CSL unitcell:
    LoopBound[0, :] = [min(np.array(transformed)[:, 0]),
                       max(np.array(transformed)[:, 0])]
    LoopBound[1, :] = [min(np.array(transformed)[:, 1]),
                       max(np.array(transformed)[:, 1])]
    LoopBound[2, :] = [min(np.array(transformed)[:, 2]),
                       max(np.array(transformed)[:, 2])]

    # Filling up the unitcell:
    Tol = 1
    x = np.arange(LoopBound[0, 0] - Tol, LoopBound[0, 1] + Tol + 1, 1)
    y = np.arange(LoopBound[1, 0] - Tol, LoopBound[1, 1] + Tol + 1, 1)
    z = np.arange(LoopBound[2, 0] - Tol, LoopBound[2, 1] + Tol + 1, 1)
    V = len(x) * len(y) * len(z)
    indice = (np.stack(np.meshgrid(x, y, z)).T).reshape(V, 3)
    Base = basis
    Atoms = []
    tol = 0.001
    if V > 5e6:
        print("Warning! It may take a very long time"
              "to produce this cell!")
    # produce Atoms:

    for i in range(V):
        for j in range(len(Base)):
            Atoms.append(indice[i, 0:3] + Base[j, 0:3])
    Atoms = np.array(Atoms)

    # Cell conditions
    Con1 = dot(Atoms, Or[0]) / norm(Or[0]) + tol
    Con2 = dot(Atoms, Or[1]) / norm(Or[1]) + tol
    Con3 = dot(Atoms, Or[2]) / norm(Or[2]) + tol
    # Application of the conditions:
    Atoms = (Atoms[(Con1 >= 0) & (Con1 <= norm(Or[0])) & (Con2 >= 0) &
             (Con2 <= norm(Or[1])) &
             (Con3 >= 0) & (Con3 <= norm(Or[2]))])

    if len(Atoms) == (round(det(Or) * len(Base), 7)).astype(int):
        return Atoms

def generate_bicrystal_atoms(ortho1, ortho2, basis):
    """
    builds the unitcells for both grains g1 and g2.
    """
    Or_1 = ortho1.T
    Or_2 = ortho2.T
    rot1 = np.array([Or_1[0, :] / norm(Or_1[0, :]),
                         Or_1[1, :] / norm(Or_1[1, :]),
                         Or_1[2, :] / norm(Or_1[2, :])])
    rot2 = np.array([Or_2[0, :] / norm(Or_2[0, :]),
                         Or_2[1, :] / norm(Or_2[1, :]),
                         Or_2[2, :] / norm(Or_2[2, :])])

    atoms1 = generate_ortho_unitcell_atoms(ortho1.copy(), basis)
    atoms2 = generate_ortho_unitcell_atoms(ortho2.copy(), basis)

    atoms1 = dot(rot1, atoms1.T).T
    atoms2 = dot(rot2, atoms2.T).T
    
    atoms2[:, 0] = atoms2[:, 0] - norm(Or_2[0, :])  # - tol
    return atoms1, atoms2


def expand_super_cell(ortho1, atoms1, atoms2, dim):
    """
    expands the smallest CSL unitcell to the given dimensions.
    """
    a = norm(ortho1[:, 0])
    b = norm(ortho1[:, 1])
    c = norm(ortho1[:, 2])
    dimX, dimY, dimZ = dim

    X = atoms1.copy()
    Y = atoms2.copy()

    X_new = []
    Y_new = []
    for i in range(dimX):
        for j in range(dimY):
            for k in range(dimZ):
                Position1 = [i * a, j * b, k * c]
                Position2 = [-i * a, j * b, k * c]
                for l in range(len(X)):
                    X_new.append(Position1[0:3] + X[l, 0:3])
                for m in range(len(Y)):
                    Y_new.append(Position2[0:3] + Y[m, 0:3])

    atoms1 = np.array(X_new)
    atoms2 = np.array(Y_new)
    return atoms1, atoms2

def find_overlapping_atoms(ortho1, atoms1, atoms2, dim, overlap=0.0):
    """
    finds the overlapping atoms.
    """
    periodic_length = norm(ortho1[:, 0]) * dim[0]
    periodic_image = atoms2 + [periodic_length * 2, 0, 0]
    # select atoms contained in a smaller window around the GB and its
    # periodic image
    IndX = np.where([(atoms1[:, 0] < 1) |
                     (atoms1[:, 0] > (periodic_length - 1))])[1]
    IndY = np.where([atoms2[:, 0] > -1])[1]
    IndY_image = np.where([periodic_image[:, 0] <
                          (periodic_length + 1)])[1]
    X_new = atoms1[IndX]
    Y_new = np.concatenate((atoms2[IndY], periodic_image[IndY_image]))
    IndY_new = np.concatenate((IndY, IndY_image))
    # create a meshgrid search routine
    x = np.arange(0, len(X_new), 1)
    y = np.arange(0, len(Y_new), 1)
    indice = (np.stack(np.meshgrid(x, y)).T).reshape(len(x) * len(y), 2)
    norms = norm(X_new[indice[:, 0]] - Y_new[indice[:, 1]], axis=1)
    indice_x = indice[norms < overlap][:, 0]
    indice_y = indice[norms < overlap][:, 1]
    X_del = X_new[indice_x]
    Y_del = Y_new[indice_y]

    if (len(X_del) != len(Y_del)):
        print("Warning! the number of deleted atoms"
              "in the two grains are not equal!")
    # print(type(IndX), len(IndY), len(IndY_image))
    return (X_del, Y_del, IndX[indice_x], IndY_new[indice_y])

def populate_gb(ortho1, ortho2, basis, 
    lattice_parameter, dim=(1,1,1), 
    overlap=0.0, rigid=False):
    
    atoms1, atoms2 = generate_bicrystal_atoms(ortho1, ortho2, basis)
    if overlap > 0.0:
        which_gb = "g1"
    elif overlap < 0.0:
        which_gb = "g2"
    else:
        which_gb = "a"

    atoms1, atoms2 = expand_super_cell(ortho1, atoms1, atoms2, dim)
    xdel, _, x_indice, y_indice = find_overlapping_atoms(ortho1, atoms1, atoms2, dim, overlap=overlap)

    if which_gb == "g1":
        atoms1 = np.delete(atoms1, x_indice, axis=0)
    elif which_gb == "g2":
        atoms2 = np.delete(atoms2, y_indice, axis=0)

    atoms1 = atoms1*lattice_parameter
    atoms2 = atoms2*lattice_parameter
    dimx, dimy, dimz = dim

    #get box bounds
    xlo = -1 * np.round(norm(ortho1[:, 0]) * dimx * lattice_parameter, 8)
    xhi = np.round(norm(ortho1[:, 0]) * dimx * lattice_parameter, 8)
    LenX = xhi - xlo
    ylo = 0.0
    yhi = np.round(norm(ortho1[:, 1]) * dimy * lattice_parameter, 8)
    LenY = yhi - ylo
    zlo = 0.0
    zhi = np.round(norm(ortho1[:, 2]) * dimz * lattice_parameter, 8)
    LenZ = zhi - zlo    

    return [[LenX, 0, 0],[0, LenY, 0],[0, 0, LenZ]], atoms1, atoms2


