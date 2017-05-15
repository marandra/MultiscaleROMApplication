import sys
import configparser
import glob
import numpy as np
import lsqnonneg
import numpy

eps = 2.22e-16    # from matlab 

def lsqnonneg(C, d, x0=None, tol=None, itmax_factor=10):
    '''Linear least squares with nonnegativity constraints.
    (x, resnorm, residual) = lsqnonneg(C,d) returns the vector x that minimizes norm(d-C*x)
    subject to x >= 0, C and d must be real
    '''
 
    eps = 2.22e-16    # from matlab

    def norm1(x):
        return abs(x).sum().max()

 
    def msize(x, dim):
        s = x.shape
        if dim >= len(s):
            return 1
        else:
            return s[dim]


    if tol is None:
        tol = 10*eps*norm1(C)*(max(C.shape)+1)
    C = numpy.asarray(C)
    (m,n) = C.shape
    P = numpy.zeros(n)
    Z = numpy.arange(1, n+1)
    if x0 is None:
        x = P
    else:
        if any(x0 < 0):
            x = P
        else:
            x = x0
    ZZ = Z
    resid = d - numpy.dot(C, x)
    w = numpy.dot(C.T, resid)
    outeriter = 0
    it = 0
    itmax = itmax_factor * n
    exitflag = 1

    # outer loop to put variables into set to hold positive coefficients
    while numpy.any(Z) and numpy.any(w[ZZ - 1] > tol):
        outeriter += 1
        t = w[ZZ - 1].argmax()
        t = ZZ[t]
        P[t - 1] = t
        Z[t - 1] = 0
        PP = numpy.where(P <> 0)[0] + 1
        ZZ = numpy.where(Z <> 0)[0] + 1
        CP = numpy.zeros(C.shape)
        CP[:, PP - 1] = C[:, PP - 1]
        CP[:, ZZ - 1] = numpy.zeros((m, msize(ZZ, 1)))
        z = numpy.dot(numpy.linalg.pinv(CP), d)
        z[ZZ - 1] = numpy.zeros((msize(ZZ, 1), msize(ZZ, 0)))

        # inner loop to remove elements from the positve set which no longer belong
        while numpy.any(z[PP-1] <= tol):
            it += 1
            if it > itmax:
                max_error = z[PP - 1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n Try raising the \
                                 tolerance tol. (max_error=%d)' % (it, max_error))
            QQ = numpy.where((z <= tol) & (P <> 0))[0]
            alpha = min(x[QQ] / (x[QQ] - z[QQ]))
            x = x + alpha*(z - x)
            ij = numpy.where((abs(x) < tol) & (P <> 0))[0] + 1
            Z[ij - 1] = ij
            P[ij - 1] = numpy.zeros(max(ij.shape))
            PP = numpy.where(P <> 0)[0] + 1
            ZZ = numpy.where(Z <> 0)[0] + 1
            CP[:, PP - 1] = C[:, PP-1]
            CP[:, ZZ - 1] = numpy.zeros((m, msize(ZZ, 1)))
            z=numpy.dot(numpy.linalg.pinv(CP), d)
            z[ZZ - 1] = numpy.zeros((msize(ZZ, 1), msize(ZZ, 0)))
        x = z
        resid = d - numpy.dot(C, x)
        w = numpy.dot(C.T, resid)
    return (x, sum(resid * resid), resid)


def ComputeJandb(Modes, weights, factorLEQ=None):
    if factorLEQ is None:
        factorLEQ = 1.0
    # Exact integral - numerical integration
    INTexact = numpy.dot(Modes.T,weights)

    # Total microscale volume
    vol = numpy.sum(weights)
    sqrtVol = numpy.sqrt(weights)

    # Matrix of modified modes (with zero integral)
    Xf = numpy.zeros(Modes.shape)

    # Loops over the initial modes
    Xf = numpy.subtract(Modes, INTexact / vol)
    Xf = numpy.multiply(Xf.T, sqrtVol).T

    # Singular Value Decomposition
    [Lambda,SValues,VValues] = numpy.linalg.svd(Xf,full_matrices=False)

    # fixed tolerance to define the reduced modified set of modes
    tol = numpy.max(Modes.shape) * eps * numpy.max(SValues)
    RankXf = sum(i > tol for i in SValues)
    Lambda = Lambda[:,0:RankXf]
    J = Lambda.T
    Jw = factorLEQ * sqrtVol / numpy.sqrt(vol)

    # Adding last row related with the sqrt of gauss integration weigths
    J = numpy.vstack([J, Jw.T])

    # Initializing the RHS vector for the optimization problem
    b = numpy.append(numpy.zeros(Lambda.T.shape[0]), factorLEQ * numpy.sqrt(vol))
    return (J, b, INTexact)

# get parameters
fname = sys.argv[1]
conf = configparser.ConfigParser()
conf.read(fname)
file_format = conf['Parameters']['file_format']
nr_elastic_snapshots = int(conf['Parameters']['nr_elastic_snapshots'])
nr_elements = int(conf['Parameters']['nr_elements'])
nr_integration_points = int(conf['Parameters']['nr_integration_points'])
nr_strain_components = int(conf['Parameters']['nr_strain_components'])
energy_file_name = conf['Parameters']['energy_file_name']
strain_file_name = conf['Parameters']['strain_file_name']

trajectory_paths = glob.glob('*_?')
energy_elastic_files = []
strain_elastic_files = []
energy_inelast_files = []
strain_inelast_files = []
for path in trajectory_paths:
    energy_elastic_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[:nr_elastic_snapshots])
    strain_elastic_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[:nr_elastic_snapshots])
    energy_inelast_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[nr_elastic_snapshots:])
    strain_inelast_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[nr_elastic_snapshots:])
    
#import pprint
#print("AAAAAAA")
#pprint.pprint(energy_elastic_files)
#print("BBBBBB")
#pprint.pprint(strain_elastic_files)
#print("CCCCCCC")
#pprint.pprint(energy_inelast_files)
#print("DDDDDDD")
#pprint.pprint(strain_inelast_files)


# first part: read and compute elastic energy modes, compute projector
nr_dofs = nr_elements * nr_integration_points
X = np.empty([nr_dofs, len(energy_elastic_files)])
for i, file in enumerate(energy_elastic_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
Ue_el = np.linalg.svd(X, full_matrices=False)[0]
print(Ue_el)
print("Fin svd strain elastic")
print(X.shape)

sys.exit()

# second part: read inelastic energy modes, remove elastic component, decomp svd
X = np.empty([nr_dofs, len(energy_inelast_files)])
for i, file in enumerate(energy_inelast_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
print("Antes de proyector")
X = X - np.dot(Ue_el, np.dot(Ue_el.T, X))
print("Desp de proyector")
Ue_in = np.linalg.svd(X, full_matrices=False)[0]
print("Fin svd energy inelastic")

# third part: read and compute elastic strain modes
nr_dofs = nr_elements * nr_integration_points * nr_strain_components
X = np.empty([nr_dofs, len(strain_elastic_files)])
for i, file in enumerate(strain_elastic_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
Us_el = np.linalg.svd(X, full_matrices=False)[0]
print("Fin svd strain elastic")

# forth part: read inelastic strain modes, remove elastic component, decomp svd
X = np.empty([nr_dofs, len(strain_inelast_files)])
for i, file in enumerate(strain_inelast_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
print("Antes de proyector")
X = X - np.dot(Us_el, np.dot(Us_el.T, X))
print("Desp de proyector")
Us_in = np.linalg.svd(X, full_matrices=False)[0]
print("fin")
