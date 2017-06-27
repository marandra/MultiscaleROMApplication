import time
import configparser
import argparse
import glob
import numpy as np
import logging

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
        tol = 10 * eps * norm1(C) * (max(C.shape) + 1)
    C = np.asarray(C)
    (m,n) = C.shape
    P = np.zeros(n)
    Z = np.arange(1, n + 1)
    if x0 is None:
        x = P
    else:
        if any(x0 < 0):
            x = P
        else:
            x = x0
    ZZ = Z
    resid = d - np.dot(C, x)
    w = np.dot(C.T, resid)
    outeriter = 0
    it = 0
    itmax = itmax_factor * n
    exitflag = 1

    # outer loop to put variables into set to hold positive coefficients
    while np.any(Z) and np.any(w[ZZ - 1] > tol):
        outeriter += 1
        t = w[ZZ - 1].argmax()
        t = ZZ[t]
        P[t - 1] = t
        Z[t - 1] = 0
        PP = np.where(P != 0)[0] + 1
        ZZ = np.where(Z != 0)[0] + 1
        CP = np.zeros(C.shape)
        CP[:, PP - 1] = C[:, PP - 1]
        CP[:, ZZ - 1] = np.zeros((m, msize(ZZ, 1)))
        z = np.dot(np.linalg.pinv(CP), d)
        z[ZZ - 1] = np.zeros((msize(ZZ, 1), msize(ZZ, 0)))

        # inner loop to remove elements from the positve set which no longer belong
        while np.any(z[PP-1] <= tol):
            it += 1
            if it > itmax:
                max_error = z[PP - 1].max()
                raise Exception('Exiting: Iteration count (=%d) exceeded\n Try raising the \
                                 tolerance tol. (max_error=%d)' % (it, max_error))
            QQ = np.where((z <= tol) & (P != 0))[0]
            alpha = min(x[QQ] / (x[QQ] - z[QQ]))
            x = x + alpha * (z - x)
            ij = np.where((abs(x) < tol) & (P != 0))[0] + 1
            Z[ij - 1] = ij
            P[ij - 1] = np.zeros(max(ij.shape))
            PP = np.where(P != 0)[0] + 1
            ZZ = np.where(Z != 0)[0] + 1
            CP[:, PP - 1] = C[:, PP - 1]
            CP[:, ZZ - 1] = np.zeros((m, msize(ZZ, 1)))
            z=np.dot(np.linalg.pinv(CP), d)
            z[ZZ - 1] = np.zeros((msize(ZZ, 1), msize(ZZ, 0)))
        x = z
        resid = d - np.dot(C, x)
        w = np.dot(C.T, resid)
    return (x, sum(resid * resid), resid)


def ComputeJandb(Modes, weights, factorLEQ=1.0):

    eps = 2.22e-16    # from matlab

    # Exact integral - numerical integration
    INTexact = np.dot(Modes.T, weights)

    # Total microscale volume
    vol = np.sum(weights)
    sqrtVol = np.sqrt(weights)

    # Matrix of modified modes (with zero integral)
    Xf = np.zeros(Modes.shape)

    # Loops over the initial modes
    Xf = np.subtract(Modes, INTexact / vol)
    Xf = np.multiply(Xf.T, sqrtVol).T

    # Singular Value Decomposition
    [Lambda,SValues,VValues] = np.linalg.svd(Xf,full_matrices=False)

    # fixed tolerance to define the reduced modified set of modes
    tol = np.max(Modes.shape) * eps * np.max(SValues)
    RankXf = sum(i > tol for i in SValues)
    Lambda = Lambda[:,0:RankXf]
    J = Lambda.T
    Jw = factorLEQ * sqrtVol / np.sqrt(vol)

    # Adding last row related with the sqrt of gauss integration weigths
    J = np.vstack([J, Jw.T])

    # Initializing the RHS vector for the optimization problem
    b = np.append(np.zeros(Lambda.T.shape[0]), factorLEQ * np.sqrt(vol))
    return J, b, INTexact


def ComputeROQ(Modes, weights, factorLEQ, nGP, tol):
    # computation of integration points weights
    [J, b, INTexact] = ComputeJandb(Modes, weights, factorLEQ)
    M = len(weights)
    y = np.arange(M)

    # Resudual vector, initial guess
    r = b

    #sys.exit()

    # Number of iterations
    it = 0
    mPOS=0
    z = []
    Jnorm = np.sqrt(sum(np.multiply(J,J),0))

    #print(J.shape)

    #sys.exit()

    # Point Selection Algorithm
    while (np.linalg.norm(r) / np.linalg.norm(b) > tol and mPOS <= nGP):

        # 1. Compute new point
        ObjFun = np.dot((J[:, y]).T, r)
        div = np.multiply(Jnorm[y], np.linalg.norm(r))
        ObjFun = np.divide(ObjFun, div)
        s = ObjFun.argmax()
        t = y[s]

        # 2. Move i from set y to set z
        z = (np.append(z, t)).astype(int)
        y=np.delete(y, s)

        # 3. solving LS conventional problem
        x = np.linalg.lstsq(J[:, z], b)[0]
        if any(x < 0):
            # 3. Determime alpha for solving a NNLS
            [x, resnorm, residual] = lsqnonneg(J[:,z], b)

        # 3. Determime alpha for solving a NNLS
        #[x, resnorm, residual] = lsqnonneg(J[:,z], b)

        # 4. Update the residual
        r = b - np.dot(J[:,z],x)

        # 5. Update mPOS and k
        mPOS = len(np.where(x>0)[0])

        # Iteration counter
        it = it + 1
        logger.info("k = {}, mPOS = {}, error = {:.2f}%".format(it, mPOS, np.linalg.norm(r)/np.linalg.norm(b) * 100))

    # 6. Postprocess of points - neglecting null weights
    INDzero = np.where(x == 0)[0]
    if any(INDzero):
        z = np.delete(z, INDzero)
    w = np.multiply(x, np.sqrt(weights[z]))

    logger.info("Reduced Weights")
    logger.info(w)
    #
    logger.info("sum of reduced weights")
    logger.info(np.sum(w))
    #
    logger.info("GP's index")
    logger.info(z)

    return(w, z)


def compute_reduced_set():
    # get parameters
    conf = configparser.ConfigParser()
    conf.read(config_filename)
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_energy_reduced_modes = int(conf['Parameters']['nr_energy_reduced_modes'])
    integration_weights_filename = conf['Parameters']['integration_weights_filename']
    energy_bases_filename = conf['Parameters']['energy_bases_filename']
    trajectory_filename = conf['Parameters']['trajectory_filename']
    roq_weights_filename = conf['Parameters']['roq_weights_filename']

    logger.info("REDUCED ORDER QUADRATURE")

    trajectory_paths = sorted(glob.glob("{}_?".format(trajectory_filename)))
    integration_weights = np.loadtxt(integration_weights_filename)
    energy_modes = np.loadtxt(energy_bases_filename)

    if nr_energy_reduced_modes > energy_modes.shape[1]:
        sys.exit("Error: number of energy modes greater than the total number of computed energy modes")
    energy_modes_reduced = energy_modes[:, 0:nr_energy_reduced_modes]

    #  subset of modes, or define a value (number of modes) as an input
    # TODO: Define a criteron to choose a
    factorLEQ = 1.0
    tol = 1e-10
    nGP = energy_modes.shape[1] #In case of use the same number of points as energy modes.
    [w, z] = ComputeROQ(energy_modes_reduced, integration_weights, factorLEQ, nGP, tol)

    # print matrix with new weigths
    roq_weigths = -1 * np.ones([nr_elements, nr_integration_points])
    for x, igg in enumerate(z):
        e = int(igg / 4)
        ig = igg % 4
        roq_weigths[e][ig] = w[x]

    with open(roq_weights_filename,'wb') as ofile:
        np.savetxt(ofile, roq_weigths, fmt='%.17f')

def create_rom_weights():
    fo = open("gauss_weights_rom", 'w')
    with open(sys.argv[1], 'r') as fi:
        for j in range(27000):
            for i in range(8):
                line = fi.readline().strip()
                fo.write('{}  '.format(line))
            fo.write('\n')

#######################################
# Main
#######################################
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
parser = argparse.ArgumentParser(description="Computes Reduced Order Quadrature (ROQ) integration weights)
parser.add_argument('config_file', help="configuration file")
parser.add_argument('-r', '--rom', action="store_true", help="compute ROM instead of HPROM")
args = parser.parse_args()
#handler = logging.FileHandler(config_filename.rsplit('.', 1)[0] + '.log')
#handler.setLevel(logging.INFO)
#logger.addHandler(handler)

if __name__ == '__main__':
    compute_reduced_set()
