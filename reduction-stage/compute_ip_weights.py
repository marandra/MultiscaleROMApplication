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


def ComputeROQ(Modes, weights, nGP, factorLEQ, tol):
    [J, b] = ComputeJandb(Modes, weights, factorLEQ)[:2]
    #M = len(weights)
    y = np.arange(len(weights))

    # resudual vector, initial guess
    r = b

    # number of iterations
    it = 0
    mPOS = 0
    z = []
    Jnorm = np.sqrt(sum(np.multiply(J, J), 0))

    # point Selection Algorithm
    while (np.linalg.norm(r) / np.linalg.norm(b) > tol) and (mPOS <= nGP):
        # 1. Compute new point
        ObjFun = np.dot((J[:, y]).T, r)
        div = np.multiply(Jnorm[y], np.linalg.norm(r))
        ObjFun = np.divide(ObjFun, div)
        s = ObjFun.argmax()
        t = y[s]
        # 2. Move i from set y to set z
        z = (np.append(z, t)).astype(int)
        y = np.delete(y, s)
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
        #TODO: is iterator really needed? Iteration counter
        it = it + 1
        logger.debug("k = {}, mPOS = {}, error = {:.2f}%".format(it, mPOS, np.linalg.norm(r)/np.linalg.norm(b) * 100))
    # 6. Postprocess of points - neglecting null weights
    INDzero = np.where(x == 0)[0]
    if any(INDzero):
        z = np.delete(z, INDzero)
    w = np.multiply(x, np.sqrt(weights[z]))
    logger.debug("Reduced Weights: {}".format(w))
    logger.debug("sum of reduced weights: {}".format(np.sum(w)))
    logger.debug("GP's index: {}".format(z))
    return w, z


def write_bases(filename, U):
    #file_format = conf['Parameters']['roq_file_format']
    #if file_format == 'ascii':
    #else:
    np.savetxt(filename, U)
    return


def compute_reduced_set(conf):
    nr_roq_points = int(conf['Parameters']['nr_roq_points'])
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    energy_bases_filename = conf['Parameters']['energy_bases_filename']
    integration_weights_filename = conf['Parameters']['integration_weights_filename']
    integration_weights = np.loadtxt(integration_weights_filename)
    energy_modes = np.load(energy_bases_filename)
    [w, z] = ComputeROQ(energy_modes, integration_weights,
                        nr_roq_points, factorLEQ=1.0, tol=1.e-10)
    roq_weigths = -1 * np.ones([nr_elements, nr_integration_points])
    for x, igg in enumerate(z):
        e = int(igg / nr_integration_points)
        ig = igg % nr_integration_points
        roq_weigths[e][ig] = w[x]
    return roq_weigths


def create_rom_weights(conf):
    integration_weights_filename = conf['Parameters']['integration_weights_filename']
    integration_weights = np.loadtxt(integration_weights_filename)
    nr_elements = int(conf['Parameters']['nr_elements'])
    return integration_weights.reshape((nr_elements, -1))


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Computes Reduced Order Quadrature (ROQ) integration weights")
parser.add_argument('config_file', help="configuration file")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
parser.add_argument('-r', '--rom', action="store_true", help="compute ROM instead of HPROM")
args = parser.parse_args()

# parse configuration file
conf = configparser.ConfigParser()
conf.read(args.config_file)

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s',
                    datefmt='%H:%M:%S',level=verbosity_level)
logger = logging.getLogger(__name__)
handler = logging.FileHandler('log_' + args.config_file.rsplit('.', 1)[0])
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)

if __name__ == '__main__':
    logger.info("Reduced Order Quadrature")
    if args.rom:
        logger.info("Computing ROM")
        roq = create_rom_weights(conf)
    else:
        logger.info("Computing HPROM")
        roq = compute_reduced_set(conf)
    filename = conf['Parameters']['roq_weights_filename']
    np.savetxt(filename, roq)
