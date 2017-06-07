import sys
import configparser
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


#######################################
# Main
#######################################
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
if __name__ == '__main__':
    handler = logging.FileHandler(sys.argv[1].rsplit('.', 1)[0] + '.log')
    handler.setLevel(logging.INFO)
    logger.addHandler(handler)

    # get parameters
    conf = configparser.ConfigParser()
    conf.read(sys.argv[1])
    snapshot_file_format = conf['Parameters']['snapshot_file_format']
    #bases_file_format = conf['Parameters']['bases_file_format']
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_strain_components = int(conf['Parameters']['nr_strain_components'])
    nr_energy_reduced_modes = int(conf['Parameters']['nr_energy_reduced_modes'])
    trajectory_filename = conf['Parameters']['trajectory_filename']
    energy_filename = conf['Parameters']['energy_filename']
    strain_filename = conf['Parameters']['strain_filename']
    integration_weights_filename = conf['Parameters']['integration_weights_filename']
    nr_elastic_snapshots_filename = conf['Parameters']['nr_elastic_snapshots_filename']
    energy_bases_filename = conf['Parameters']['energy_bases_filename']
    strain_bases_filename = conf['Parameters']['strain_bases_filename']
    roq_weights_filename = conf['Parameters']['roq_weights_filename']
    tolerance_svd_elastic_strain = float(conf['Parameters']['tolerance_svd_elastic_strain'])
    tolerance_svd_elastic_energy = float(conf['Parameters']['tolerance_svd_elastic_energy'])
    tolerance_svd_inelastic_strain = float(conf['Parameters']['tolerance_svd_inelastic_strain'])
    tolerance_svd_inelastic_energy = float(conf['Parameters']['tolerance_svd_inelastic_energy'])

    # get files
    integration_weights_files = glob.glob("{}_?/{}".format(trajectory_filename, integration_weights_filename))
    elastic_snapshots_files = glob.glob("{}_?/{}".format(trajectory_filename, nr_elastic_snapshots_filename))
    trajectory_paths = sorted(glob.glob("{}_?".format(trajectory_filename)))

    # TODO : for the future, take into account loadin/unloading trajectory cases
    elastic_mode_traj = []
    for filename in elastic_snapshots_files:
        with open(filename, "r") as f:
            elastic_mode_traj.append(int(f.readline().strip()))
    nr_elastic_snapshots = min(elastic_mode_traj)
    logger.info("Nr of elastic snapshots: {}".format(nr_elastic_snapshots))

    energy_elastic_files = []
    energy_inelast_files = []
    strain_elastic_files = []
    strain_inelast_files = []
    for path in trajectory_paths:
        energy_elastic_files.extend(sorted(glob.glob("{}/{}*".format(path, energy_filename)))[:nr_elastic_snapshots])
        energy_inelast_files.extend(sorted(glob.glob("{}/{}*".format(path, energy_filename)))[nr_elastic_snapshots:])
        strain_elastic_files.extend(sorted(glob.glob("{}/{}*".format(path, strain_filename)))[:nr_elastic_snapshots])
        strain_inelast_files.extend(sorted(glob.glob("{}/{}*".format(path, strain_filename)))[nr_elastic_snapshots:])

    logger.info("STRAIN SNAPSHOTS")

    logger.info("Step 01: SVD of strain elastic snapshots")
    nr_dofs = nr_elements * nr_integration_points * nr_strain_components
    X = np.empty([nr_dofs, len(strain_elastic_files)])
    for i, file in enumerate(strain_elastic_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
    # TODO: incluir el SVD del paquete scipy, parece que es mas optimo, ver si esta instalada una version reciente de scipy en el cluster
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    logger.info(S)

    logger.info("Step 02: selection process of strain elastic modes")
    cont = 1
    Us_el = []
    for iValue, singular_value in enumerate(S):
        if singular_value > tolerance_svd_elastic_strain:
            if cont == 1:
                Us_el = U[:, iValue]
            else:
                Us_el = np.column_stack((Us_el, U[:, iValue]))
            cont = cont + 1
    logger.info(Us_el.shape)

    logger.info("Step 03: projection of strain inelastic snapshots")
    X = np.empty([nr_dofs, len(strain_inelast_files)])
    for i, file in enumerate(strain_inelast_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        for j in range(Us_el.shape[1]):
            X[:,i] = X[:,i] - np.multiply(np.dot(Us_el[:,j],X[:,i]),Us_el[:,j])

    logger.info("Step 04: SVD of strain inelastic modified snapshots")
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    logger.info(S)


    logger.info("Step 05: Selection process of strain inelastic modes")
    cont = 1
    Us_in = []
    for iValue, singular_value in enumerate(S):
        if singular_value > tolerance_svd_inelastic_strain:
            if cont == 1:
                Us_in = U[:, iValue]
            else:
                Us_in = np.column_stack((Us_in, U[:, iValue]))
            cont = cont + 1
    logger.info(Us_in.shape)


    logger.info("Step 06: assembly of global matrix of strain modes")
    Us = np.hstack([Us_el, Us_in])


    logger.info("ENERGY SNAPSHOTS")


    logger.info("Step 07: SVD of elastic energy snapshots")
    nr_dofs = nr_elements * nr_integration_points
    X = np.empty([nr_dofs, len(energy_elastic_files)])
    for i, file in enumerate(energy_elastic_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    logger.info(S)

    logger.info("Step 08: Selection process of elastic energy modes")
    cont = 1
    Ue_el = []
    for iValue, sin_val in enumerate(S):
        if sin_val > tolerance_svd_elastic_energy:
            if cont == 1:
                Ue_el = U[:, iValue]
            else:
                Ue_el = np.column_stack((Ue_el, U[:, iValue]))
            cont = cont + 1
    logger.info(Ue_el.shape)

    logger.info("Step 09: projection of energy inelastic snapshots")
    X = np.empty([nr_dofs, len(energy_inelast_files)])
    for i, file in enumerate(energy_inelast_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        for j in range(Ue_el.shape[1]):
            X[:,i] = X[:,i] - np.multiply(np.dot(Ue_el[:,j],X[:,i]),Ue_el[:,j])

    logger.info("Step 10: svd of inelastic energy modified snapshots")
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    logger.info(S)


    logger.info("Step 11: selection process of inelastic energy modes")
    cont = 1
    Ue_in = []
    for iValue, sin_val in enumerate(S):
        if sin_val > tolerance_svd_inelastic_energy:
            if cont == 1:
                Ue_in = U[:, iValue]
            else:
                Ue_in = np.column_stack((Ue_in, U[:, iValue]))
            cont = cont + 1
    logger.info(Ue_in.shape)

    logger.info("Step 12: assembly of global matrix of energy modes")
    Ue=np.hstack([Ue_el, Ue_in])

    # reading the gauss weights for computing the ROQ
    nr_dofs = nr_elements * nr_integration_points
    gauss_weights = np.loadtxt(integration_weights_files[0])

    # TODO: change the print format of the matrix in order to avoid wrong tabulation because of the minus (-) sign.
    logger.info("Printing data to files")
    with open(strain_bases_filename, 'wb') as ofile:
        np.savetxt(ofile, Us, fmt='%.13f')

    if nr_energy_reduced_modes > Ue.shape[1]:
        sys.exit("Error: number of energy modes greater than the total number of computed energy modes")

    Ue_red = Ue[:, 0:nr_energy_reduced_modes]

    with open(energy_bases_filename, 'wb') as ofile:
       np.savetxt(ofile, Ue_red, fmt='%.13f')

    logger.info("COMPUTING REDUCED ORDER QUADRATURE (ROQ)")
    factorLEQ = 1.0
    tol = 1e-10
    nGP = nr_energy_reduced_modes
    [w, z] = ComputeROQ(Ue_red, gauss_weights, factorLEQ, nGP, tol)
    roq_weigths = np.empty([nr_elements, nr_integration_points])
    for i in range(nr_elements):
        for j in range(nr_integration_points):
            i_elem = nr_integration_points * i + j
            d = np.where(z==i_elem)[0]
            if not d:
                roq_weigths[i][j] = -1
            else:
                roq_weigths[i][j] = w[d[0]]
    with open(roq_weights_filename, 'wb') as ofile:
        np.savetxt(ofile, roq_weigths, fmt='%.13f')