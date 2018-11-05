import configparser
import argparse
import numpy as np
import logging
import json


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


def UpdateWeightsInverse(H, alpha, bases_current, base_new, r):
    c = np.dot(bases_current.T, base_new)
    d = np.dot(H, c).reshape(-1, 1)
    s = np.dot(base_new.T, base_new) - np.dot(c.T, d)
    aux1 = np.hstack([H + np.outer(d, d)/s, -d/s])
    aux2 = np.hstack([np.squeeze(-d.T/s), 1/s])
    H_new = np.vstack([aux1, aux2])
    v = np.dot(base_new.T, r) / s
    alpha = np.vstack([(alpha - d * v), v])
    return H_new, alpha


def UpdateInverseHermitian(Binv, jrow):
    if jrow == np.shape(Binv)[1]:
        aux = (Binv[0:-1, -1] * Binv[-1, 0:-1]) / Binv(-1, -1)
        Ahinv = Binv[:-1, :-1] - aux
    else:
        aux1 = np.hstack([Binv[:, 0:jrow], Binv[:, jrow+1:], Binv[:, jrow].reshape(-1,1)])
        aux2 = np.vstack([aux1[0:jrow, :], aux1[jrow+1:, :], aux1[jrow, :]])
        Ahinv = aux2[0:-1, 0:-1] - np.outer(aux2[0:-1, -1], aux2[-1,0:-1])/aux2[-1, -1]
    return Ahinv


#def MultiUpdateInverseHermitian(Binv, jrowMAT):
#    jrowMAT = np.sort(jrowMAT)
#    BinvOLD = Binv
#    for i in range(np.size(jrowMAT)):
#        jrow = jrowMAT[i] - i
#        Ahinv = UpdateInverseHermitian(BinvOLD, jrow)
#        BinvOLD = Ahinv
#    return Ahinv

def MultiUpdateInverseHermitian(H, jrowMAT):
    jrowMAT = np.sort(jrowMAT)
    for i in range(np.size(jrowMAT)):
        jrow = jrowMAT[i] - i
        H = UpdateInverseHermitian(H, jrow)
    return H


def ComputeROQ(Modes, weights, nGP, factorLEQ, tol):
    [J, bT] = ComputeJandb(Modes, weights, factorLEQ)[:2]
    b = bT.reshape(-1, 1)

    y = np.arange(len(weights))
    r = b  # residual vector, initial guess
    it = 0  # number of iterations
    mPOS = 0  # number of non-zero weights
    z = []
    Jnorm = np.sqrt(sum(np.multiply(J, J), 0))

    # point selection algorithm
    while (np.linalg.norm(r) / np.linalg.norm(b) > tol) and (mPOS <= nGP):
        # 1. Compute new point
        ObjFun = np.dot((J[:, y]).T, r)
        div = np.multiply(Jnorm[y], np.linalg.norm(r)).reshape(-1,1)
        ObjFun = np.divide(ObjFun, div)
        s = ObjFun.argmax()
        i = y[s]
        # 2. Update alpha and H (unrestricted least squares)
        if it == 0:
            # complies with newer versions of numpy
            #alpha = np.linalg.lstsq(J[:, [i]], b, rcond=None)[0]
            alpha = np.linalg.lstsq(J[:, [i]], b)[0]
            H = 1 / np.dot((J[:, i]).T, J[:, i])
        else:
            H, alpha = UpdateWeightsInverse(H, alpha, J[:, z], J[:, i], r)
        # 3. Move i from set y to set z
        z = (np.append(z, i)).astype(int)
        y = np.delete(y, s)
        # 4. Find possible negative weights
        if any(alpha < 0):
            print("WARNING: NEGATIVE weight found")
            indexes_neg_weight = np.where(alpha <= 0.)[0]
            y = np.append(y, (z[indexes_neg_weight]).T)
            z = np.delete(z, indexes_neg_weight)
            H = MultiUpdateInverseHermitian(H, indexes_neg_weight)
            alpha = H @ np.dot(J[:, z].T, b)

        # 6. Update the residual
        r = b - np.dot(J[:, z], alpha)
        # 7. Update mPOS and k
        mPOS = np.size(z)
        it = it + 1
        logger.debug("k = {}, mPOS = {}, error = {:.2f}%".format(it, mPOS, np.linalg.norm(r)/np.linalg.norm(b) * 100))
    # 6. Postprocess of points - neglecting null weights
    w = np.multiply(alpha, np.sqrt(weights[z]).reshape(-1, 1))
    logger.debug("Reduced Weights: {}".format(w.T))
    logger.debug("sum of reduced weights: {}".format(np.sum(w)))
    logger.debug("GP's index (elems and ip starting from zero): {}".format(z))
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

    bases_file_format = conf['Parameters']['bases_file_format']
    if bases_file_format == 'ascii':
        energy_modes = np.loadtxt(energy_bases_filename)[:,:nr_roq_points]
    else:
        energy_modes = np.load(energy_bases_filename)[:,:nr_roq_points]

    [w, z] = ComputeROQ(energy_modes, integration_weights,
                        nr_roq_points, factorLEQ=1.0, tol=1.e-14)
    roq_weigths = -1 * np.ones([nr_elements, nr_integration_points])
    roq_list = []
    for x, igg in enumerate(z):
        e = int(igg / nr_integration_points)
        ig = igg % nr_integration_points
        roq_weigths[e][ig] = w[x]
        roq_list.append([e, ig, w[x]])
    return roq_weigths, roq_list


def create_rom_weights(conf):
    integration_weights_filename = conf['Parameters']['integration_weights_filename']
    integration_weights = np.loadtxt(integration_weights_filename)
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])

    roq_list = []
    for x, ipw in enumerate(integration_weights):
        e = int(x / nr_integration_points)
        ip = x % nr_integration_points
        roq_list.append([e, ip, ipw])

    roq_weights = integration_weights.reshape((nr_elements, -1))

    return roq_weights, roq_list


def generate_rve_params(conf, iw_list):
    strain_bases_filename = conf['Parameters']['strain_bases_filename']
    nr_ip = int(conf['Parameters']['nr_integration_points'])
    nr_comps = int(conf['Parameters']['nr_strain_components'])
    nr_modes = int(conf['Parameters']['nr_active_modes'])
    rve_mdpa_filename = conf['Parameters']['rve_mdpa_filename']
    nr_dofs = nr_ip * nr_comps
    strain_bases = np.load(strain_bases_filename, mmap_mode='r')
    strain_bases = strain_bases[:,:nr_modes]

    # read model materials
    material = {}
    flag_elements = False
    with open(rve_mdpa_filename, 'r') as fi:
        for line in fi.readlines():
            if not flag_elements:
                if "Begin Elements" not in line:
                    continue
                else:
                    flag_elements = True
            else:
                if "End Elements" in line:
                    flag_elements = False
                    continue
                else:
                    material[int(line.split()[0]) - 1] = int(line.split()[1])
    out = {}
    out_B = []
    out_w = []
    out_prop = []
    B = np.empty((nr_comps, nr_modes))
    for list in iw_list:
        e = int(list[0])
        i = int(list[1])
        w = float(list[2])

        # get B
        index = e * nr_ip * nr_comps + i * nr_comps
        B = strain_bases[index:index + nr_comps, :]

        out_B.append(B.tolist())
        out_w.append(w)
        out_prop.append(material[e])

    out['props_id'] = out_prop
    out['w'] = out_w
    out['B'] = out_B

    return out

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
        roq_mask, roq_list = create_rom_weights(conf)
    else:
        logger.info("Computing HPROM")
        roq_mask, roq_list = compute_reduced_set(conf)
    logger.info("Generating RVE parameters for CL")
    rve_params = generate_rve_params(conf, roq_list)
    # output
    logging.debug("ROQ mask size {}".format(np.shape(roq_mask)))
    logging.debug("ROQ list size {}".format(np.shape(roq_list)))
    #filename = conf['Parameters']['roq_weights_filename']
    #np.savetxt(filename, roq_mask)
    #np.savetxt("roq_list.dat", roq_list)
    with open("rve.json", 'w') as fo:
        json.dump(rve_params, fo, indent=2)
