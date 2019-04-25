import configparser
import argparse
import numpy as np
import logging


def remove_exact_integral_energy(modes, weights):
    eps = np.finfo(float).eps # 2.22044604925e-16
    # Total microscale volume
    total_weight = np.sum(weights)
    sqrt_total_weight = np.sqrt(total_weight)
    sqrt_weights = np.sqrt(weights)
    # Normalized exact integral
    norm_exact_integral = modes.T @ weights / total_weight
    # Matrix of modified modes (with zero integral)
    modified_modes = (modes - norm_exact_integral) * sqrt_weights.reshape(-1,1)
    [modified_modes, bases_weights] = np.linalg.svd(modified_modes, full_matrices=False)[:2]
    # filter the reduced modified set of modes
    tolerance = np.max(modes.shape) * eps * np.max(bases_weights)
    rank_mod_modes = sum(i > tolerance for i in bases_weights)
    modified_modes = modified_modes[:, 0:rank_mod_modes]
    # Adding last row related with the sqrt of gauss integration weigths
    # and initializing the RHS vector for the optimization problem
    J = np.vstack([modified_modes.T, (sqrt_weights / sqrt_total_weight).T])
    b = np.vstack([np.zeros((modified_modes.T.shape[0], 1)), sqrt_total_weight])
    return J, b


def update_weights_inverse(H, alpha, bases_current, base_new, r):
    c = np.dot(bases_current.T, base_new)
    d = np.dot(H, c).reshape(-1, 1)
    s = np.dot(base_new.T, base_new) - np.dot(c.T, d)
    aux1 = np.hstack([H + np.outer(d, d)/s, -d/s])
    aux2 = np.hstack([np.squeeze(-d.T/s), 1/s])
    H_new = np.vstack([aux1, aux2])
    v = np.dot(base_new.T, r) / s
    alpha = np.vstack([(alpha - d * v), v])
    return H_new, alpha


def update_inverse_hermitian(invH, neg_index):
    if neg_index == np.shape(invH)[1]:
        aux = (invH[0:-1, -1] * invH[-1, 0:-1]) / invH(-1, -1)
        invH_new = invH[:-1, :-1] - aux
    else:
        aux1 = np.hstack([invH[:, 0:neg_index], invH[:, neg_index + 1:], invH[:, neg_index].reshape(-1, 1)])
        aux2 = np.vstack([aux1[0:neg_index, :], aux1[neg_index + 1:, :], aux1[neg_index, :]])
        invH_new = aux2[0:-1, 0:-1] - np.outer(aux2[0:-1, -1], aux2[-1,0:-1])/aux2[-1, -1]
    return invH_new


def multiupdate_inverse_hermitian(invH, neg_indexes):
    neg_indexes = np.sort(neg_indexes)
    for i in range(np.size(neg_indexes)):
        neg_index = neg_indexes[i] - i
        invH = update_inverse_hermitian(invH, neg_index)
    return invH


def compute_roq(Modes, weights, nGP, tol):
    J, b = remove_exact_integral_energy(Modes, weights)
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
            H, alpha = update_weights_inverse(H, alpha, J[:, z], J[:, i], r)
        # 3. Move i from set y to set z
        z = (np.append(z, i)).astype(int)
        y = np.delete(y, s)
        # 4. Find possible negative weights
        if any(alpha < 0):
            print("WARNING: NEGATIVE weight found")
            indexes_neg_weight = np.where(alpha <= 0.)[0]
            y = np.append(y, (z[indexes_neg_weight]).T)
            z = np.delete(z, indexes_neg_weight)
            H = multiupdate_inverse_hermitian(H, indexes_neg_weight)
            alpha = np.dot(H, np.dot(J[:, z].T, b))

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


def compute_hprom_weights(nr_elemental_ip, integration_weights, nr_roq_points, energy_bases_filename):
    logger.info("Computing reduced set of integration points (HPROM)")
    energy_modes = np.load(energy_bases_filename)[:,:nr_roq_points]

    [w, z] = compute_roq(energy_modes, np.array(integration_weights),
                         nr_roq_points, tol=1.e-14)
    roq_list = []
    for x, igg in enumerate(z):
        e = int(igg / nr_elemental_ip)
        ig = igg % nr_elemental_ip
        roq_list.append([e, ig, w[x][0], igg])
    return roq_list


def compute_rom_weights(nr_elemental_ip, integration_weights):
    logger.info("Computing complete set of integration points (ROM)")

    roq_list = []
    for x, ipw in enumerate(integration_weights):
        e = int(x / nr_elemental_ip )
        ip = x % nr_elemental_ip
        roq_list.append([e, ip, ipw])
    return roq_list


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Computes Reduced Order Quadrature (ROQ) integration weights")
#parser.add_argument('config_file', help="configuration file")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
parser.add_argument('-r', '--rom', action="store_true", help="compute ROM instead of HPROM")
args = parser.parse_args()

# parse configuration file
#conf = configparser.ConfigParser()
#conf.read(args.config_file)

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s',
                    datefmt='%H:%M:%S',level=verbosity_level)
logger = logging.getLogger(__name__)
#handler = logging.FileHandler('log_' + args.config_file.rsplit('.', 1)[0])
#handler.setLevel(logging.DEBUG)
#logger.addHandler(handler)

if __name__ == '__main__':
    logger.info("Reduced Order Quadrature")
    nr_ip_per_element = 8
    integration_weights = np.loadtxt("integration_weight")
    nr_roq_points = 50
    energy_bases_filename = 'bases_energy.npy'
    if args.rom:
        roq_list = compute_rom_weights(nr_ip_per_element, integration_weights)
    else:
        roq_list = compute_hprom_weights(nr_ip_per_element, integration_weights, nr_roq_points, energy_bases_filename)

    logging.debug("ROQ list size {}".format(np.shape(roq_list)))
    np.savetxt("roq_list.dat", roq_list)
