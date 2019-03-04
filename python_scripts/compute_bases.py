import time
import glob
import numpy as np
import scipy.sparse.linalg as sp
import logging


def list_of_snapshots(trajectory_filename, nr_e_snap_filename, filename):
    logger.debug("Scanning trajectories")
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))
    e_files = []
    i_files = []
    for path in trajectory_paths:
        with open("{}/{}".format(path, nr_e_snap_filename), "r") as f:
            nr_e_snap = int(f.readline().strip())
        logger.debug("  {} - elastic snapshots: {}".format(path, nr_e_snap))
        e_files.extend(sorted(glob.glob("{}/{}*".format(path, filename)))[:nr_e_snap])
        i_files.extend(sorted(glob.glob("{}/{}*".format(path, filename)))[nr_e_snap:])
    logger.debug("")
    return e_files, i_files


def compute_modes(nr_integration_points, nr_elements, files, nr_modes, nr_components, Ue=None):
    logger.info("Loading snapshots")
    if Ue is not None:
        logger.info("and removing elastic component")
    nr_dofs = nr_elements * nr_integration_points * nr_components
    X = np.empty([nr_dofs, len(files)])
    total = len(files)
    batch_size = int(total / 10 + .5)
    counter = 1
    for i, file in enumerate(files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        if Ue is not None:
            for j in range(Ue.shape[1]):
                X[:,i] -= np.multiply(np.dot(Ue[:, j], X[:, i]), Ue[:, j])
        if not counter % batch_size:
            logger.debug("    {}/{} snapshots processed".format(counter, total))
        counter = counter + 1
    logger.info("")

    #SVD stage
    svd_iterative = True
    if svd_iterative:
        logger.info("Computing SVD using ITERATIVE algorithm")
        #[U, S] = sp.svds(X, k=nr_modes + 4)[:2]
        [U, S] = sp.svds(X, k=nr_modes)[:2]
        # to order values in decreasing order (svds returns them in increasing order)
        S= S[::-1]
        U = U[:,::-1]
        U = U[:, :nr_modes]
    else:
        logger.info("Computing SVD using standard algorithm")
        [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
        U = U[:,:nr_modes]

    logger.info("    - singular value of selected modes:")
    logger.info("      {}".format(S[:nr_modes]))
    logger.info("      validation: following singular values (excluded):")
    logger.info("      {}".format(S[nr_modes: nr_modes + 4]))
    logger.info("    - nr and size of modes: {}, {}".format(U.shape[1], U.shape[0]))
    logger.info("")
    # temp change default print option, as we want to have all the values printed out
    #np.set_printoptions(threshold=np.inf)
    #logger.debug("    - all singular values:")
    #logger.debug("      {}".format(S))
    #logger.debug("")
    #np.set_printoptions(threshold=1000)

    return U

def generate_bases(nr_elements, nr_ip, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, bases_fname):
        t0 = time.time()
        logger.info("Generating bases ENERGY")
        if nr_elastic_modes > 0:
            logger.info("Processing elastic snapshots")
            Ue = compute_modes(nr_ip, nr_elements, e_files, nr_elastic_modes, nr_components=nr_strain_components)
            logger.info("Processing inelastic snapshots")
            Ui = compute_modes(nr_ip, nr_elements, i_files, nr_inelastic_modes, nr_components=nr_strain_components, Ue=Ue)
            U = np.hstack([Ue, Ui])
        else:
            logger.info("Nr of elastic modes set to zero -> Not discriminating elastic/inelastic snapshots")
            U = compute_modes(nr_ip, nr_elements, e_files + i_files, nr_inelastic_modes, nr_components=nr_strain_components)
        t1 = time.time()
        np.save(bases_fname, U)
        logger.info("  SVD time: {:.1f}s".format(time.time() - t0))
        logger.info("  Writing time: {:.1f}s".format(time.time() - t1))
        logger.info("")

#######################################
# Main
#######################################

logging.basicConfig(format='[%(asctime)s] %(message)s',
                    datefmt='%H:%M:%S', level=logging.DEBUG)
logger = logging.getLogger(__name__)

if __name__ == '__main__':

    trajectory_filename = conf['Parameters']['trajectory_filename']
    nr_e_snap_filename = conf['Parameters']['nr_elastic_snapshots_filename']

    if flag_comp_energy:
        nr_elements = 320
        nr_ip = 8
        nr_strain_components = 1
        nr_elastic_modes = 21
        nr_inelastic_modes = 100
        bases_fname = "bases_energy.npy"
        snapshot_filename = conf['Parameters']['energy_filename']
        e_files, i_files = list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
        generate_bases(nr_elements, nr_ip, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, bases_fname)

    if flag_comp_strain:
        nr_elements = 320
        nr_ip = 8
        nr_strain_components = 6
        nr_elastic_modes = 6
        nr_inelastic_modes = 100
        bases_fname = "bases_strain.npy"
        snapshot_filename = conf['Parameters']['strain_filename']
        e_files, i_files = list_of_snapshots(trajectory_filename, nr_e_snap_filename, snapshot_filename)
        generate_bases(nr_elements, nr_ip, nr_strain_components, nr_elastic_modes, nr_inelastic_modes, e_files, i_files, bases_fname)
