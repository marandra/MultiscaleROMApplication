import os
import time
import glob
import numpy as np
import scipy.sparse.linalg as sp
import logging

try:
    import sklearn.decomposition
except ImportError:
    pass


logging.basicConfig(
    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S", level=logging.DEBUG
)
logger = logging.getLogger(__name__)


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


def compute_modes(
    nr_integration_points,
    files,
    nr_modes,
    nr_components,
    Ue=None,
    svd_algorithm="standard",
):
    logger.info("Loading snapshots")
    if Ue is not None:
        logger.info("and removing elastic component")
    nr_dofs = nr_integration_points * nr_components
    X = np.empty([nr_dofs, len(files)])
    total = len(files)
    batch_size = int(total / 10 + 0.5)
    counter = 1
    for i, file in enumerate(files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        if Ue is not None:
            for j in range(Ue.shape[1]):
                X[:, i] -= np.multiply(np.dot(Ue[:, j], X[:, i]), Ue[:, j])
        if not counter % batch_size:
            logger.debug("    {}/{} snapshots processed".format(counter, total))
        counter = counter + 1
    logger.info("")

    # SVD stage  # svd_algorithm = standard, iterative, arpack, randomized, auto?
    t0 = time.time()
    if "randomized" in svd_algorithm:
        t0 = time.time()
        logger.info("Computing SVD using RANDOMIZED algorithm")
        svd = sklearn.decomposition.TruncatedSVD(
            n_components=nr_modes, algorithm="randomized"
        )
        svd.fit(X.T)
        U = svd.components_.T
        S = svd.singular_values_.T
        logger.info("SVD time: {:.1f}s".format(time.time() - t0))
    #    #elif svd_algorithm is "arpack":
    #        t0 = time.time()
    #        logger.info("Computing SVD using ARPACK algorithm")
    #        svd = sklearn.decomposition.TruncatedSVD(n_components=nr_modes, algorithm="arpack")
    #        svd.fit(X.T)
    #        U = svd.components_.T
    #        S = svd.singular_values_.T
    #        print("DEBUG:")
    #        print("SVD time: {:.1f}s".format(time.time() - t0))
    #        print(svd)
    #        #print(U)
    #        #print(S)
    #    if True:
    #    #elif svd_algorithm is "iterative":
    #        t0 = time.time()
    #        logger.info("Computing SVD using ITERATIVE algorithm")
    #        #[U, S] = sp.svds(X, k=nr_modes + 4)[:2]
    #        [U, S] = sp.svds(X, k=nr_modes)[:2]
    #        # to order values in decreasing order (svds returns them in increasing order)
    #        S= S[::-1]
    #        U = U[:,::-1]
    #        U = U[:, :nr_modes]
    #        print("DEBUG:")
    #        print("SVD time: {:.1f}s".format(time.time() - t0))
    #        #print(U)
    #        #print(S)
    elif "standard" in svd_algorithm:
        t0 = time.time()
        logger.info("Computing SVD using STANDARD algorithm")
        [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
        U = U[:, :nr_modes]
        logger.info("SVD time: {:.1f}s".format(time.time() - t0))

    logger.info("    - SVD time: {:.1f}s".format(time.time() - t0))
    logger.info("    - singular value of selected modes:")
    logger.info("      {}".format(S[:nr_modes]))
    #logger.info("      validation: following singular values (excluded):")
    #logger.info("      {}".format(S[nr_modes: nr_modes + 4]))
    logger.info("    - nr and size of modes: {}, {}".format(U.shape[1], U.shape[0]))
    logger.info("")
    np.savetxt("singular_values.dat", S)

    return U


def generate_bases(
    nr_ips,
    nr_strain_components,
    nr_elastic_modes,
    nr_inelastic_modes,
    e_files,
    i_files,
    bases_fname,
    svd_algorithm="standard",
):
    t0 = time.time()
    if nr_elastic_modes > 0:
        logger.info("Processing elastic snapshots")
        Ue = compute_modes(
            nr_ips,
            e_files,
            nr_elastic_modes,
            nr_components=nr_strain_components,
            svd_algorithm=svd_algorithm,
        )
        os.rename("singular_values.dat", "singular_values_elastic.dat")
        logger.info("Processing inelastic snapshots")
        Ui = compute_modes(
            nr_ips,
            i_files,
            nr_inelastic_modes,
            nr_components=nr_strain_components,
            Ue=Ue,
            svd_algorithm=svd_algorithm,
        )
        os.rename("singular_values.dat", "singular_values_inelastic.dat")
        U = np.hstack([Ue, Ui])
    else:
        logger.info(
            "Nr of elastic modes set to zero -> Not discriminating elastic/inelastic snapshots"
        )
        U = compute_modes(
            nr_ips,
            e_files + i_files,
            nr_inelastic_modes,
            nr_components=nr_strain_components,
            svd_algorithm=svd_algorithm,
        )
    t1 = time.time()
    np.save(bases_fname, U)
    logger.info("  Writing time: {:.1f}s".format(time.time() - t1))
    logger.info("  Total time: {:.1f}s".format(time.time() - t0))
    logger.info("")
