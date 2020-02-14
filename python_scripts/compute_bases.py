import os
import time
import glob
import logging
import h5py
import numpy

try:
    import sklearn.decomposition
except ImportError:
    pass


logging.basicConfig(
    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S", level=logging.DEBUG
)
logger = logging.getLogger(__name__)


def list_of_snapshots(trajectory_filename, value_name):
    logger.debug("Scanning trajectories")
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))
    e_arrays = []
    i_arrays = []
    for path in trajectory_paths:
        logger.debug("    trajectory {}".format(path))
        with h5py.File(path + "/" + "snapshots.hdf5", "r") as f:
            try:
                d = f["ELASTIC"][value_name]
                for k, v in d.items():
                    e_arrays.append(numpy.array(v))
            except KeyError:
                pass
            d = f["INELASTIC"][value_name]
            for k, v in d.items():
                i_arrays.append(numpy.array(v))
    logger.debug("Nr of elastic arrays loaded: {}".format(len(e_arrays)))
    logger.debug("Nr of inelastic arrays loaded: {}".format(len(i_arrays)))
    return e_arrays, i_arrays


def compute_modes(list_of_arrays, nr_modes, Ue=None, svd_algorithm="standard"):
    logger.info("Loading snapshots")
    if Ue is not None:
        logger.info("and removing elastic component")
    X = numpy.empty([len(list_of_arrays[0]), len(list_of_arrays)])
    total = len(list_of_arrays)
    batch_size = int(total / 10 + 0.5)
    counter = 1
    for i, array in enumerate(list_of_arrays):
        X[:, i] = array
        if Ue is not None:
            for j in range(Ue.shape[1]):
                X[:, i] -= numpy.multiply(numpy.dot(Ue[:, j], X[:, i]), Ue[:, j])
        if not counter % batch_size:
            logger.info("    {}/{} snapshots processed".format(counter, total))
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
    elif "standard" in svd_algorithm:
        t0 = time.time()
        logger.info("Computing SVD using STANDARD algorithm")
        [U, S] = numpy.linalg.svd(X, full_matrices=False)[:2]
        U = U[:, :nr_modes]
        logger.info("SVD time: {:.1f}s".format(time.time() - t0))

    logger.info("    - SVD time: {:.1f}s".format(time.time() - t0))
    logger.info("    - singular value of selected modes:")
    logger.info("      {}".format(S[:nr_modes]))
    logger.info("      validation: following singular values (excluded):")
    logger.info("      {}".format(S[nr_modes : nr_modes + 4]))
    logger.info("    - nr and size of modes: {}, {}".format(U.shape[1], U.shape[0]))
    logger.info("")
    numpy.savetxt("singular_values.dat", S)

    return U


def generate_bases(
    nr_elastic_modes,
    nr_inelastic_modes,
    e_arrays,
    i_arrays,
    bases_fname,
    svd_algorithm="standard",
):
    t0 = time.time()
    if nr_elastic_modes > 0:
        logger.info("Processing elastic snapshots")
        Ue = compute_modes(e_arrays, nr_elastic_modes, svd_algorithm=svd_algorithm)
        os.rename("singular_values.dat", "singular_values_elastic.dat")
        logger.info("Processing inelastic snapshots")

        Ui = compute_modes(
            i_arrays, nr_inelastic_modes, Ue=Ue, svd_algorithm=svd_algorithm
        )
        os.rename("singular_values.dat", "singular_values_inelastic.dat")
        U = numpy.hstack([Ue, Ui])

    # Needs modificacion of e_files + i_files to work
    else:
        logger.info(
            "Nr of elastic modes set to zero -> Not discriminating elastic/inelastic snapshots"
        )
        U = compute_modes(i_arrays, nr_inelastic_modes, svd_algorithm=svd_algorithm)
    t1 = time.time()
    numpy.save(bases_fname, U)
    logger.info("  Writing time: {:.1f}s".format(time.time() - t1))
    logger.info("  Total time: {:.1f}s".format(time.time() - t0))
    logger.info("")
