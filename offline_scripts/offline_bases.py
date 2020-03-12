import os
import numpy
import h5py
import time
import glob
import logging
import sklearn.decomposition
from offline_common import Common


"""
TODO: pending description here.
"""

logging.basicConfig(
    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S", level=logging.DEBUG
)
fh = logging.FileHandler("offline_bases.log")
# fh.setLevel(logging.DEBUG)
fh.setFormatter(logging.Formatter("[%(asctime)s] %(message)s"))
# ch = logging.StreamHandler()
# ch.setLevel(logging.DEBUG)
# ch.setFormatter(logging.Formatter("[%(asctime)s] %(message)s"))
logger = logging.getLogger(__name__)
logger.addHandler(fh)
# logger.addHandler(ch)


def skip_calculation(filename, flag_reuse):
    try:
        with open(filename):
            flag_exists = True
    except IOError:
        flag_exists = False
    return flag_exists and flag_reuse


def read_snapshots(trajectory_filename, group, field):
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))

    # count snapshots so we know the size of X
    logger.info("Getting number and size of snapshots to allocate array")
    nr_snapshots = 0
    for path in trajectory_paths:
        with h5py.File(path + "/" + "snapshots.hdf5", "r") as f:
            try:
                d = f[group][field]
                for k, v in d.items():
                    nr_snapshots += 1
            except KeyError:
                pass
    for path in trajectory_paths:
        with h5py.File(path + "/" + "snapshots.hdf5", "r") as f:
            try:
                d = f[group][field]
                for k, v in d.items():
                    len_snapshot = len(v)
            except KeyError:
                pass
    logger.info("    - {} snapshots size {}".format(nr_snapshots, len_snapshot))

    # start loading snapshots
    logger.info("Loading snapshots")
    arrays = numpy.empty([len_snapshot, nr_snapshots])
    batch_size = int(len(trajectory_paths) / 10 + 0.5)
    counter = 1
    column = 0
    for path in trajectory_paths:
        with h5py.File(path + "/" + "snapshots.hdf5", "r") as f:
            try:
                d = f[group][field]
                for k, v in d.items():
                    arrays[:, column] = v
                    column += 1
            except KeyError:
                logger.debug(
                    "Skipping {}/{} of {} (dataset not present)".format(group, field, path)
                )
        if not counter % batch_size:
            logger.info(
                "    {}/{} trajectories processed".format(
                    counter, len(trajectory_paths)
                )
            )
        counter += 1
    return arrays


def remove_elastic_modes(X, Ue):
    logger.info("Removing elastic component")
    for i in range(numpy.shape(X)[1]):
        for j in range(Ue.shape[1]):
            X[:, i] -= numpy.multiply(numpy.dot(Ue[:, j], X[:, i]), Ue[:, j])
    return X


def compute_svd(X, nr_modes, svd_algorithm="standard"):
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
    elif "standard" in svd_algorithm:
        t0 = time.time()
        logger.info("Computing SVD using STANDARD algorithm")
        [U, S] = numpy.linalg.svd(X, full_matrices=False)[:2]
        U = U[:, :nr_modes]

    logger.info("    - SVD time: {:.1f}s".format(time.time() - t0))
    logger.info("    - singular value of selected modes:")
    logger.info("      {}".format(S[:nr_modes]))
    logger.info("      validation: following singular values (excluded):")
    logger.info("      {}".format(S[nr_modes : nr_modes + 4]))
    logger.info("    - nr and size of modes: {}, {}".format(U.shape[1], U.shape[0]))
    logger.info("")
    numpy.savetxt("singular_values.dat", S)
    return U


def create_bases(
    field_name,
    nr_elastic_modes,
    nr_inelastic_modes,
    trajectory_filename,
    bases_fname,
    svd_algorithm,
    reuse_files,
):
    if skip_calculation(bases_fname, reuse_files):
        logger.info("File {} exists. Skipping calculation".format(bases_fname))
        return
    logger.info("Generating {} bases".format(field_name))

    t0 = time.time()
    # Snapshots splitted in elastic and inelastic groups
    if nr_elastic_modes > 0:
        logger.info("Processing elastic snapshots")
        X = read_snapshots(trajectory_filename, "ELASTIC", field_name)
        Ue = compute_svd(X, nr_elastic_modes, svd_algorithm=svd_algorithm)
        os.rename(
            "singular_values.dat",
            "sv_{}_elastic_{}.dat".format(field_name, nr_elastic_modes),
        )

        logger.info("Processing inelastic snapshots")
        X = read_snapshots(trajectory_filename, "INELASTIC", field_name)
        # backup inelastic snapshot in case we run out of memory
        snapshots_fname = "auxiliar_snapshots_{}_{}.npy".format("INELASTIC", field_name)
        numpy.save(snapshots_fname, X)
        X = remove_elastic_modes(X, Ue)
        # backup inelastic snapshot in case we run out of memory
        snapshots_fname = "auxiliar_snapshots_{}_{}_removed_elastic.npy".format("INELASTIC", field_name)
        numpy.save(snapshots_fname, X)
        try:
            Ui = compute_svd(X, nr_inelastic_modes, svd_algorithm=svd_algorithm)
        except MemoryError:
            # We don't have enough RAM for SVD. We free up memory moving snapshots array to disk
            logger.warning("Run out of memory. Reloading snapshots from disk as memory map.")
            del X
            X = numpy.load(snapshots_fname, mmap_mode="r")
            Ui = compute_svd(X, nr_inelastic_modes, svd_algorithm=svd_algorithm)
        os.rename(
            "singular_values.dat",
            "sv_{}_inelastic_{}.dat".format(field_name, nr_inelastic_modes),
        )

        U = numpy.hstack([Ue, Ui])

    # No splitting of elastic and inelastic snapshots
    else:
        logger.info(
            "Nr of elastic modes set to zero -> "
            "Not discriminating elastic/inelastic snapshots"
        )
        X = read_snapshots(trajectory_filename, "INELASTIC", field_name)
        U = compute_svd(X, nr_inelastic_modes, svd_algorithm=svd_algorithm)
        os.rename(
            "singular_values.dat",
            "sv_{}_{}.dat".format(field_name, nr_elastic_modes + nr_inelastic_modes),
        )

    numpy.save(bases_fname, U)
    logger.info("  Elapsed time: {:.1f}s".format(time.time() - t0))
    logger.info("")


#######################################################################
#######################################################################

if __name__ == "__main__":

    logger.info("Begining bases calculation -----------------------")

    case_basename = "../training/trajectory"

    #
    # compute energy bases
    #
    create_bases(
        Common().energy_name,
        Common().energy_elastic_modes,
        Common().energy_inelastic_modes,
        case_basename,
        Common().energy_bases_fname,
        Common().svd_algorithm,
        Common().reuse_existing_files,
    )

    #
    # compute strain bases
    #
    create_bases(
        Common().strain_name,
        Common().strain_elastic_modes,
        Common().strain_inelastic_modes,
        case_basename,
        Common().strain_bases_fname,
        Common().svd_algorithm,
        Common().reuse_existing_files,
    )

    #
    # compute R_VALUE bases
    #
    create_bases(
        Common().rvalue_name,
        Common().rvalue_elastic_modes,
        Common().rvalue_inelastic_modes,
        case_basename,
        Common().rvalue_bases_fname,
        Common().svd_algorithm,
        Common().reuse_existing_files,
    )

    logger.info("Finished -----------------------------------------")
