import os
import numpy
import h5py
import time
import glob
import logging
import sklearn.decomposition
from offline_common import Common
from pathlib import Path
import multiprocessing


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


def get_shape_of_snapshots(path, group, field):
    ns = 0
    ls = 0
    fpath = Path(path, "snapshots.hdf5")
    with h5py.File(fpath, "r") as f:
        try:
            d = f[group][field]
            for k, v in d.items():
                ns += 1
                ls = len(v)  # TODO not optimal. read only once
        except KeyError:
            pass
    return ns, ls


def get_snapshots(path, group, field):
    ns, ls = get_shape_of_snapshots(path, group, field)
    arrays = numpy.empty([ls, ns])
    column = 0
    fpath = Path(path, "snapshots.hdf5")
    with h5py.File(fpath, "r") as f:
        try:
            d = f[group][field]
            for k, v in d.items():
                arrays[:, column] = v
                column += 1
        except KeyError:
            logger.debug(
                "Skipping {}/{} of {} (dataset not present)".format(group, field, path)
            )
    return arrays


def read_local_svd(trajectory_filename, field, cutoff_tol):
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))

    logger.info("Getting size of local svd")

    nr_modes = 0
    for path in trajectory_paths:
        a = numpy.load(
            path + "/" + "bases_inelastic_local_{}.npy".format(field), mmap_mode="r"
        )
        sv = numpy.loadtxt(path + "/" + "sv_inelastic_local_{}".format(field))
        idx = numpy.where(sv > cutoff_tol)[0]
        nr_col = len(idx)
        nr_modes += nr_col
        len_mode = numpy.shape(a)[0]
    logger.info("    - {} modes size {}".format(nr_modes, len_mode))

    # start loading bases
    logger.info("Loading local bases")
    arrays = numpy.empty([len_mode, nr_modes])
    batch_size = int(len(trajectory_paths) / 10 + 0.5)
    counter = 1
    column = 0
    for path in trajectory_paths:
        array = numpy.load(
            path + "/" + "bases_inelastic_local_{}.npy".format(field), mmap_mode="r"
        )
        sv = numpy.loadtxt(path + "/" + "sv_inelastic_local_{}".format(field))
        idx = numpy.where(sv > cutoff_tol)[0]
        nr_col = len(idx)
        arrays[:, column : column + nr_col] = array[:, idx] * sv[idx]
        column += nr_col
        #
        if not counter % batch_size:
            logger.info(
                "    {}/{} trajectories processed".format(
                    counter, len(trajectory_paths)
                )
            )
        counter += 1
        #
    return arrays


def read_snapshots(trajectory_filename, group, field):
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))

    # count snapshots so we know the size of X
    logger.info("Getting number and size of snapshots to allocate array")

    nr_snapshots = 0
    for path in trajectory_paths:
        ns, len_snapshot = get_shape_of_snapshots(path, group, field)
        nr_snapshots += ns
    logger.info("    - {} snapshots size {}".format(nr_snapshots, len_snapshot))

    # start loading snapshots
    logger.info("Loading snapshots")
    arrays = numpy.empty([len_snapshot, nr_snapshots])
    batch_size = int(len(trajectory_paths) / 10 + 0.5)
    counter = 1
    column = 0
    for path in trajectory_paths:
        array = get_snapshots(path, group, field)
        arrays[:, column : column + numpy.shape(array)[1]] = array
        column += numpy.shape(array)[1]
        #
        if not counter % batch_size:
            logger.info(
                "    {}/{} trajectories processed".format(
                    counter, len(trajectory_paths)
                )
            )
        counter += 1
        #
    return arrays


def remove_elastic_modes(X, Ue):
    logger.info("Removing elastic componennt")
    t0 = time.time()
    for i in range(numpy.shape(X)[1]):
        projection = X[:, i] @ Ue
        X[:, i] -= numpy.sum(projection * Ue, axis=1)
    logger.debug("    elapsed time: {:.1f}s".format(time.time() - t0))
    return X


def compute_svd(X, nr_modes, svd_algorithm="standard"):
    # SVD stage  # svd_algorithm = standard, iterative, arpack, randomized, auto?
    t0 = time.time()
    if nr_modes > -1:
        logger.info("Computing SVD using RANDOMIZED algorithm")
        svd = sklearn.decomposition.TruncatedSVD(
            n_components=nr_modes, algorithm="randomized"
        )
        svd.fit(X.T)
        U = svd.components_.T
        S = svd.singular_values_.T
    else:
        logger.info("Computing SVD using STANDARD algorithm")
        [U, S] = numpy.linalg.svd(X, full_matrices=False)[:2]
        # U = U[:, :nr_modes]

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
    reuse_files,
    cutoff_tol,
):
    if Common().skip_calculation(bases_fname.format(field_name, "*"), reuse_files):
        logger.info(
            "File {} exists. Skipping calculation".format(
                bases_fname.format(field_name, "*")
            )
        )
        return
    logger.info("Generating {} bases".format(field_name))

    t0 = time.time()
    # Snapshots splitted in elastic and inelastic groups
    if nr_elastic_modes > 0:
        logger.info("Processing elastic snapshots")
        X = read_snapshots(trajectory_filename, "ELASTIC", field_name)
        Ue = compute_svd(X, nr_elastic_modes)
        os.rename(
            "singular_values.dat", "sv_{}_elastic.dat".format(field_name),
        )

        logger.info("Processing inelastic LOCAL MODES")
        X = read_local_svd(trajectory_filename, field_name, cutoff_tol)
        X = remove_elastic_modes(X, Ue)
        Ui = compute_svd(X, nr_inelastic_modes)
        os.rename(
            "singular_values.dat", "sv_{}_inelastic.dat".format(field_name),
        )

        U = numpy.hstack([Ue, Ui])

    # No splitting of elastic and inelastic snapshots
    else:
        logger.info(
            "Nr of elastic modes set to zero -> "
            "Not discriminating elastic/inelastic snapshots"
        )
        X = read_snapshots(trajectory_filename, "INELASTIC", field_name)
        U = compute_svd(X, nr_inelastic_modes)
        os.rename(
            "singular_values.dat", "sv_{}.dat".format(field_name),
        )

    numpy.save(bases_fname.format(field_name, numpy.shape(U)[1]), U)
    logger.info("  Elapsed time: {:.1f}s".format(time.time() - t0))
    logger.info("")


def generate_local_bases(case_path, field):
    logger.debug("   - missing {} {}".format(case_path, field))
    case_path = Path(case_path)  # we can accept strings or Paths
    base_path = case_path / Common().local_bases_fname.format(field)
    X = get_snapshots(case_path, "INELASTIC", field)
    [U, S] = numpy.linalg.svd(X, full_matrices=False)[:2]
    numpy.save(base_path, U)
    path = Path(case_path, Common().local_sv_fname.format(field))
    numpy.savetxt(path, S)


def generate_missing_local_bases(training_path, field, threads=1):
    logger.info("Looking for missing local bases {}".format(field))
    cases_path = Path(training_path).glob("trajectory_*")
    missing = []
    for case in cases_path:
        bases = case / Common().local_bases_fname.format(field)
        if Common().skip_calculation(bases, Common().reuse_existing_files):
            continue
        missing.append(case)
    if not missing:
        return
    # There are missing bases files. Let's generate them.

    # Testing: version with Pool
    with multiprocessing.Pool(processes=threads) as pool:
        logger.debug("   - generating bases")
        logger.debug("   - multiprocessing {} threads".format(threads))
        pool.starmap(generate_local_bases, zip(missing, [field] * len(missing)))

    # Testing: version with Process
    # processes = []
    # for case in missing:
    #    semaphore.acquire()
    #    p = multiprocessing.Process(target=generate_local_bases, args=(case, field, semaphore))
    #    processes.append(p)
    #    p.start()
    # for p in processes:
    #    p.join()


#######################################################################
#######################################################################

if __name__ == "__main__":

    logger.info("Beginning bases calculation -----------------------")

    case_basename = "../training"

    #
    # generate missing local bases
    #
    generate_missing_local_bases(
        case_basename, Common().energy_name,
    )
    generate_missing_local_bases(
        case_basename, Common().strain_name,
    )
    generate_missing_local_bases(
        case_basename, Common().rvalue_name,
    )

    #
    # compute bases
    #
    create_bases(
        Common().energy_name,
        Common().energy_elastic_modes,
        Common().energy_inelastic_modes,
        case_basename,
        Common().bases_fname,
        Common().reuse_existing_files,
        Common().svd_cutoff[Common().energy_name],
    )
    create_bases(
        Common().strain_name,
        Common().strain_elastic_modes,
        Common().strain_inelastic_modes,
        case_basename,
        Common().bases_fname,
        Common().reuse_existing_files,
        Common().svd_cutoff[Common().strain_name],
    )
    create_bases(
        Common().rvalue_name,
        Common().rvalue_elastic_modes,
        Common().rvalue_inelastic_modes,
        case_basename,
        Common().bases_fname,
        Common().reuse_existing_files,
        Common().svd_cutoff[Common().rvalue_name],
    )

    logger.info("Finished -----------------------------------------")
