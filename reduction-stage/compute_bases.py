import time
import configparser
import argparse
import glob
import numpy as np
import logging


def write_bases(filename, U):
    bases_file_format = conf['Parameters']['bases_file_format']
    if bases_file_format == 'ascii':
        np.savetxt(filename, U)
    else:
        np.save(filename, U)
    return


def make_list_of_files(conf):
    logger.debug("Scanning trajectories")
    trajectory_filename = conf['Parameters']['trajectory_filename']
    nr_e_snap_filename = conf['Parameters']['nr_elastic_snapshots_filename']
    energy_filename = conf['Parameters']['energy_filename']
    strain_filename = conf['Parameters']['strain_filename']
    trajectory_paths = sorted(glob.glob("{}_?".format(trajectory_filename)))
    ene_e_files = []
    ene_i_files = []
    str_e_files = []
    str_i_files = []
    for path in trajectory_paths:
        with open("{}/{}".format(path, nr_e_snap_filename), "r") as f:
            nr_e_snap = int(f.readline().strip())
        logger.debug("  {} - elastic snapshots: {}".format(path, nr_e_snap))
        ene_e_files.extend(sorted(glob.glob("{}/{}*".format(path, energy_filename)))[:nr_e_snap])
        ene_i_files.extend(sorted(glob.glob("{}/{}*".format(path, energy_filename)))[nr_e_snap:])
        str_e_files.extend(sorted(glob.glob("{}/{}*".format(path, strain_filename)))[:nr_e_snap])
        str_i_files.extend(sorted(glob.glob("{}/{}*".format(path, strain_filename)))[nr_e_snap:])
    logger.debug("")
    return ene_e_files, ene_i_files, str_e_files, str_i_files


def compute_inelastic_modes(conf, files, Ue, nr_components):
    logger.info("  Projection of inelastic snapshots")
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_modes = int(conf['Parameters']['max_nr_reduced_modes'])
    nr_dofs = nr_elements * nr_integration_points * nr_components
    X = np.empty([nr_dofs, len(files)])
    for i, file in enumerate(files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        for j in range(Ue.shape[1]):
            X[:,i] = X[:,i] - np.multiply(np.dot(Ue[:, j], X[:, i]), Ue[:, j])
    logger.info("  SVD of inelastic snapshots")
    U = np.linalg.svd(X, full_matrices=False)[0]
    Ur = U[:,:nr_modes]
    logger.info("  - nr of modes: {}".format(Ur.shape[1]))
    logger.info("  - size of mode: {}".format(Ur.shape[0]))
    return Ur


def compute_elastic_modes(conf, files, nr_components):
    logger.info("  SVD of elastic snapshots")
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_strain_components = int(conf['Parameters']['nr_strain_components'])
    nr_modes = 2 * nr_strain_components
    nr_dofs = nr_elements * nr_integration_points * nr_components
    X = np.empty([nr_dofs, len(files)])
    for i, file in enumerate(files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    Ur = U[:,:nr_modes]
    logger.debug("  - singular values of selected modes: {}".format(S[:nr_modes]))
    logger.info("  - nr of modes: {}".format(Ur.shape[1]))
    logger.info("  - size of mode: {}".format(Ur.shape[0]))
    return Ur


#######################################
# Main
#######################################
parser = argparse.ArgumentParser(description="Computes energy and strain reduced bases.")
parser.add_argument('config_file', help="configuration file")
group = parser.add_mutually_exclusive_group()
group.add_argument('-e', '--only-energy', action="store_true", help="compute only energy reduced bases")
group.add_argument('-s', '--only-strain', action="store_true", help="compute only strain reduced bases")
args = parser.parse_args()
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
#handler = logging.FileHandler(config_filename.rsplit('.', 1)[0] + '.log')
#handler.setLevel(logging.INFO)
#logger.addHandler(handler)
conf = configparser.ConfigParser()
conf.read(args.config_file)
flag_comp_energy = True
flag_comp_strain =  True
if args.only_energy:
    flag_comp_strain =  False
elif args.only_strain:
    flag_comp_energy = False

if __name__ == '__main__':
    ene_bases_fname = conf['Parameters']['energy_bases_filename']
    str_bases_fname = conf['Parameters']['strain_bases_filename']
    bases_file_format = conf['Parameters']['bases_file_format']
    ene_e_files, ene_i_files, str_e_files, str_i_files = make_list_of_files(conf)

    if flag_comp_energy:
        t0 = time.time()
        logger.info("Generating energy bases")
        Ue = compute_elastic_modes(conf, ene_e_files, nr_components=1)
        Ui = compute_inelastic_modes(conf, ene_i_files, Ue, nr_components=1)
        U = np.hstack([Ue, Ui])
        t1 = time.time()
        write_bases(ene_bases_fname, U)
        logger.info("  SVD time: {:.1f}s".format(time.time() - t0))
        logger.info("  Writing time: {:.1f}s".format(time.time() - t1))
        logger.info("")

    if flag_comp_strain:
        t0 = time.time()
        logger.info("Generation strain bases")
        nr_strain_components = int(conf['Parameters']['nr_strain_components'])
        Ue = compute_elastic_modes(conf, str_e_files, nr_components=nr_strain_components)
        Ui = compute_inelastic_modes(conf, str_i_files, Ue, nr_components=nr_strain_components)
        U = np.hstack([Ue, Ui])
        t1 = time.time()
        write_bases(str_bases_fname, U)
        logger.info("  SVD time: {:.1f}s".format(time.time() - t0))
        logger.info("  Writing time: {:.1f}s".format(time.time() - t1))
        logger.info("")
