import time
import configparser
import argparse
import glob
import numpy as np
import scipy.sparse.linalg as sp
import logging
import pprint as pp


def write_bases(filename, U):
    bases_file_format = conf['Parameters']['bases_file_format']
    if bases_file_format == 'ascii':
        np.savetxt(filename, U)
    else:
        np.save(filename, U)
    return


def _make_list_of_files(conf):
    logger.debug("Scanning trajectories")
    trajectory_filename = conf['Parameters']['trajectory_filename']
    nr_e_snap_filename = conf['Parameters']['nr_elastic_snapshots_filename']
    energy_filename = conf['Parameters']['energy_filename']
    strain_filename = conf['Parameters']['strain_filename']
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))
    ene_files = []
    str_files = []
    for path in trajectory_paths:
        ene_files.extend(sorted(glob.glob("{}/{}*".format(path, energy_filename))))
        str_files.extend(sorted(glob.glob("{}/{}*".format(path, strain_filename))))
    return ene_files, str_files


def split_files_elastic_inelastic(conf, ene_files, str_files):
    logger.debug("Scanning trajectories")
    trajectory_filename = conf['Parameters']['trajectory_filename']
    nr_e_snap_filename = conf['Parameters']['nr_elastic_snapshots_filename']
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))
    ene_e_files = []
    ene_i_files = []
    str_e_files = []
    str_i_files = []
    nr_elastic_snap = []
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


def make_list_of_files(conf):
    logger.debug("Scanning trajectories")
    trajectory_filename = conf['Parameters']['trajectory_filename']
    nr_e_snap_filename = conf['Parameters']['nr_elastic_snapshots_filename']
    energy_filename = conf['Parameters']['energy_filename']
    strain_filename = conf['Parameters']['strain_filename']
    trajectory_paths = sorted(glob.glob("{}_*".format(trajectory_filename)))
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


def compute_modes(conf, files, nr_modes, nr_components, Ue=None):
    logger.info("Loading snapshots")
    if Ue is not None:
        logger.info("and removing elastic component")
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
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
    if args.iterative:
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


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Computes energy and strain reduced bases.")
parser.add_argument('config_file', help="configuration file")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
parser.add_argument('-i', '--iterative', action="store_true", help="performs iterative svd (svds algorithm)")
group = parser.add_mutually_exclusive_group()
group.add_argument('-e', '--only-energy', action="store_true", help="compute only energy reduced bases")
group.add_argument('-s', '--only-strain', action="store_true", help="compute only strain reduced bases")
args = parser.parse_args()

# parse configuration file
conf = configparser.ConfigParser()
conf.read(args.config_file)
flag_comp_energy = True
flag_comp_strain =  True
if args.only_energy:
    flag_comp_strain =  False
elif args.only_strain:
    flag_comp_energy = False

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s',
                    datefmt='%H:%M:%S', level=verbosity_level)
logger = logging.getLogger(__name__)
handler = logging.FileHandler('log_' + args.config_file.rsplit('.', 1)[0], mode='w')
handler.setLevel(logging.DEBUG)
handler.setFormatter(logging.Formatter('[%(asctime)s] %(message)s'))
logger.addHandler(handler)

if __name__ == '__main__':
    ene_bases_fname = conf['Parameters']['energy_bases_filename']
    str_bases_fname = conf['Parameters']['strain_bases_filename']
    bases_file_format = conf['Parameters']['bases_file_format']
    ene_e_files, ene_i_files, str_e_files, str_i_files = make_list_of_files(conf)

    ##filtering out snapshots
    #min_nr_elast_snp = 13
    #ind1 = list(np.linspace(0, min_nr_elast_snp, num=3, dtype=np.dtype('uint8')))
    #import pprint
    #print(ind1)
    #pprint.pprint(ene_e_files)
    #ene_e_files = [ene_e_files[i] for i in ind1]
    #pprint.pprint(ene_e_files)
    #err
    ##ene_i_files =
    ##str_e_files =
    ##str_i_files =

    logger.debug('Config parameters"')
    # TODO: not working properly
    logger.debug(pp.pprint(conf.items('Parameters')))
    logger.debug('')

    if flag_comp_energy:
        t0 = time.time()
        logger.info("Generating bases ENERGY")
        nr_elastic_modes = int(conf['Parameters']['nr_elastic_modes_energy'])
        nr_inelastic_modes = int(conf['Parameters']['nr_inelastic_modes_energy'])
        if nr_elastic_modes > 0:
            logger.info("Processing elastic snapshots")
            Ue = compute_modes(conf, ene_e_files, nr_elastic_modes, nr_components=1)
            logger.info("Processing inelastic snapshots")
            Ui = compute_modes(conf, ene_i_files, nr_inelastic_modes, nr_components=1, Ue=Ue)
            U = np.hstack([Ue, Ui])
        else:
            logger.info("Nr of elastic modes set to zero -> Not discriminating elastic/inelastic snapshots")
            U = compute_modes(conf, ene_e_files + ene_i_files, nr_inelastic_modes, nr_components=1)
        t1 = time.time()
        write_bases(ene_bases_fname, U)
        logger.info("  SVD time: {:.1f}s".format(time.time() - t0))
        logger.info("  Writing time: {:.1f}s".format(time.time() - t1))
        logger.info("")

    if flag_comp_strain:
        t0 = time.time()
        logger.info("Generation bases STRAIN")
        nr_strain_components = int(conf['Parameters']['nr_strain_components'])
        nr_inelastic_modes = int(conf['Parameters']['nr_inelastic_modes_strain'])
        nr_elastic_modes = int(conf['Parameters']['nr_elastic_modes_strain'])
        # discriminate elastic snapshots if nr_elastic modes is True
        if nr_elastic_modes > 0:
            logger.info("Processing elastic snapshots")
            Ue = compute_modes(conf, str_e_files, nr_elastic_modes, nr_components=nr_strain_components)
            logger.info("Processing inelastic snapshots")
            Ui = compute_modes(conf, str_i_files, nr_inelastic_modes, nr_components=nr_strain_components, Ue=Ue)
            U = np.hstack([Ue, Ui])
        else:
            logger.info("Nr of elastic modes set to zero -> Not discriminating elastic/inelastic snapshots")
            U = compute_modes(conf, str_e_files + str_i_files, nr_inelastic_modes, nr_components=nr_strain_components)
        t1 = time.time()
        write_bases(str_bases_fname, U)
        logger.info("  SVD time: {:.1f}s".format(time.time() - t0))
        logger.info("  Writing time: {:.1f}s".format(time.time() - t1))
        logger.info("")
