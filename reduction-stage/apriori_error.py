import time
import configparser
import argparse
import glob
import numpy as np
import scipy.sparse.linalg as sp
import logging
import pprint as pp


def make_list_of_files(conf, traj):
    logger.debug("Scanning trajectories")
    trajectory_filename = conf['Parameters']['trajectory_filename']
    nr_e_snap_filename = conf['Parameters']['nr_elastic_snapshots_filename']
    energy_filename = conf['Parameters']['energy_filename']
    strain_filename = conf['Parameters']['strain_filename']
    trajectory_paths = sorted(glob.glob("{}_{}".format(trajectory_filename, traj)))
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


def load_snapshots(conf, files, nr_components):
    logger.info("Loading snapshots")
    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_dofs = nr_elements * nr_integration_points * nr_components
    X = np.empty([nr_dofs, len(files)])
    total = len(files)
    batch_size = int(total / 10 + .5)
    counter = 1
    for i, file in enumerate(files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        if not counter % batch_size:
            logger.debug("    {}/{} snapshots processed".format(counter, total))
        counter = counter + 1
    logger.info("")
    return X

def compute_error(X, Ue):
    aux = np.dot(Ue.T, X)
    err = X - np.dot(Ue, aux)
    nerr = np.linalg.norm(err)
    normX = np.linalg.norm(X)
    return nerr / normX


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Computes energy and strain reduced bases.")
parser.add_argument('config_file', help="configuration file")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
parser.add_argument('-n', '--nmodes', help="performs iterative svd (svds algorithm)")
parser.add_argument('-t', help="performs iterative svd (svds algorithm)")
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

nr_modes = int(args.nmodes)
traj = args.t

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
    ene_e_files, ene_i_files, str_e_files, str_i_files = make_list_of_files(conf, traj)

    logger.debug('Config parameters')
    # TODO: not working properly
    logger.debug(pp.pprint(conf.items('Parameters')))
    logger.debug('')
    fo = open("apriori_error.out", 'w')
    fo.write("#{}   {}\n".format("modes", "norm_err"))

    if flag_comp_energy:
        logger.info("apriori error ENERGY")
        nr_components = 1
        files = ene_e_files + ene_i_files
        bases_fname = ene_bases_fname

    elif flag_comp_strain:
        logger.info("apriori error STRAIN")
        nr_components = int(conf['Parameters']['nr_strain_components'])
        files = str_e_files + str_i_files
        bases_fname = str_bases_fname
    else:
        print("error: select -s or -e")
        exit()

    Ut = np.load(bases_fname, mmap_mode='r') 
    snapshots = load_snapshots(conf, files, nr_components)
    for n in range(nr_modes):
        U = Ut[:, :n] 
        nerr = compute_error(snapshots, U)
        logger.info("modes: {} norm error: {}".format(n, nerr))
        fo.write("{}   {}\n".format(n, nerr))
    fo.close()

