import sys
import configparser
import argparse
import glob
import numpy as np
import logging

def make_list_of_files(conf):

    trajectory_filename = conf['Parameters']['trajectory_filename']
    nr_elastic_snapshots_filename = conf['Parameters']['nr_elastic_snapshots_filename']
    energy_filename = conf['Parameters']['energy_filename']
    strain_filename = conf['Parameters']['strain_filename']

    trajectory_paths = sorted(glob.glob("{}_?".format(trajectory_filename)))
    energy_elastic_files = []
    energy_inelast_files = []
    strain_elastic_files = []
    strain_inelast_files = []
    for path in trajectory_paths:
        with open("{}/{}".format(path, nr_elastic_snapshots_filename), "r") as f:
            nr_elastic_snapshots = int(f.readline().strip())
        logger.info("Nr of elastic snapshots: {}".format(nr_elastic_snapshots))
        energy_elastic_files.extend(sorted(glob.glob("{}/{}*".format(path, energy_filename)))[:nr_elastic_snapshots])
        energy_inelast_files.extend(sorted(glob.glob("{}/{}*".format(path, energy_filename)))[nr_elastic_snapshots:])
        strain_elastic_files.extend(sorted(glob.glob("{}/{}*".format(path, strain_filename)))[:nr_elastic_snapshots])
        strain_inelast_files.extend(sorted(glob.glob("{}/{}*".format(path, strain_filename)))[nr_elastic_snapshots:])

    return energy_elastic_files, energy_inelast_files, strain_elastic_files, strain_inelast_files


def compute_strain_elastic(conf, strain_elastic_files):
    logger.info("SVD of strain elastic snapshots")

    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_strain_components = int(conf['Parameters']['nr_strain_components'])
    tol_svd_elastic_strain = float(conf['Parameters']['tolerance_svd_elastic_strain'])

    nr_dofs = nr_elements * nr_integration_points * nr_strain_components
    X = np.empty([nr_dofs, len(strain_elastic_files)])
    for i, file in enumerate(strain_elastic_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
    # TODO: incluir el SVD del paquete scipy, parece que es mas optimo, ver si esta instalada una version reciente de scipy en el cluster
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    logger.info(S)

    logger.info("selection of strain elastic modes")
    cont = 1
    Us_el = []
    for iValue, singular_value in enumerate(S):
        if singular_value > tol_svd_elastic_strain:
            if cont == 1:
                Us_el = U[:, iValue]
            else:
                Us_el = np.column_stack((Us_el, U[:, iValue]))
            cont = cont + 1
    logger.info(Us_el.shape)

    return Us_el


def compute_strain_inelastic(conf, strain_inelast_files, Us_el):
    logger.info("projection of strain inelastic snapshots")

    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_strain_components = int(conf['Parameters']['nr_strain_components'])
    tol_svd_inelastic_strain = float(conf['Parameters']['tolerance_svd_inelastic_strain'])

    nr_dofs = nr_elements * nr_integration_points * nr_strain_components

    X = np.empty([nr_dofs, len(strain_inelast_files)])
    for i, file in enumerate(strain_inelast_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        for j in range(Us_el.shape[1]):
            X[:,i] = X[:,i] - np.multiply(np.dot(Us_el[:,j],X[:,i]),Us_el[:,j])

    logger.info("SVD of strain inelastic modified snapshots")
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    logger.info(S)

    logger.info("selection of strain inelastic modes")
    cont = 1
    Us_in = []
    for iValue, singular_value in enumerate(S):
        if singular_value > tol_svd_inelastic_strain:
            if cont == 1:
                Us_in = U[:, iValue]
            else:
                Us_in = np.column_stack((Us_in, U[:, iValue]))
            cont = cont + 1
    logger.info(Us_in.shape)

    return Us_in


def compute_energy_elastic(conf, energy_elastic_files):

    logger.info("SVD of elastic energy snapshots")

    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_strain_components = int(conf['Parameters']['nr_strain_components'])
    nr_modes = 2 * nr_strain_components

    nr_dofs = nr_elements * nr_integration_points
    X = np.empty([nr_dofs, len(energy_elastic_files)])
    for i, file in enumerate(energy_elastic_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]

    Ue_el = U[:,:nr_modes]
    logger.info(S[:nr_modes])
    logger.info(Ue_el.shape)

    return Ue_el


def compute_energy_inelastic(conf, energy_inelast_files, Ue_el):

    logger.info("projection of energy inelastic snapshots")

    nr_elements = int(conf['Parameters']['nr_elements'])
    nr_integration_points = int(conf['Parameters']['nr_integration_points'])
    nr_modes = int(conf['Parameters']['max_nr_reduced_modes'])

    nr_dofs = nr_elements * nr_integration_points

    X = np.empty([nr_dofs, len(energy_inelast_files)])
    for i, file in enumerate(energy_inelast_files):
        X[:, i] = np.fromfile(file, dtype=np.float32)
        for j in range(Ue_el.shape[1]):
            X[:,i] = X[:,i] - np.multiply(np.dot(Ue_el[:,j],X[:,i]),Ue_el[:,j])

    logger.info("SVD of inelastic energy modified snapshots")
    [U, S] = np.linalg.svd(X, full_matrices=False)[:2]
    Ue_in = U[:,:nr_modes]

    logger.info(Ue_in.shape)

    return Ue_in


#######################################
# Main
#######################################
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

parser = argparse.ArgumentParser(description="Computes energy and strain reduced bases.")
parser.add_argument('config_file', help="configuration file")
group = parser.add_mutually_exclusive_group()
group.add_argument('-e', '--only-energy', action="store_true", help="compute only energy reduced bases")
group.add_argument('-s', '--only-strain', action="store_true", help="compute only strain reduced bases")
args = parser.parse_args()

flag_comp_energy = True
flag_comp_strain =  True
if args.only_energy:
    flag_comp_strain =  False
elif args.only_strain:
    flag_comp_energy = False

if __name__ == '__main__':
    conf = configparser.ConfigParser()
    conf.read(args.config_file)
    energy_bases_filename = conf['Parameters']['energy_bases_filename']
    strain_bases_filename = conf['Parameters']['strain_bases_filename']
    bases_file_format = conf['Parameters']['bases_file_format']

    energy_elastic_files, energy_inelast_files, strain_elastic_files, strain_inelast_files = make_list_of_files(conf)

    if flag_comp_energy:
        logger.info("Computing energy snapshots")
        Ue_el = compute_energy_elastic(conf, energy_elastic_files)
        Ue_in = compute_energy_inelastic(conf, energy_inelast_files, Ue_el)
        Ue = np.hstack([Ue_el, Ue_in])
        if bases_file_format == 'ascii':
            with open(energy_bases_filename, 'wb') as ofile:
                np.savetxt(ofile, Ue, fmt='%.17f')
        else:
            logger.info("binary output not implemented yet. writing in ascii")
            with open(energy_bases_filename, 'wb') as ofile:
                np.savetxt(ofile, Ue, fmt='%.17f')
    if flag_comp_strain:
        logger.info("Computing strain bases")
        Us_el = compute_strain_elastic(conf, strain_elastic_files)
        Us_in = compute_strain_inelastic(conf, strain_inelast_files, Us_el)
        Us = np.hstack([Us_el, Us_in])
        if bases_file_format == 'ascii':
            with open(strain_bases_filename, 'wb') as ofile:
                np.savetxt(ofile, Us, fmt='%.17f')
        else:
            logger.info("binary output not implemented yet. writing in ascii")
            with open(strain_bases_filename, 'wb') as ofile:
                np.savetxt(ofile, Us, fmt='%.17f')

