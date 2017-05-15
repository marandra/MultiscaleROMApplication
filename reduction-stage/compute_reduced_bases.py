import sys
import configparser
import glob
import numpy as np

'''
def read_config(fname):
    conf = configparser.ConfigParser()
    file_format = config['Parameters']['file_format']
    nr_elastic_snapshots = config['Parameters']['nr_elastic_snapshots']
    nr_elements = config['Parameters']['nr_elements']
    nr_integration_points = config['Parameters']['nr_integration_points']
    return file_format, nr_elastic_snapshots


def get_list_of_files(fname):
    glog.glob['fname*']


def read_file_to_numpy(fname, nfields, dtype=np.float32):
    val = np.fromfile(fname, dtype=dtype)
    return val
'''

# get parameters
fname = sys.argv[1]
conf = configparser.ConfigParser()
conf.read(fname)
file_format = conf['Parameters']['file_format']
nr_elastic_snapshots = int(conf['Parameters']['nr_elastic_snapshots'])
nr_elements = int(conf['Parameters']['nr_elements'])
nr_integration_points = int(conf['Parameters']['nr_integration_points'])
nr_strain_components = int(conf['Parameters']['nr_strain_components'])
energy_file_name = conf['Parameters']['energy_file_name']
strain_file_name = conf['Parameters']['strain_file_name']
energy_elastic_files = sorted(glob.glob(energy_file_name + '*'))[:nr_elastic_snapshots]
strain_elastic_files = sorted(glob.glob(strain_file_name + '*'))[:nr_elastic_snapshots]
energy_inelast_files = sorted(glob.glob(energy_file_name + '*'))[nr_elastic_snapshots:]
strain_inelast_files = sorted(glob.glob(strain_file_name + '*'))[nr_elastic_snapshots:]
print(file_format)
print(nr_elastic_snapshots)
print(nr_elements)
print(nr_integration_points)
print(energy_elastic_files)
print(strain_elastic_files)
print(energy_inelast_files)
print(strain_inelast_files)

# first part: read and compute elastic energy modes, compute projector
nr_dofs = nr_elements * nr_integration_points
X = np.empty([nr_dofs, len(energy_elastic_files)])
for i, file in enumerate(energy_elastic_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
Ue_el = np.linalg.svd(X, full_matrices=False)[0]
P = np.identity(nr_dofs, dtype=np.float32) - np.dot(Ue_el, Ue_el.T)

# second part: read inelastic energy modes, remove elastic component, decomp svd
X = np.empty([nr_dofs, len(energy_inelast_files)])
for i, file in enumerate(energy_inelast_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
X = np.dot(P, X)
Ue_in = np.linalg.svd(X, full_matrices=False)[0]

# third part: read and compute elastic strain modes
nr_dofs = nr_elements * nr_integration_points * nr_strain_components
X = np.empty([nr_dofs, len(strain_elastic_files)])
for i, file in enumerate(strain_elastic_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
Us_el = np.linalg.svd(X, full_matrices=False)[0]
P = np.identity(nr_dofs, dtype=np.float32) - np.dot(Us_el, Us_el.T)

# forth part: read inelastic strain modes, remove elastic component, decomp svd
X = np.empty([nr_dofs, len(strain_inelast_files)])
for i, file in enumerate(strain_inelast_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
X = np.dot(P, X)
Us_in = np.linalg.svd(X, full_matrices=False)[0]

