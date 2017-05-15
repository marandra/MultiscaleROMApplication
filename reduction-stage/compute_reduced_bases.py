import sys
import configparser
import glob
import numpy as np

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

trajectory_paths = glob.glob('*_?')
energy_elastic_files = []
strain_elastic_files = []
energy_inelast_files = []
strain_inelast_files = []
for path in trajectory_paths:
    energy_elastic_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[:nr_elastic_snapshots])
    strain_elastic_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[:nr_elastic_snapshots])
    energy_inelast_files.extend(sorted(glob.glob(path + '/' + energy_file_name + '*'))[nr_elastic_snapshots:])
    strain_inelast_files.extend(sorted(glob.glob(path + '/' + strain_file_name + '*'))[nr_elastic_snapshots:])
    
#import pprint
#print("AAAAAAA")
#pprint.pprint(energy_elastic_files)
#print("BBBBBB")
#pprint.pprint(strain_elastic_files)
#print("CCCCCCC")
#pprint.pprint(energy_inelast_files)
#print("DDDDDDD")
#pprint.pprint(strain_inelast_files)


# first part: read and compute elastic energy modes, compute projector
nr_dofs = nr_elements * nr_integration_points
X = np.empty([nr_dofs, len(energy_elastic_files)])
for i, file in enumerate(energy_elastic_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
Ue_el = np.linalg.svd(X, full_matrices=False)[0]
print(Ue_el)
print("Fin svd strain elastic")
print(X.shape)

sys.exit()

# second part: read inelastic energy modes, remove elastic component, decomp svd
X = np.empty([nr_dofs, len(energy_inelast_files)])
for i, file in enumerate(energy_inelast_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
print("Antes de proyector")
X = X - np.dot(Ue_el, np.dot(Ue_el.T, X))
print("Desp de proyector")
Ue_in = np.linalg.svd(X, full_matrices=False)[0]
print("Fin svd energy inelastic")

# third part: read and compute elastic strain modes
nr_dofs = nr_elements * nr_integration_points * nr_strain_components
X = np.empty([nr_dofs, len(strain_elastic_files)])
for i, file in enumerate(strain_elastic_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
Us_el = np.linalg.svd(X, full_matrices=False)[0]
print("Fin svd strain elastic")

# forth part: read inelastic strain modes, remove elastic component, decomp svd
X = np.empty([nr_dofs, len(strain_inelast_files)])
for i, file in enumerate(strain_inelast_files):
    X[:, i] = np.fromfile(file, dtype=np.float32)
print("Antes de proyector")
X = X - np.dot(Us_el, np.dot(Us_el.T, X))
print("Desp de proyector")
Us_in = np.linalg.svd(X, full_matrices=False)[0]
print("fin")
