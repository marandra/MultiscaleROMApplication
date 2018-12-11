import configparser
import argparse
import numpy as np
import logging
import json


def generate_rve_params(conf, iw_list):
    strain_bases_filename = conf['Parameters']['strain_bases_filename']
    nr_ip = int(conf['Parameters']['nr_integration_points'])
    nr_comps = int(conf['Parameters']['nr_strain_components'])
    nr_modes = int(conf['Parameters']['nr_active_modes'])
    rve_mdpa_filename = conf['Parameters']['rve_mdpa_filename']
    nr_dofs = nr_ip * nr_comps
    strain_bases = np.load(strain_bases_filename, mmap_mode='r')
    strain_bases = strain_bases[:,:nr_modes]

    # read model materials
    material = {}
    flag_elements = False
    with open(rve_mdpa_filename, 'r') as fi:
        for line in fi.readlines():
            if not flag_elements:
                if "Begin Elements" not in line:
                    continue
                else:
                    flag_elements = True
            else:
                if "End Elements" in line:
                    flag_elements = False
                    continue
                else:
                    material[int(line.split()[0]) - 1] = int(line.split()[1])
    out = {}
    out_B = []
    out_w = []
    out_prop = []
    B = np.empty((nr_comps, nr_modes))
    for list in iw_list:
        e = int(list[0])
        i = int(list[1])
        w = float(list[2])

        # get B
        index = e * nr_ip * nr_comps + i * nr_comps
        B = strain_bases[index:index + nr_comps, :]

        out_B.append(B.tolist())
        out_w.append(w)
        out_prop.append(material[e])

    out['props_id'] = out_prop
    out['w'] = out_w
    out['B'] = out_B

    return out

#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Create reduced RVE data set")
parser.add_argument('config_file', help="configuration file")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
args = parser.parse_args()

# parse configuration file
conf = configparser.ConfigParser()
conf.read(args.config_file)

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s',
                    datefmt='%H:%M:%S',level=verbosity_level)
logger = logging.getLogger(__name__)
handler = logging.FileHandler('log_' + args.config_file.rsplit('.', 1)[0])
handler.setLevel(logging.DEBUG)
logger.addHandler(handler)

if __name__ == '__main__':
    logger.info("Generating RVE parameters for CL")
    np.loadtxt("roq_list.dat", roq_list)
    rve_params = generate_rve_params(conf, roq_list)
    logging.debug("ROQ list size {}".format(np.shape(roq_list)))
    #filename = conf['Parameters']['roq_weights_filename']
    #np.savetxt(filename, roq_mask)
    with open("rve.json", 'w') as fo:
        json.dump(rve_params, fo, indent=2)
