import configparser
import argparse
import numpy
import logging
import postprocess_utilities as util


def get_properties(rve_mdpa_filename, iw_list):
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
    out_prop = []
    for list in iw_list:
        e = int(list[0])
        out_prop.append(material[e])
    return out_prop


def get_elements_info(iw_list):
    out_e = []
    out_ip = []
    out_w = []
    out_gip = []
    for list in iw_list:
        out_e.append(int(list[0]))
        out_ip.append(int(list[1]))
        out_w.append(float(list[2]))
        out_gip.append(int(list[3]))
    return out_e, out_ip, out_w, out_gip


def parse_strain_bases(strain_bases_filename, iw_list, nr_modes):
    strain_bases = numpy.load(strain_bases_filename, mmap_mode='r')
    strain_bases = strain_bases[:,:nr_modes]
    nr_ip = 8
    nr_comps = 6
    out_B = []
    for list in iw_list:
        e = int(list[0])
        i = int(list[1])
        index = e * nr_ip * nr_comps + i * nr_comps
        B = strain_bases[index:index + nr_comps, :]
        out_B.append(B.tolist())
    return out_B


def create_rve_params_structure(strain_bases_filename, rve_mdpa_filename, nr_modes, reduced_ip_set):
    rve_params = {}
    # pack elements, ip and weights of reduced_ip
    out_e, out_lip, out_w, out_gip = get_elements_info(reduced_ip_set)
    rve_params['ip_element_id'] = out_e
    rve_params['ip_local_id'] = out_lip
    rve_params['ip_global_id'] = out_gip
    rve_params['ip_weight'] = out_w  # reduced ip weight
    # pack material (CL) index
    out_properties = get_properties(rve_mdpa_filename, reduced_ip_set)
    rve_params['ip_property_id'] = out_properties
    # read strain modes
    out_B = parse_strain_bases(strain_bases_filename, reduced_ip_set, nr_modes)
    rve_params['ip_strain_modes'] = out_B
    # add extra info
    rve_params['nr_modes'] = nr_modes
    rve_params['nr_reduced_ip'] = len(reduced_ip_set)
    return rve_params


#######################################
# Main
#######################################

# parse command line arguments
#parser = argparse.ArgumentParser(description="Create reduced RVE data set")
#parser.add_argument('config_file', help="configuration file")
#args = parser.parse_args()
# parse configuration file
#conf = configparser.ConfigParser()
#conf.read(args.config_file)
# configure logger
#verbosity_level = logging.INFO
#logging.basicConfig(format='[%(asctime)s] %(message)s',
#                    datefmt='%H:%M:%S',level=verbosity_level)
#logger = logging.getLogger(__name__)

if __name__ == '__main__':
    strain_bases_filename = "bases_strain.npy"
    rve_mdpa_filename = "../model.mdpa"
    reduced_ip_set = numpy.loadtxt("../roq_list.dat")
    nr_reduced_ip = 50
    nr_modes = 15

    rve_params = create_rve_params_structure(
        strain_bases_filename, rve_mdpa_filename, nr_modes, reduced_ip_set)
    filename = "rve_{}m_{}ip.json".format(nr_modes, nr_reduced_ip)
    util.write_json(filename, rve_params)
