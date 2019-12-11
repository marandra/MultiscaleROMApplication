import numpy
import KratosMultiphysics.MultiscaleROMApplication.io_utilities as io_utilities


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
    for l in iw_list:
        e = int(l[0])
        out_prop.append(material[e])
    return out_prop


def unpack_ip_data(iw_list):
    out_e = []
    out_ip = []
    out_w = []
    out_gip = []
    for l in iw_list:
        out_e.append(int(l[0]))
        out_ip.append(int(l[1]))
        out_w.append(float(l[2]))
        out_gip.append(int(l[3]))
    return out_e, out_ip, out_w, out_gip


def parse_strain_bases(strain_bases_filename, iw_list, nr_modes):
    strain_bases = numpy.load(strain_bases_filename, mmap_mode='r')
    strain_bases = strain_bases[:, :nr_modes]
    nr_comps = 6
    out_B = []
    for l in iw_list:
        gip = int(l[3])
        index = gip * nr_comps
        B = strain_bases[index:index + nr_comps, :]
        out_B.append(B.tolist())
    return out_B


def create_rve_params_structure(strain_bases_filename, rve_mdpa_filename,
                                rve_materials_filename, nr_modes, reduced_ip_set):
    rve_params = {}
    # pack elements, ip and weights of reduced_ip
    out_e, out_lip, out_w, out_gip = unpack_ip_data(reduced_ip_set)
    rve_params['ip_element_id'] = out_e
    rve_params['ip_local_id'] = out_lip
    rve_params['ip_global_id'] = out_gip
    rve_params['ip_weight'] = out_w  # reduced ip weight
    # RVE material properties
    rve_params['material_parameters'] = io_utilities.read_json(rve_materials_filename)
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
