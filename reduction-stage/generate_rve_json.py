import configparser
import numpy as np
import sys
import json


#################################33
if __name__ == "__main__":

    conf = configparser.ConfigParser()
    conf.read("reduced_bases.cfg")
    strain_bases_filename = conf['Parameters']['strain_bases_filename']
    nr_ip = int(conf['Parameters']['nr_integration_points'])
    voigtsize = int(conf['Parameters']['nr_strain_components'])

    nr_dofs = nr_ip * voigtsize
    nr_modes = 4

    strain_bases = np.load(strain_bases_filename, mmap_mode='r')
    strain_bases = strain_bases[:,:nr_modes]

    # read roq list
    iw_list = np.loadtxt("roq_list.dat")

    # read model materials
    material = []
    flag_elements = False
    with open("model.mdpa", "r") as fi:
        for line in fi.readlines():
            if not flag_elements:
                if "Begin Elements" not in line:
                    continue
                else:
                    flag_elements = True
            else:
                if "End Elements" in line:
                    break
                else:
                    material.append(int(line.split()[1]))

    out = {}
    out_B = []
    out_w = []
    out_prop = []
    B = np.empty((voigtsize, nr_modes))
    for list in iw_list:
        e = int(list[0])
        i = int(list[1])
        w = float(list[2])

        # get B
        index = e * nr_ip * voigtsize + i * voigtsize
        B = strain_bases[index:index + voigtsize, :]

        out_B.append(B.tolist())
        out_w.append(w)
        out_prop.append(material[e])
        

    out['props_id'] = out_prop
    out['w'] = out_w
    out['B'] = out_B

    with open("rve.json", 'w') as fo:
        json.dump(out, fo, indent=2)
