import json
import logging
import os

import numpy

logger = logging.getLogger(__name__)


def read_numpy_file(filename, file_format):
    if file_format == 'ascii':
        logger.debug("Reading numpy file in ASCII format")
        modes = numpy.loadtxt(filename)
    elif file_format == 'binary':
        logger.debug("Reading numpy file in BINARY format")
        modes = numpy.load(filename)
    elif file_format == 'auto':
        logger.debug("Autodetecting numpy file format")
        logger.debug("Autodetection currently not implemented")
        modes = numpy.empty([1, 1])
    return modes


def write_numpy_file(filename, file_format, U):
    if file_format == 'ascii':
        numpy.savetxt(filename, U)
    else:
        numpy.save(filename, U)
    return


def write_json(filename, data_dict):
    with open(filename, 'w') as fo:
        json.dump(data_dict, fo, indent=2)


def read_json(filename):
    with open(filename) as f:
        data_dict = json.load(f)
    return data_dict


def write_gid_header(filename):
    header = '''GiD Post Results File 1.0
GaussPoints "hex8_element_gp" ElemType Hexahedra
Number Of Gauss Points: 8
Natural Coordinates: Internal
End GaussPoints
'''
    fo = open(filename, 'w')
    fo.write(header)
    return fo


def write_gid_vector_field(fo, field_name, field, time):
    fo.write('Result "{}" "{}" {} Vector OnNodes\n'.format(
        field_name, "Kratos", float(time)))
    fo.write('Values\n')
    for i, e in enumerate(field):
        fo.write('{} {} {} {}\n'.format(i + 1, e[0], e[1], e[2]))
    fo.write('End Values\n')


def write_gid_matrix_field(fo, field_name, field, time):
    fo.write('Result "{}" "{}" {} Matrix OnGaussPoints "hex8_element_gp" \n'.format(
        field_name, "Kratos", float(time)))
    fo.write('Values\n')
    nr_ip = 8
    nr_elems = int(numpy.shape(field)[0] / nr_ip)
    row = 0
    for e in range(nr_elems):
        fo.write('{}'.format(e + 1))
        for ip in range(8):
            for c in field[row, :]:
                fo.write(' {}'.format(c))
            fo.write('\n')
            row = row + 1
    fo.write('End Values\n')


def read_gid_msh_nodes(filename):
    fo = open(filename, 'r')
    # parse coordinates
    nodes = []
    fo.readline()  # Header
    fo.readline()  # "Coordinates"
    for line in fo.readlines():
        if "End Coordinates" in line:
            break
        coordinates = line.strip().split()[1:]
        nodes.append(coordinates)
    return numpy.array(nodes).astype(numpy.double)


def write_strain_stress_header(filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    with open(filename, "w") as fo:
        fo.write("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
            "1", "2", "3", "4", "5", "6",  # strain
            "7", "8", "9", "10", "11", "12"))  # stress
        fo.write("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
            "strain XX", "YY", "ZZ", "XY", "YZ", "XZ",
            "stress XX", "YY", "ZZ", "XY", "YZ", "XZ"))


def write_strain_stress(filename, strain, stress):
    line = "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}\n".format(
        strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
        stress[0], stress[1], stress[2], stress[3], stress[4], stress[5])
    with open(filename, 'a') as ofile:
        ofile.write(line)


def write_strain_stress_ct_header(filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    with open(filename, "w") as fo:
        fo.write("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12}\n".format(
            "2", "3", "4", "5", "6", "7",  # strain
            "8", "9", "10", "11", "12", "13",  # stress
            "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
            "25", "26", "27", "28", "29", "30", "31", "32", "33", "34"))
        fo.write("{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "  # strain
                 "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "  # stress
                 "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                 "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
            "strain XX", "YY", "ZZ", "XY", "YZ", "XZ",
            "stress XX", "YY", "ZZ", "XY", "YZ", "XZ",
            "CT 11", "12", "13", "14", "15", "16", "22", "23", "24", "25", "26",
            "33", "34", "35", "36", "44", "45", "46", "55", "56", "66"))


def write_strain_stress_ct(filename, strain, stress, ct):
    line = "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}  {:<+1.4e}  " \
           "{:<+1.4e}\n".format(
        strain[0], strain[1], strain[2], strain[3], strain[4], strain[5],
        stress[0], stress[1], stress[2], stress[3], stress[4], stress[5],
        ct[0], ct[1], ct[2], ct[3], ct[4], ct[5],
        ct[7], ct[8], ct[9], ct[10], ct[11],
        ct[14], ct[15], ct[16], ct[17],
        ct[21], ct[22], ct[23],
        ct[28], ct[29],
        ct[35])
    with open(filename, 'a') as ofile:
        ofile.write(line)
