import numpy
import logging
import json

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
        modes = numpy.empty([1,1])
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
    nr_elems = numpy.shape(field)[0] / nr_ip
    row = 0
    for e in range(nr_elems):
        fo.write('{}'.format(e + 1))
        for ip in range(8):
            for c in field[row, :]:
                fo.write(' {}'.format(c))
                row = row + 1
        fo.write('\n')
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
