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


def write_gid_vector_field(fo, field, time):
    fo.write('Result "{}" "{}" {} Vector OnNodes\n'.format("DISPLACEMENT", "Kratos", float(time), 'Vector OnNodes'))
    fo.write('Values\n')
    for i, e in enumerate(field):
        fo.write('{} {} {} {}\n'.format(i + 1, e[0], e[1], e[2]))
    fo.write('End Values\n')

