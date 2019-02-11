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
