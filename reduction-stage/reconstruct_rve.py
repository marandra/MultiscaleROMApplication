import argparse
import logging
import numpy
#import meshio


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


#######################################
# Main
#######################################

# parse command line arguments
parser = argparse.ArgumentParser(description="Reconstructs RVE displacement"
                                             " field from strain mode weights")
parser.add_argument('gid_msh_file', help="gid output .msh file")
parser.add_argument('mode_node_file', help="modes_displacement file")
#parser.add_argument('macro_strain_file', help="macro strain file per timestep")
parser.add_argument('weights_file', help="mode weights file."
                         "Each row is a timestep, columns are mode weights")
parser.add_argument('-v', '--verbose', action="store_true", help="shows debug information")
args = parser.parse_args()

# configure logger
verbosity_level = logging.INFO
if args.verbose:
    verbosity_level = logging.DEBUG
logging.basicConfig(format='[%(asctime)s] %(message)s',
                    datefmt='%H:%M:%S',level=verbosity_level)
logger = logging.getLogger(__name__)

if __name__ == '__main__':
    logger.info("Reconstruct RVE fluctuant displacement field")

    logger.info("Loading modes_displacement files")
    mode_node_matrix = numpy.loadtxt(args.mode_node_file)
    #mode_node_matrix = numpy.load(args.mode_node_file)

    logger.info("Loading modes_weights files")
    #macro_strain = numpy.loadtxt(args.weights_file)[:, 0]
    print(args.weights_file)
    weights = numpy.loadtxt(args.weights_file)
    nr_timesteps = numpy.shape(weights)[0]
    nr_modes = numpy.shape(weights)[1]
    logger.debug("Number of timesteps detected: {}".format(nr_timesteps))
    logger.debug("Number of modes detected: {}".format(nr_modes))
    #logger.debug("Macro strain: ")
    #logger.debug(macro_strain)
    logger.debug("Mode weights: ")
    logger.debug(weights)

    logger.info("Solving fluctuant displacement")
    filename = args.gid_msh_file.rsplit(".", 1)[0] + ".res"
    f = write_gid_header(filename)
    for t in range(nr_timesteps):
        logger.info("Timestep {}".format(t))
        displacement = numpy.dot(mode_node_matrix, weights[t, :])
        #print(displacement)
        print(mode_node_matrix)
        displacement_form = numpy.reshape(displacement, (-1, 3))
        nnode = displacement_form.shape[0]
        #print(nnode)
        #gid_output = numpy.hstack([numpy.arange(1, nnode+1).reshape(-1,1).astype(int), displacement_form])
        #numpy.savetxt(filename, gid_output)

        write_gid_vector_field(f, displacement_form, t)
