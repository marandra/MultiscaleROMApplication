import sys
import struct


def create_header(filename, nfields, field_format):
    with open(filename + '.hdr', 'w+') as fo:
        fo.write('snapshot_fields {}\n'.format(nfields))
        fo.write('field_format {}\n'.format(field_format))
        fo.write('elastic_modes {}\n'.format(0))


def create_snapshots_file_ascii(filename):
    return open(filename + '.out', 'w+')


def create_snapshots_file_binary(filename):
    return open(filename + '.bin', 'wb+')


def append_to_file_ascii(fo, values):
    for v in values:
        fo.write("{} ".format(v))
    fo.write("\n")


def append_to_file_binary(fo, values, field_format):
    if field_format == 'float16':
        fformat = 'e'
    if field_format == 'float32':
        fformat = 'f'
    for v in values:
        fo.write(struct.pack(fformat, v))


def close_snapshots_file_ascii(fo):
   fo.write('\n')
   fo.close()

def close_snapshots_file_binary(fo):
   fo.write(b'\n')
   fo.close()


def parse_msh(fname):
    with open(fname + '.msh', 'r') as fi:
        fi.readline()
        count_nodes_flag = False
        nnode = 0
        elem_mat = {}
        for line_raw in fi.readlines():
            line = line_raw.strip()
            if line == 'Coordinates':
                count_nodes_flag = True
                continue
            if line == 'End Coordinates':
                count_nodes_flag = False
                continue
            if count_nodes_flag:
                nnode = nnode + 1
                continue
            if 'Elem' in line:
                continue
            elem = int(line.split()[0])
            material = int(line.split()[-1])
            elem_mat[elem] = material
    return nnode


def parse_res_and_write(fname, nnode, fo, var_name='DISPLACEMENT',
                        field_format='float32', mode='binary'):
    with open(fname + '.res', 'r') as fi:
        fi.readline()
        fi.readline()
        ngauss = int(fi.readline().strip().split(':')[-1])
        fi.readline()
        fi.readline()
        fi.readline()
        fi.readline()
        fi.readline()
        fi.readline()
        for line in iter(fi):
            if var_name in line:
                fi.readline()
                values = []
                for i in range(nnode):
                    values.extend([float(x) for x in fi.readline().strip().split()[1:]])
                if mode == 'ascii':
                    append_to_file_ascii(fo, values)
                else:
                    append_to_file_binary(fo, values, field_format)
                fi.readline()
    return len(values), field_format


########
fname = sys.argv[1]
fn_snapshots_mode = 'binary'
fn_snapshots = 'snapshots'
nnode = parse_msh(fname)
if fn_snapshots_mode == 'ascii':
    fo = create_snapshots_file_ascii(fn_snapshots)
    nfield, field_format = parse_res_and_write(fname, nnode, fo, mode='ascii')
    close_snapshots_file_ascii(fo)
else:
    fo = create_snapshots_file_binary(fn_snapshots)
    nfield, field_format = parse_res_and_write(fname, nnode, fo, mode='binary')
    close_snapshots_file_binary(fo)
create_header(fn_snapshots, nfield, field_format)
