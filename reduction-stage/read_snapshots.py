import sys
import struct
import numpy as np


def read_header_file(fname):
    with open(fname, 'r') as fi:
        nfields = int(fi.readline().strip().split()[-1])
        ffmt = fi.readline().strip().split()[-1]
    return nfields, ffmt


def read_file_np(fname, nfields, dtype=np.float32):
    # TODO change nfield to nnode and ndofs
    val = np.fromfile(fname, dtype=dtype)
    val = val.reshape((-1, nfields)).T
    return val


def create_lumped_mass_matrix(fname):
    # TODO change nfield to nnode and ndofs
    # TODO save connectivity, read it here
    nnode = 9
    m_lump = [0] * nnode * 3  # 3 is dofs
    # change this from here
    connect = {}
    connect[1] = [2, 5, 9, 6]
    connect[2] = [6, 9, 7, 3]
    connect[3] = [5, 1, 8, 9]
    connect[4] = [9, 8, 4, 7]
    # to here

    with open(fname, 'r') as fi:
        for line in fi.readlines():
            elem = int(line.strip().split()[0])
            nlocnods = int(line.strip().split()[-1].split(',')[0].split('[')[1])
            mat = line.split('((')[1].split('))')[0].split('),(')
            ldiag = []
            for i in range(0, nlocnods, 2):
                ldiag.append(float(mat[i].split(',')[i]))
            for e_l, e_g in enumerate(connect[elem]):
                # make a for loop over dofs:
                # for i in range(ndof):
                # m_lump[(e_g - 1) * ndof + i] += ldiag[e_l]
                m_lump[(e_g - 1) * 3 + 0] += ldiag[e_l]
                m_lump[(e_g - 1) * 3 + 1] += ldiag[e_l]
                m_lump[(e_g - 1) * 3 + 2] += ldiag[e_l]

    return m_lump


nfields, field_fmt = read_header_file(sys.argv[1] + '.hdr')
X = read_file_np(sys.argv[1] + '.bin', nfields)
print(X)
# M = create_lumped_mass_matrix('M.out')
# M = np.diag(M)
# M_bar = np.linalg.cholesky(M)
M_bar = np.identity(nfields)
X_bar = np.dot(M_bar, X)
U, s, V = np.linalg.svd(X_bar, full_matrices=False)
# print(U.shape, V.shape, s.shape)
reduced_base = np.dot(np.linalg.inv(M_bar), U)
# print('B_r')
# print(B_r)
