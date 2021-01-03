"""Compute error between datasets of HF and reconstructed fields.
"""
import sys
import numpy
import meshio

F1 = sys.argv[1]
F2 = sys.argv[2]

with meshio.xdmf.TimeSeriesReader(F1) as reader:
    points, cells = reader.read_points_cells()
    for k in range(reader.num_steps):
        t_1, point_data_1, cell_data_1 = reader.read_data(k)

with meshio.xdmf.TimeSeriesReader(F2) as reader:
    points, cells = reader.read_points_cells()
    for k in range(reader.num_steps):
        t_2, point_data_2, cell_data_2 = reader.read_data(k)

npoint = len(points)
ncells = 0
for c in cells:
    ncells += len(c.data)
print("Nr points: ", npoint)
print("Nr cells: ", ncells)

#ord = numpy.Inf
ord = 1
for key in point_data_1:
    a = point_data_1[key]
    b = point_data_2[key]
    ren = numpy.linalg.norm(a - b, ord=ord, axis=0)
    red = numpy.linalg.norm(a, ord=ord)
    re = numpy.max(ren / red * 100.0)
    print(key, "{:.3}%".format(re))

#print(cell_data_1)
for key in cell_data_1:
    A = cell_data_1[key]
    #print(len(A))
    #print(numpy.shape(A[0][0]))
    #print((A[0][0]))
    #print(numpy.shape(A[1]))
    a = numpy.array(numpy.concatenate(cell_data_1[key]))
    b = numpy.array(numpy.concatenate(cell_data_2[key]))
    ren = numpy.linalg.norm(a - b, ord=ord, axis=0)
    red = numpy.linalg.norm(a, ord=ord)
    re = numpy.max(ren / red * 100.0)
    #re = ren / red * 100.0
    #print(key, "{:.3}%".format(ren))
    #print(key, "{:.3}%".format(red))
    print(key, "{:.3}%".format(re))
    #print(numpy.shape(a))
    #print(numpy.shape(b))

