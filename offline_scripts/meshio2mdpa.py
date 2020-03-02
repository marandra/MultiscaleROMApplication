import meshio
import sys
import numpy
import skin_detect

def write_header(fo):
    fo.write("Begin ModelPartData\n")
    fo.write("End ModelPartData\n")
    fo.write("\n")
    fo.write("Begin Properties 0\n")
    fo.write("End Properties\n")
    fo.write("\n")


def write_points(fo, points):
    fo.write("Begin Nodes\n")
    for i, p in enumerate(points):
        fo.write("{:6d}   {:6f} {:6f} {:6f}\n".format(i + 1, p[0], p[1], p[2]))
    fo.write("End Nodes\n")
    fo.write("\n")


def write_elements_hexahedron(fo, elems, offset=0):
    fo.write("Begin Elements SmallDisplacementElement3D8N\n")
    i = 0
    for i0, p0 in enumerate(elems):
        i = i0 + 1  # We start elements by 1
        p = p0 + 1  # We start nodes by 1
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d} {:6d} {:6d} {:6d} {:6d} {:6d}\n".format(
                offset + i, *p
            )
        )
    fo.write("End Elements\n")
    fo.write("\n")
    return i


def write_elements_wedge(fo, elems, offset=0):
    fo.write("Begin Elements SmallDisplacementElement3D6N\n")
    i = 0
    for i0, p0 in enumerate(elems):
        i = i0 + 1  # We start elements by 1
        p = p0 + 1  # We start nodes by 1
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d} {:6d} {:6d} {:6d}\n".format(offset + i, *p)
        )
    fo.write("End Elements\n")
    fo.write("\n")
    return i


def write_conditions_quad(fo, elems, offset=0):
    fo.write("Begin Conditions SurfaceCondition3D4N\n")
    i = 0
    for i0, p0 in enumerate(elems):
        i = i0 + 1  # We start elements by 1
        p = p0 + 1  # We start elements by 1
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d} {:6d}\n".format(
                offset + i, *p
            )
        )
    fo.write("End Conditions\n\n")
    return i


def write_conditions_triangle(fo, elems, offset=0):
    fo.write("Begin Conditions SurfaceCondition3D3N\n")
    i = 0
    for i0, p0 in enumerate(elems):
        i = i0 + 1  # We start elements by 1
        p = p0 + 1  # We start elements by 1
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d}\n".format(
                offset + i, *p
            )
        )
    fo.write("End Conditions\n\n")
    return i


def write_submodelpart(fo, group_name, points=[], cells=[], conditions=[]):
    fo.write("Begin SubModelPart {}\n".format(group_name))
    fo.write("    Begin SubModelPartNodes\n")
    for p0 in points:
        p = p0 + 1
        fo.write("        {:6d}\n".format(p))
    fo.write("    End SubModelPartNodes\n")
    fo.write("    Begin SubModelPartElements\n")
    for c0 in cells:
        c = c0 + 1
        fo.write("        {:6d}\n".format(c))
    fo.write("    End SubModelPartElements\n")
    fo.write("    Begin SubModelPartConditions\n")
    for c0 in conditions:
        c = c0 + 1
        fo.write("        {:6d}\n".format(c))
    fo.write("    End SubModelPartConditions\n")
    fo.write("End SubModelPart\n\n")


######################################################
# main
######################################################

i_filename = sys.argv[1]
o_filename = sys.argv[2]
mesh = meshio.read(i_filename)
print(mesh)
print("")

print("*********** WIP notice *************")
print("Elements currently implemented:")
print("   hexahedron")
print("   wedge")
print("Conditions currently implemented:")
print("   quad")
print("   triangle")
print("************************************")
with open(o_filename, "w") as fo:

    #  Header
    write_header(fo)

    #  Nodes
    write_points(fo, mesh.points)

    #  Elements
    element_offset = []
    offset = 0
    for cell_block in mesh.cells:
        element_type = cell_block[0]
        element_data = cell_block[1]
        if "hexa" in element_type:
            element_offset.append(offset)
            offset += write_elements_hexahedron(fo, element_data, offset=offset)
        elif "wedge" in element_type:
            element_offset.append(offset)
            offset += write_elements_wedge(fo, element_data, offset=offset)
    nr_cells = offset

    # Conditions. Surface elements obtained exploiding volume elements
    # my script: generating skin
    quads = []
    triangles = []
    for cell_block in mesh.cells:
        element_type = cell_block[0]
        element_data = cell_block[1]
        if "hexa" in element_type:
            for e in element_data:
                quads.extend(skin_detect.explode_hexa(e))
        elif "wedge" in element_type:
            for e in element_data:
                qs, ts = skin_detect.explode_prism(e)
                quads.extend(qs)
                triangles.extend(ts)
    skin_quads = skin_detect.filter_faces(quads)
    skin_triangles = skin_detect.filter_faces(triangles)
    cb_q = meshio.CellBlock("quads", numpy.array(skin_quads))
    cb_t = meshio.CellBlock("triangles", numpy.array(skin_triangles))
    skin_cells = [cb_q, cb_t]
    # end of my scrip

    condition_offset = []
    offset = 0
    for cell_block in skin_cells:
        condition_type = cell_block[0]
        condition_data = cell_block[1]
        if "quad" in condition_type:
            condition_offset.append(offset)
            offset += write_conditions_quad(fo, condition_data, offset=offset)
        elif "triangle" in condition_type:
            condition_offset.append(offset)
            offset += write_conditions_triangle(fo, condition_data, offset=offset)
    nr_conditions = offset

    # Groups (as submodelparts)
    for group_name, cell_arrays in mesh.cell_sets.items():
        group_cells = []
        for i, e in enumerate(cell_arrays):
            g = element_offset[i]
            elements = e + g
            group_cells.extend(elements)
        write_submodelpart(fo, group_name, cells=group_cells)

    #  Custom groups
    # TODO: add following submodelparts as groups in mesh.cell_sets, so we avoid these write_submodelparts()
    p_origin = numpy.where((mesh.points[:,0]==0) * (mesh.points[:,1]==0) * (mesh.points[:,2]==0))[0]
    write_submodelpart(fo, "PINNED", p_origin)
    points = [x for x in range(len(mesh.points))]
    cells = [x for x in range(nr_cells)]
    write_submodelpart(fo, "RVE", points=points, cells=cells)
    aux_points = []
    for cb in skin_cells:
        for a in cb[1]:
            aux_points.extend(a)
    skin_points = set(aux_points)
    points = [x for x in skin_points]
    conditions = [x for x in range(nr_conditions)]
    write_submodelpart(fo, "SKIN", points=points, conditions=conditions)
