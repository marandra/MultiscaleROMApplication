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
    i = -1
    for i, r in enumerate(elems):
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d} {:6d}\n".format(
                offset + i + 1, r[0] + 1, r[1] + 1, r[2] + 1, r[3] + 1
            )
        )
    fo.write("End Conditions\n\n")
    return i + 1


def write_conditions_triangle(fo, elems, offset=0):
    fo.write("Begin Conditions SurfaceCondition3D3N\n")
    i = -1
    for i, r in enumerate(elems):
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d}\n".format(
                offset + i + 1, r[0] + 1, r[1] + 1, r[2] + 1
            )
        )
    fo.write("End Conditions\n\n")
    return i + 1


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
    faces = []
    for cell_block in mesh.cells:
        element_type = cell_block[0]
        element_data = cell_block[1]
        if "hexa" in element_type:
            for e in element_data:
                faces.extend(skin_detect.explode_hexa(e))
        elif "wedge" in element_type:
            for e in element_data:
                faces.extend(skin_detect.explode_prism(e))
        skin_faces = skin_detect.filter_faces(faces)
        print(skin_faces)

    # end of my scrip
    condition_offset = []
    offset = 0
    for cell_block in mesh.cells:
        condition_type = cell_block[0]
        condition_data = cell_block[1]
        if "quad" in condition_type:
            condition_offset.append(offset)
            offset += write_conditions_quad(fo, condition_data, offset=offset)
        elif "triangle" in condition_type:
            condition_offset.append(offset)
            offset += write_conditions_triangle(fo, condition_data, offset=offset)

    # Groups (as submodelparts)
    for group_name, cell_arrays in mesh.cell_sets.items():
        group_cells = []
        for i, e in enumerate(cell_arrays):
            g = element_offset[i]
            elements = e + g
            group_cells.extend(elements)
        write_submodelpart(fo, group_name, cells=group_cells)

    #  Custom groups
    write_submodelpart(fo, "PINNED", points=[0])
    points = [x for x in range(len(mesh.points))]
    cells = [x for x in range(nr_cells)]
    write_submodelpart(fo, "RVE", points=points, cells=cells)
    points = [x for x in range(len(mesh.points))]  # FIX
    conditions = [x for x in range(nr_cells)]  # FIX
    write_submodelpart(fo, "SKIN", points=points, conditions=conditions)

    # fo.write("Begin SubModelPart SKIN\n")
    # fo.write("    Begin SubModelPartNodes\n")
    # # if we have triangles use numpy concatenate quad and triangles arrays
    # skin_nodes = set(mesh.cells["quad"].flatten())
    # for n in sorted(skin_nodes):
    #     fo.write("{:6d}\n".format(int(n)))
    # fo.write("    End SubModelPartNodes\n")
    # fo.write("    Begin SubModelPartElements\n")
    # fo.write("    End SubModelPartElements\n")
    # fo.write("    Begin SubModelPartConditions\n")
    # for i in range(max_i):
    #     fo.write("{:6d}\n".format(i + 1))
    # fo.write("    End SubModelPartConditions\n")
    # fo.write("End SubModelPart\n\n")
    #
    # for name, idx in mesh.field_data.items():
    #     if "SKIN" in name:
    #         continue
    #     group_array = mesh.cell_data["hexahedron"]["gmsh:physical"]
    #     fo.write("Begin SubModelPart {}\n".format(name))
    #     fo.write("    Begin SubModelPartNodes\n")
    #     fo.write("    End SubModelPartNodes\n")
    #     fo.write("    Begin SubModelPartElements\n")
    #     for e in numpy.where(group_array == idx[0])[0]:
    #         fo.write("{:6d}\n".format(e + 1))
    #     fo.write("    End SubModelPartElements\n")
    #     fo.write("    Begin SubModelPartConditions\n")
    #     fo.write("    End SubModelPartConditions\n")
    #     fo.write("End SubModelPart\n")
    #     fo.write("\n")
