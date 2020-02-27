import meshio
import sys
import numpy


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
    i = -1
    for i, p in enumerate(elems):
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d} {:6d} {:6d} {:6d} {:6d} {:6d}\n".format(
                offset + i + 1, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7]
            )
        )
    fo.write("End Elements\n")
    fo.write("\n")
    return i + 1


def write_elements_wedge(fo, elems, offset=0):
    fo.write("Begin Elements SmallDisplacementElement3D6N\n")
    i = -1
    for i, p in enumerate(elems):
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d} {:6d} {:6d} {:6d}\n".format(
                offset + i + 1, p[0], p[1], p[2], p[3], p[4], p[5]
            )
        )
    fo.write("End Elements\n")
    fo.write("\n")
    return i + 1


def write_conditions_quad(fo, elems, offset=0):
    fo.write("Begin Conditions SurfaceCondition3D4N\n")
    i = -1
    for i, r in enumerate(elems):
        fo.write(
            "{:6d}  0  {:6d} {:6d} {:6d} {:6d}\n".format(offset + i + 1, r[0], r[1], r[2], r[3])
        )
    fo.write("End Conditions\n\n")
    return i + 1


def write_conditions_triangle(fo, elems, offset=0):
    fo.write("Begin Conditions SurfaceCondition3D3N\n")
    i = -1
    for i, r in enumerate(elems):
       fo.write(
           "{:6d}  0  {:6d} {:6d} {:6d}\n".format(offset + i + 1, r[0], r[1], r[2])
       )
    fo.write("End Conditions\n\n")
    return i + 1

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


    write_header(fo)

    #write_points(fo, mesh.points)

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

    for group_name, element_list in mesh.cell_sets.items():
        name = cell_group
        write_submodelpart_elements(fo, group_name)



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
    #
    # fo.write("Begin SubModelPart PINNED\n")
    # fo.write("    Begin SubModelPartNodes\n")
    # fo.write("        {:6d}\n".format(0))
    # fo.write("    End SubModelPartNodes\n")
    # fo.write("    Begin SubModelPartElements\n")
    # fo.write("    End SubModelPartElements\n")
    # fo.write("    Begin SubModelPartConditions\n")
    # fo.write("    End SubModelPartConditions\n")
    # fo.write("End SubModelPart\n\n")
    #
    # fo.write("Begin SubModelPart RVE\n")
    # fo.write("    Begin SubModelPartNodes\n")
    # for i, p in enumerate(mesh.points):
    #     fo.write("        {:6d}\n".format(i))
    # fo.write("    End SubModelPartNodes\n")
    # fo.write("    Begin SubModelPartElements\n")
    # for i, p in enumerate(mesh.cells["hexahedron"]):
    #     fo.write("        {:6d}\n".format(i + 1))
    # fo.write("    End SubModelPartElements\n")
    # fo.write("    Begin SubModelPartConditions\n")
    # # for i in range(len(skin_eq) + len(skin_et)):
    # #    fo.write("       {:6d}\n".format(i + 1))
    # fo.write("    End SubModelPartConditions\n")
    # fo.write("End SubModelPart\n\n")
