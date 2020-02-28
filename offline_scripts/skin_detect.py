import numpy
import meshio
import sys
import pprint as pp


def explode_hexa(e):
    """return quadrilaterals of the input hexahedron"""
    q0 = [e[0], e[1], e[5], e[4]]
    q1 = [e[1], e[2], e[6], e[5]]
    q2 = [e[2], e[3], e[7], e[6]]
    q3 = [e[3], e[0], e[4], e[7]]
    q4 = [e[4], e[5], e[6], e[7]]
    q5 = [e[0], e[3], e[2], e[1]]
    return [q0, q1, q2, q3, q4, q5]


def explode_prism(e):
    """return quadrilaterals and triangles of the input prism"""
    q0 = [e[0], e[1], e[4], e[3]]
    q1 = [e[1], e[2], e[5], e[4]]
    q2 = [e[2], e[0], e[3], e[5]]
    t0 = [e[3], e[4], e[5]]
    t1 = [e[0], e[2], e[1]]
    return [q0, q1, q2], [t0, t1]


def explode_tetra(e):
    """return triangles of the input tetrahedron"""
    t0 = [e[0], e[2], e[1]]
    t1 = [e[0], e[1], e[3]]
    t2 = [e[1], e[2], e[3]]
    t3 = [e[2], e[0], e[3]]
    return [t0, t1, t2, t3]


def filter_faces(candidates):
    """ deletes repeted faces, as skin faces will be unique"""
    # we use a dictionary to detect duplicated faces, as they
    # will have the same key regardless of their orientation
    skin = {}
    for es in candidates:
        key = "{}".format(sorted(es))
        # print("debug try: ", key, es)
        if key in skin:
            del skin[key]
            print('debug:    deleted faces with key "{}"'.format(key))
        else:
            skin[key] = es
    # we return a list, as we do not need keys anymore
    return list(skin.values())

####################################################    
# main
####################################################    
if __name__=="__main":
    i_filename = sys.argv[1]
    o_filename = sys.argv[2]
    s_filename = sys.argv[3]
    # o_name = sys.argv[1].rsplit(".", 1)[0]
    # o_ext = sys.argv[1].rsplit(".")[-1]
    # o_filename = "{}_skin.{}".format(o_name, o_ext)
    mesh = meshio.read(i_filename)
    meshio.write(o_filename, mesh)
    
    candidates_t = []
    candidates_q = []
    # hexas = [[1, 2, 3, 4, 5, 6, 7, 8], [4, 3, 9, 10, 8, 7, 11, 12]]
    # prism1 = [[21, 22, 23, 24, 25, 26], [27, 28, 29, 21, 22, 23]]
    # prism2 = [[31, 32, 33, 34, 35, 36], [37, 31, 33, 38, 34, 36]]
    # tetras = [[41, 42, 43, 44], [42, 45, 43, 44]]
    hexas = mesh.cells["hexahedron"]
    print(hexas)
    stop
    prisms = mesh.cells["wedge"]
    print(prisms)
    stop
    tetras = mesh.cells["tetra"]
    print(tetras)
    stop
    
    for ev in hexas:
        quads = explode_hexa(ev)
        candidates_q.extend(quads)
    
    for ev in prisms:
        quads, triangles = explode_prism(ev)
        candidates_q.extend(quads)
        candidates_t.extend(triangles)
    
    for ev in tetras:
        triangles = explode_tetra(ev)
        candidates_t.extend(triangles)
    
    skin_eq = filter_faces(candidates_q)
    skin_et = filter_faces(candidates_t)
    
    pp.pprint(skin_eq)
    pp.pprint(skin_et)
    
    skin_cells = {"triangle": numpy.array(skin_et), "quad": numpy.array(skin_eq)}
    
    meshio.write(s_filename, mesh.points, skin_cells)
