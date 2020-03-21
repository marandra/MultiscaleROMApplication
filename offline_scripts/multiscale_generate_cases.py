import KratosMultiphysics
import json
from pathlib import Path, PurePath
from offline_common import Common


def create_properties_file(props_fname, training_path, traj_name, case_path):
    props_path = Path(props_fname)
    training_props_path = PurePath.joinpath(
        training_path, Path("trajectory_{}".format(traj_name[1:][:-1])), props_path
    )
    with training_props_path.open() as fi:
        props = json.load(fi)
        strain_versor = props["processes"]["loads_process_list"][0]["Parameters"][
            "initial_strain"
        ]
        ampl = props["processes"]["loads_process_list"][0]["Parameters"][
            "lookuptable_mult"
        ][-1]
    with props_path.open() as fi:
        props = json.load(fi)
        # compute displacements u = E * x
        s0, s1, s2, s3, s4, s5 = strain_versor
        x0 = 1.0 * s0 * ampl
        y0 = 0.5 * s3 * ampl
        z0 = 0.5 * s5 * ampl
        x1 = 0.5 * s3 * ampl
        y1 = 1.0 * s1 * ampl
        z1 = 0.5 * s4 * ampl
        x2 = 0.5 * s5 * ampl
        y2 = 0.5 * s4 * ampl
        z2 = 1.0 * s2 * ampl
        props["processes"]["list_boundary_processes"][1]["Parameters"]["value"] = [
            "{}*t".format(x0),
            "{}*t".format(x1),
            "{}*t".format(x2),
        ]
        props["processes"]["list_boundary_processes"][2]["Parameters"]["value"] = [
            "{}*t".format(y0),
            "{}*t".format(y1),
            "{}*t".format(y2),
        ]
        props["processes"]["list_boundary_processes"][3]["Parameters"]["value"] = [
            "{}*t".format(z0),
            "{}*t".format(z1),
            "{}*t".format(z2),
        ]
    with PurePath.joinpath(case_path, props_path).open("w") as fo:
        json.dump(props, fo, indent=4)
    return


def create_materials_file(materials_fname, offline_path, rve_name, case_path):
    materials_path = Path(materials_fname)
    rve_data_path = PurePath.joinpath(offline_path, Path("rve{}.json".format(rve_name)))
    with materials_path.open() as fi:
        materials = json.load(fi)
        materials["properties"][0]["Material"]["constitutive_law"]["Parameters"][
            "rve_data_filename"
        ] = str(rve_data_path)
    with open(PurePath.joinpath(case_path, materials_path), "w") as fo:
        json.dump(materials, fo, indent=4)
    return


def copy_file(filename, case_path):
    src = Path(filename)
    dest = PurePath.joinpath(case_path, src)
    dest.write_text(src.read_text())
    return


def create_case_dir(traj_name, rve_name, training_path, offline_path):
    case_path = Path("case" + traj_name + rve_name)
    case_path.mkdir(exist_ok=True)
    offline_path = PurePath.joinpath(Path(".."), offline_path)
    # populate case dir
    create_materials_file("macro_materials.json", offline_path, rve_name, case_path)
    create_properties_file(
        "ProjectParameters.json", training_path, traj_name, case_path
    )
    copy_file("MainKratos.py", case_path)
    copy_file("macro_model.mdpa", case_path)
    return


#######################################
# main
#######################################


if __name__ == "__main__":

    with open("../configuration.json", "r") as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
    # configuration
    trajectories = [99]
    points = Common().ip_subsets
    modes = Common().reduced_nr_modes
    #modes = [10, 20, 30, 40, 50, 60, 70, 80]
    #points = [50, 100, 200, "ROM"]
    validation_path = Path("../training/validation")
    offline_path = Path("../offline_data")

    for t in trajectories:
        for mp in modes:
            m = mp.GetInt()
            for pp in points:
                p = pp.GetInt()
                if p == -1:
                    p = "ROM"
                rve_name = "_{}m_{}ip".format(m, p)
                traj_name = "_{}t".format(t)
                print(rve_name)
                create_case_dir(traj_name, rve_name, validation_path, offline_path)
