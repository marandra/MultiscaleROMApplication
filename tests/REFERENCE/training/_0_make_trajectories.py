import json
from pathlib import Path, PurePath


def create_properties_file(props_fname, training_path, traj_name, case_path, strain):
    props_path = Path(props_fname)
    case_props_path = PurePath.joinpath(training_path, case_path, props_path)
    with props_path.open() as fi:
        props = json.load(fi)
        props["processes"]["loads_process_list"][0]["Parameters"][
            "initial_strain"
        ] = strain
    with PurePath.joinpath(case_path, props_path).open("w") as fo:
        json.dump(props, fo, indent=4)
    return


def copy_file(filename, case_path):
    src = Path(filename)
    dest = PurePath.joinpath(case_path, src)
    dest.write_text(src.read_text())
    return


def create_case_dir(traj_name, training_path, strain):
    case_path = Path("trajectory" + traj_name)
    case_path.mkdir(exist_ok=True)
    # populate case dir
    create_properties_file(
        "ProjectParameters.json", training_path, traj_name, case_path, strain
    )
    copy_file("MainKratos.py", case_path)
    copy_file("model.mdpa", case_path)
    copy_file("materials.json", case_path)
    return


######################################
# main
######################################

if __name__ == "__main__":
    strain_set_fname = "_training_strain_set.dat"
    with open(strain_set_fname, "r") as fi:
        strain_set = fi.readlines()
        nr_cases = len(strain_set)
        #  size of the case number string (e.g. nr_cases=100, id: 00..99, nr_id=2)
        nr_id = len(str(nr_cases - 1))
        for i, line in enumerate(strain_set):
            traj_name = "_{:0{}d}".format(i, nr_id)
            strain = [float(x) for x in line.split()]
            print(traj_name, strain)
            create_case_dir(traj_name, Path("."), strain)
    validation_path = Path("validation")
    validation_path.mkdir(exist_ok=True)
