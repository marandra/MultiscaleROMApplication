"""Offline pipeline.

Usage:
    program.py [-h]
    program.py [-c FILE] [COMMAND]

Options:
-h --help           Show this
-c --config FILE    Specify configuration file [default: ../ProjectsParameters.json]
                    The location of the configuration file is taken as root path

Commands:
    dump_config     Dump a configuration file with default values, overwritten by
                    user values if option "-c FILE" is passed.
    test            Write configuration values generated automatically.

"""

# Check schema for input validation https://github.com/keleshev/schema
# Check fire for exposing objects to CLI https://github.com/google/python-fire
# Check docopt for args parsing https://github.com/docopt/docopt

import json
from pathlib import Path
import fire
from docopt import docopt


def validate_context(default, user):
    """
    Validates and merges defaults and user configurations.
    Received configeration dictionaries
    """
    # Ideas for validation:
    # number of base modes < number of snapshots
    # number of base mode > number of requested modes

    # all keys in user must be present in default
    d_k = [dk for dk in default.keys()]
    for u_k in user.keys():
        if u_k not in d_k:
            raise SystemExit(
                'Not recognized key "{}" in user configuration. Exit.'.format(u_k)
            )

    # create context by merging default and user
    context = {**default, **user}
    return context


#######################################################################
#######################################################################


class Common:
    """
    TODO add docstrings
    """

    def __init__(self, config_fname="../ProjectParameters.json"):
        try:
            context_user = json.loads(Path(config_fname).read_text())["config_data"]
        except FileNotFoundError:
            print("WARNING: No such configuration file: '{}'".format(config_fname))
            context_user = {}
        context_defaults = {
            # most frequently set
            "cases_test_dataset": [0],
            "rve_data_points": [100],
            "rve_data_points_range": [100, 500, 100],
            "rve_data_points_rom": True,
            "rve_data_modes": [20],
            #
            "energy_name": "ENERGY_FREE",
            "energy_elastic_modes": 21,
            "energy_inelastic_modes": -1,
            "energy_svd_cutoff": 1e-4,
            "strain_name": "STRAIN_FLUCTUANT",
            "strain_elastic_modes": 6,
            "strain_inelastic_modes": -1,
            "strain_svd_cutoff": 1e-4,
            "rvalue_name": "R_VALUE",
            "rvalue_elastic_modes": 1,
            "rvalue_inelastic_modes": 30,
            "rvalue_svd_cutoff": 1e-4,
            "reuse_existing_files": True,
            # training files stuff
            "training_path": "training",
            "training_rve_materials_fname": "materials.json",
            "case_path_pattern": "case_{}",
            "snapshots_fname": "snapshots.hdf5",
            "training_strain_fname": "_training_strain_set.dat",
            # offline files stuff
            "offline_path": "offline_data",
            "bases_fname_pattern": "bases_{}_{}m.npy",
            "local_bases_fname_pattern": "bases_inelastic_local_{}.npy",
            "local_sv_fname_pattern": "sv_inelastic_local_{}.dat",
            "roc_fname_pattern": "roc_{}ip",
            "rve_fname_pattern": "rve_{}m_{}ip.json",
            # multiscale files stuff
            "multiscale_path": "multiscale_1ip",
            # other files stuff
        }
        config = validate_context(context_defaults, context_user)
        self.context = config
        # TODO:
        # load defaults
        # set default user location, overwrite with args
        # load user config
        # update default config with user config

        # file management
        self.root_path = Path.cwd() / Path(config_fname).parent
        self.training_path = self.root_path / config["training_path"]
        self.offline_path = self.root_path / config["offline_path"]
        self.multiscale_path = self.root_path / config["multiscale_path"]

        self.reuse_existing_files = config["reuse_existing_files"]
        self.bases_fname = config["bases_fname_pattern"]
        self.local_bases_fname = config["local_bases_fname_pattern"]
        self.local_sv_fname = config["local_sv_fname_pattern"]

        # bases generation
        self.svd_cutoff = {}

        # self.energy_name = config["energy_name"]
        self.energy_elastic_modes = config["energy_elastic_modes"]
        self.energy_inelastic_modes = config["energy_inelastic_modes"]
        self.svd_cutoff[config["energy_name"]] = config["energy_svd_cutoff"]

        # self.strain_name = config["strain_name"]
        self.strain_elastic_modes = config["strain_elastic_modes"]
        self.strain_inelastic_modes = config["strain_inelastic_modes"]
        self.svd_cutoff[config["strain_name"]] = config["strain_svd_cutoff"]

        # self.rvalue_name = config["rvalue_name"]
        self.rvalue_elastic_modes = config["rvalue_elastic_modes"]
        self.rvalue_inelastic_modes = config["rvalue_inelastic_modes"]
        self.svd_cutoff[config["rvalue_name"]] = config["rvalue_svd_cutoff"]

        # points
        self.ip_subsets = config["rve_data_points"]
        for i in range(*config["rve_data_points_range"]):
            self.ip_subsets.append(i)
        self.ip_subsets.sort()
        if config["rve_data_points_rom"]:
            self.ip_subsets.append("ROM")

        self.roc_fname_pattern = config["roc_fname_pattern"]

        # modes
        # self.reduced_nr_modes = config["rve_data_modes"]

        self.materials_fname = self.training_path / Path(
            config["training_rve_materials_fname"]
        )
        self.rve_fname_pattern = config["rve_fname_pattern"]

    # def get_cases_rve_paths(self, case_id):
    #     """
    #     Returns list of paths to 1ip multiscale cases for case: case_id
    #     """
    #     paths = []
    #     for mode in self.context["rve_data_modes"]:
    #         for point in self.ip_subsets:
    #             paths.append(
    #                 (
    #                     self.multiscale_path
    #                     / self.case_name(case_id)
    #                     / "_{}m_{}ip".format(mode, point)
    #                 ).resolve()
    #             )

    def dump_config(self, fname="ProjectParameters.json"):
        Path(fname).write_text(json.dumps(self.context, indent=2))
        print("Written configuration file {}. Move it to the root path.".format(fname))

    def parse_training_strain_set(self):
        """
        Returns list of strain vectors used for trainig, read from file defined in configuration
        """
        fpath = self.training_path / self.context["training_strain_fname"]
        return fpath.read_text().splitlines()

    def case_name(self, c_id):
        """
        Returns case name with corresponging leading zeros
        e.g. if nr_cases:100, id:00..99, nr_id: 2 -> case_01..case_99
        """
        strain_set = self.parse_training_strain_set()
        nr_cases = len(strain_set)
        len_id = len(str(nr_cases - 1))  # size of the case number string
        case_id = "{:0{}d}".format(c_id, len_id)
        case_name = self.context["case_path_pattern"].format(case_id)
        return case_name

    def roc_fname(self, points):
        """
        docstrings here
        """
        return self.roc_fname_pattern.format(points)

    def rve_fname(self, modes, points):
        """
        docstrings here
        """
        return self.rve_fname_pattern.format(modes, points)

    def skip_calculation(self, fname):
        """ 
        Generates a list of files following filename pattern.
        Length of list is used as flag (False if empty, True otherwise)
        """
        fpath = Path.cwd() / fname  # converts filename to absolute Path
        flag_exists = len([f for f in fpath.parent.glob(fpath.name)])
        flag_reuse = self.context["reuse_existing_files"]
        return flag_exists and flag_reuse

    def get_bases_fname(self, field):
        """
        docstrings here
        """
        filename = self.bases_fname.format(field, "*")
        fpath = self.offline_path / filename
        files = [f for f in fpath.parent.glob(fpath.name)]
        if len(files) == 0:
            return None
        if len(files) > 1:
            print(
                "Warning: More than one {} bases file detected. "
                "Picking first in the list: {}".format(field, files[0].name)
            )
        return files[0]


#####################################################################
# main
#####################################################################

if __name__ == "__main__":
    arguments = docopt(__doc__, version="2.0")
    # import pprint
    # pprint.pprint(arguments)

    if "--config" in arguments:
        C = Common(config_fname=arguments["--config"])
    else:
        C = Common()

    # parse command line commands
    if arguments["COMMAND"] is not None:
        if "dump_config" in arguments["COMMAND"]:
            C.dump_config()
            exit()
        if "test" in arguments["COMMAND"]:
            print("Test:")
            print(C.bases_fname)
            print(C.local_bases_fname)
            print(C.roc_fname("1"))
            print(C.roc_fname("100"))
            print(C.roc_fname("1000"))
            print(C.roc_fname(1000))
            print(C.roc_fname("ROM"))
            print(C.rve_fname(20, "1"))
            print(C.rve_fname(20, "100"))
            print(C.rve_fname(200, "1000"))
            print(C.rve_fname(200, 1000))
            print(C.rve_fname("2000", "ROM"))
            print(C.ip_subsets)
            # print(C.case_name(5))
