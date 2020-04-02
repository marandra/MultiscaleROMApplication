import glob
from pathlib import Path
import KratosMultiphysics


"""
TODO: pending description here.
"""


def check_consistent_config_values(config):
    # TODO: Use this function
    # Ideas:
    # number of base modes < number of snapshots
    # number of base mode > number of requested modes
    check = True
    return check


#######################################################################
#######################################################################


class Common:
    def __init__(self, config_fname="../configuration.json"):
        with open(config_fname, "r") as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        config_defaults = KratosMultiphysics.Parameters(
            """{
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
        "reuse_existing_files": true,
        "rve_data_points": [200, -1],
        "rve_data_modes": [20],
        "training_path": "../training/",
        "training_rve_materials_fname": "materials.json",
        "training_case_pattern": "trajectory_{}",
        "offline_path": "../offline/",
        "bases_fname_pattern": "bases_{}_{}m.npy",
        "local_bases_fname_pattern": "bases_inelastic_local_{}.npy",
        "local_sv_fname_pattern": "sv_inelastic_local_{}.dat",
        "roc_fname_pattern": "roc_{}ip",
        "rve_fname_pattern": "rve_{}m_{}ip.json",
        "skip_cases_from_training": []
        }
        """
        )
        config = parameters["config_data"]
        config.ValidateAndAssignDefaults(config_defaults)

        if not check_consistent_config_values(config):
            exit()

        # file management
        self.reuse_existing_files = config["reuse_existing_files"].GetBool()
        self.bases_fname = config["bases_fname_pattern"].GetString()
        self.local_bases_fname = config["local_bases_fname_pattern"].GetString()
        self.local_sv_fname = config["local_sv_fname_pattern"].GetString()
        self.skip_cases = []
        for i in config["skip_cases_from_training"]:
            self.skip_cases.append(i.GetInt())

        # bases generation
        self.svd_cutoff = {}

        self.energy_name = config["energy_name"].GetString()
        self.energy_elastic_modes = config["energy_elastic_modes"].GetInt()
        self.energy_inelastic_modes = config["energy_inelastic_modes"].GetInt()
        self.svd_cutoff[self.energy_name] = config["energy_svd_cutoff"].GetDouble()

        self.strain_name = config["strain_name"].GetString()
        self.strain_elastic_modes = config["strain_elastic_modes"].GetInt()
        self.strain_inelastic_modes = config["strain_inelastic_modes"].GetInt()
        self.svd_cutoff[self.strain_name] = config["strain_svd_cutoff"].GetDouble()

        self.rvalue_name = config["rvalue_name"].GetString()
        self.rvalue_elastic_modes = config["rvalue_elastic_modes"].GetInt()
        self.rvalue_inelastic_modes = config["rvalue_inelastic_modes"].GetInt()
        self.svd_cutoff[self.rvalue_name] = config["rvalue_svd_cutoff"].GetDouble()

        # ROC and the rest
        self.ip_subsets = config["rve_data_points"]
        self.roc_fname_pattern = config["roc_fname_pattern"].GetString()

        self.materials_fname = (
            config["training_path"].GetString()
            + config["training_rve_materials_fname"].GetString()
        )
        self.reduced_nr_modes = config["rve_data_modes"]
        self.rve_fname_pattern = config["rve_fname_pattern"].GetString()

    def roc_fname(self, p):
        if p == -1:
            return self.roc_fname_pattern.format("ROM")
        else:
            return self.roc_fname_pattern.format(p)

    def rve_fname(self, m, p):
        if p == -1:
            return self.rve_fname_pattern.format(m, "ROM")
        else:
            return self.rve_fname_pattern.format(m, p)

    def skip_calculation(self, fname, flag_reuse):
        """ 
        Generates a list of files following filename pattern.
        Length of list is used as flag (False if empty, True otherwise)
        """
        fpath = Path.cwd() / fname  # converts filename to absolute Path
        flag_exists = len([f for f in fpath.parent.glob(fpath.name)])
        return flag_exists and flag_reuse

    def get_bases_fname(self, field):
        filename = self.bases_fname.format(field, "*")
        files = glob.glob(filename)
        if not len(files):
            return None
        if len(files) > 1:
            print(
                "Warning: More than one {} bases file detected. Picking first in the list: {}".format(
                    field, files[0]
                )
            )
        return files[0]


#####################################################################
#
#####################################################################

if __name__ == "__main__":
    print("Test:")
    print(Common().energy_name)
    print(Common().strain_name)
    print(Common().rvalue_name)
    print(Common().bases_fname)
    print(Common().local_bases_fname)
    print(Common().roc_fname(1))
    print(Common().roc_fname(100))
    print(Common().roc_fname(1000))
    print(Common().roc_fname(-1))
    print(Common().rve_fname(20, 1))
    print(Common().rve_fname(20, 100))
    print(Common().rve_fname(200, 1000))
    print(Common().rve_fname(2000, -1))
    print(Common().skip_cases)
    print(Common().get_bases_fname(Common().energy_name))
    print(Common().get_bases_fname(Common().strain_name))
    print(Common().get_bases_fname(Common().rvalue_name))
