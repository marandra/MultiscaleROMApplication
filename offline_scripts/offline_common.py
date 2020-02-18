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


def skip_calculation(filename, flag_reuse):
    try:
        with open(filename):
            flag_exists = True
    except IOError:
        flag_exists = False
    return flag_exists and flag_reuse


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
        "energy_inelastic_modes": 600,
        "strain_name": "STRAIN_FLUCTUANT",
        "strain_elastic_modes": 6,
        "strain_inelastic_modes": 100,
        "rvalue_name": "R_VALUE",
        "rvalue_elastic_modes": 1,
        "rvalue_inelastic_modes": 30,
        "rve_data_points": [200, -1],
        "rve_data_modes": [20],
        "reuse_existing_files": true,
        "svd_algorithm": "randomized",
        "training_path": "../training/",
        "training_rve_materials_fname": "materials.json",
        "training_case_pattern": "trajectory_{}",
        "offline_path": "../offline/",
        "bases_fname_pattern": "bases_{}_{}m.npy",
        "roc_fname_pattern": "roc_{}ip",
        "rve_fname_pattern": "rve_{}m_{}ip"
        }
        """
        )
        config = parameters["config_data"]
        config.ValidateAndAssignDefaults(config_defaults)

        if not check_consistent_config_values(config):
            exit()

        self.energy_name = config["energy_name"].GetString()
        self.energy_elastic_modes = config["energy_elastic_modes"].GetInt()
        self.energy_inelastic_modes = config["energy_inelastic_modes"].GetInt()
        self.energy_bases_fname = (
            config["bases_fname_pattern"]
            .GetString()
            .format(
                config["energy_name"].GetString(),
                config["energy_elastic_modes"].GetInt()
                + config["energy_inelastic_modes"].GetInt(),
            )
        )
        self.strain_name = config["strain_name"].GetString()
        self.strain_elastic_modes = config["strain_elastic_modes"].GetInt()
        self.strain_inelastic_modes = config["strain_inelastic_modes"].GetInt()
        self.strain_bases_fname = (
            config["bases_fname_pattern"]
            .GetString()
            .format(
                config["strain_name"].GetString(),
                config["strain_elastic_modes"].GetInt()
                + config["strain_inelastic_modes"].GetInt(),
            )
        )
        self.rvalue_name = config["rvalue_name"].GetString()
        self.rvalue_elastic_modes = config["rvalue_elastic_modes"].GetInt()
        self.rvalue_inelastic_modes = config["rvalue_inelastic_modes"].GetInt()
        self.rvalue_bases_fname = (
            config["bases_fname_pattern"]
            .GetString()
            .format(
                config["rvalue_name"].GetString(),
                config["rvalue_elastic_modes"].GetInt()
                + config["rvalue_inelastic_modes"].GetInt(),
            )
        )
        self.svd_algorithm = config["svd_algorithm"].GetString()
        self.reuse_existing_files = config["reuse_existing_files"].GetBool()

        self.ip_subsets = config["rve_data_points"]
        self.roc_fname_pattern = config["roc_fname_pattern"].GetString()

    def roc_fname(self, p):
        if p==-1:
            return self.roc_fname_pattern.format("ROM")
        else:
            return self.roc_fname_pattern.format(p)

if __name__ == "__main__":
    print(Common().energy_name)
    print(Common().energy_bases_fname)
    print(Common().strain_name)
    print(Common().strain_bases_fname)
    print(Common().rvalue_name)
    print(Common().rvalue_bases_fname)
    print(Common().roc_fname(1))
    print(Common().roc_fname(100))
    print(Common().roc_fname(1000))
    print(Common().roc_fname(-1))
