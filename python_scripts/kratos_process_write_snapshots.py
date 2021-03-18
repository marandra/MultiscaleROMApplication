from pathlib import Path
import h5py
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as sm
from hprfe2 import bases
from hprfe2.common import Common


def Factory(settings, model):
    return WriteSnapshots(settings["Parameters"], model)


class WriteSnapshots(km.Process):
    def __init__(self, settings, model):
        km.Process.__init__(self)

        default_settings = km.Parameters(
            """
            {"model_part_name": "unset_model_part_name",
             "material_root_path": "."}
            """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        # The right root_path is not critical in this case, if we only use default
        # values, i.e., don't have modified field or files names in the config file
        C = Common(Path(settings["material_root_path"].GetString()))
        self.config = C.config

        self.model_part = model[settings["model_part_name"].GetString()]
        self.timestep_counter = 1
        self.inelastic_flag = False

    def has_damaged_elements(self):
        for elem in self.model_part.Elements:
            flag = elem.CalculateOnIntegrationPoints(
                # Only for damage CLs for now
                km.DAMAGE_VARIABLE, self.model_part.ProcessInfo
            )
            if True in [x > 0. for x in flag]:
                return True
        return False

    def write_strain(self, group, filename, timestep):
        data_list = []
        for elem in self.model_part.Elements:
            strain = elem.CalculateOnIntegrationPoints(
                km.INITIAL_STRAIN_VECTOR,
                self.model_part.ProcessInfo,
            )
            break
        strain_macro = strain[0] * (-1)  # init strain was passed as negative
        for elem in self.model_part.Elements:
            strain_vectors = elem.CalculateOnIntegrationPoints(
                km.STRAIN,
                self.model_part.ProcessInfo,
            )
            for strain_ip in strain_vectors:
                for i, strain_i in enumerate(strain_ip):
                    strain_fluctuant_i = strain_i - strain_macro[i]
                    data_list.append(strain_fluctuant_i)
        bases.write_field_to_hdf5(
            filename, group, self.config["strain_name"], timestep, data_list
        )

    def write_energy(self, group, filename, timestep):
        data_list = []
        for elem in self.model_part.Elements:
            strain_energy_values = elem.CalculateOnIntegrationPoints(
                km.STRAIN_ENERGY, self.model_part.ProcessInfo
            )
            for strain_energy_ip in strain_energy_values:
                data_list.append(strain_energy_ip)
        bases.write_field_to_hdf5(
            filename, group, self.config["energy_name"], timestep, data_list
        )

    def write_rvalue(self, group, filename, timestep):
        data_list = []
        for elem in self.model_part.Elements:
            # TODO: Check if just skipping elements without internal variables
            # (e.g. LinearElastic3D) works
            #if not elem.Has(km.INTERNAL_VARIABLES):
            #    continue
            values = elem.CalculateOnIntegrationPoints(
                km.INTERNAL_VARIABLES, self.model_part.ProcessInfo
            )
            for value_ip in values:
                data_list.append(value_ip[0])
        bases.write_field_to_hdf5(
            filename, group, self.config["rvalue_name"], timestep, data_list
        )

    ###########################################################
    ###########################################################

    def ExecuteInitialize(self):
        # Create new file
        fname = self.config["snapshots_fname"]
        h5py.File(fname, "w").close()

    def ExecuteFinalizeSolutionStep(self):
        if not self.inelastic_flag:
            self.inelastic_flag = self.has_damaged_elements()

        if not self.inelastic_flag:
            group = "ELASTIC"
        else:
            group = "INELASTIC"
        print("[WriteSnapshotProcess] Snapshot group: ", group)

        fname = self.config["snapshots_fname"]
        self.write_strain(group, fname, self.timestep_counter)
        self.write_energy(group, fname, self.timestep_counter)
        self.write_rvalue(group, fname, self.timestep_counter)

        self.timestep_counter += 1

    def ExecuteFinalize(self):
        fname = self.config["snapshots_fname"]
        field = self.config["strain_name"]
        bases.generate_local_bases(
            Path.cwd(),
            field,
            fname,
            self.config["local_bases_fname_pattern"].format(field),
            self.config["local_sv_fname_pattern"].format(field),
        )
        field = self.config["energy_name"]
        bases.generate_local_bases(
            Path.cwd(),
            field,
            fname,
            self.config["local_bases_fname_pattern"].format(field),
            self.config["local_sv_fname_pattern"].format(field),
        )
        field = self.config["rvalue_name"]
        bases.generate_local_bases(
            Path.cwd(),
            field,
            fname,
            self.config["local_bases_fname_pattern"].format(field),
            self.config["local_sv_fname_pattern"].format(field),
        )
