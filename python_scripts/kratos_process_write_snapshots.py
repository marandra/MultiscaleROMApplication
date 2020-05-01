from pathlib import Path
import h5py
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from offline_common import Common
import offline_bases


def Factory(settings, model):
    return WriteSnapshots(settings["Parameters"], model)


class WriteSnapshots(KratosMultiphysics.Process):
    def __init__(self, settings, model):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters(
            """
        {
            "model_part_name": "unset_model_part_name",
            "filename": "snapshots.hdf5",
            "svd": true,
            "config_path": "../../configuration.json"
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.filename = settings["filename"].GetString()
        self.svd = settings["svd"].GetBool()
        self.config = settings["config_path"].GetString()
        self.timestep_counter = 1
        self.inelastic_flag = False

    def has_damaged_elements(self):
        for elem in self.model_part.Elements:
            flag = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.DAMAGE_VARIABLE, self.model_part.ProcessInfo
            )
            if True in [x > 0.0 for x in flag]:
                return True

            flag = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.StructuralMechanicsApplication.ACCUMULATED_PLASTIC_STRAIN,
                self.model_part.ProcessInfo,
            )
            if True in [x > 0.0 for x in flag]:
                return True

    def write_strain(self, group, filename, timestep):
        data_list = []
        strain_macro = self.model_part.ProcessInfo[KratosMultiphysics.INITIAL_STRAIN]
        for elem in self.model_part.Elements:
            strain_vectors = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR,
                self.model_part.ProcessInfo,
            )
            for strain_ip in strain_vectors:
                for i, strain_i in enumerate(strain_ip):
                    strain_fluctuant_i = strain_i - strain_macro[i]
                    data_list.append(strain_fluctuant_i)
        # offline_bases.write_strain_to_hdf5(group, filename, timestep, data_list)
        offline_bases.write_field_to_hdf5(
            filename, group, "STRAIN_FLUCTUANT", timestep, data_list
        )

    def write_energy(self, group, filename, timestep):
        data_list = []
        for elem in self.model_part.Elements:
            strain_energy_values = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.STRAIN_ENERGY, self.model_part.ProcessInfo
            )
            for strain_energy_ip in strain_energy_values:
                data_list.append(strain_energy_ip)
        # offline_bases.write_energy_to_hdf5(group, filename, timestep_counter, data_list)
        offline_bases.write_field_to_hdf5(
            filename, group, "ENERGY_FREE", timestep, data_list
        )

    def write_rvalue(self, group, filename, timestep):
        data_list = []
        for elem in self.model_part.Elements:
            values = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.INTERNAL_VARIABLES, self.model_part.ProcessInfo
            )
            for value_ip in values:
                data_list.append(value_ip[0])
        # offline_bases.write_rvalue_to_hdf5(group, filename, timestep_counter, data_list)
        offline_bases.write_field_to_hdf5(
            filename, group, "R_VALUE", timestep, data_list
        )

    ###########################################################
    ###########################################################

    def ExecuteInitialize(self):
        # Create new file
        h5py.File(self.filename, "w").close()

        # Write global index array
        # with h5py.File(self.filename, "a") as f:
        #    f.create_dataset(
        #        "DATA/GLOBAL_ELEMENT_INDEX", data=data_list
        #    )

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        if not self.inelastic_flag:
            self.inelastic_flag = self.has_damaged_elements()

        if not self.inelastic_flag:
            group = "ELASTIC"
        else:
            group = "INELASTIC"
        self.write_strain(group, self.filename, self.timestep_counter)
        self.write_energy(group, self.filename, self.timestep_counter)
        self.write_rvalue(group, self.filename, self.timestep_counter)

        self.timestep_counter += 1

    def ExecuteFinalize(self):
        co = Common(self.config)
        if self.svd:
            field = "STRAIN_FLUCTUANT"
            offline_bases.generate_local_bases(
                Path.cwd(),
                field,
                co.context["snapshots_fname"],
                co.local_bases_fname.format(field),
                co.local_sv_fname.format(field),
            )
            field = "ENERGY_FREE"
            offline_bases.generate_local_bases(
                Path.cwd(),
                field,
                co.context["snapshots_fname"],
                co.local_bases_fname.format(field),
                co.local_sv_fname.format(field),
            )
            field = "R_VALUE"
            offline_bases.generate_local_bases(
                Path.cwd(),
                field,
                co.context["snapshots_fname"],
                co.local_bases_fname.format(field),
                co.local_sv_fname.format(field),
            )
