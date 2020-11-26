from pathlib import Path
import h5py
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as sm
from offline_bases import Bases


def Factory(settings, model):
    return WriteSnapshots(settings["Parameters"], model)


class WriteSnapshots(km.Process):
    def __init__(self, settings, model):
        km.Process.__init__(self)

        default_settings = km.Parameters(
            """{"model_part_name": "unset_model_part_name"}"""
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.bases = Bases(config_fname="../../configuration.json")
        self.timestep_counter = 1
        self.inelastic_flag = False

    def has_damaged_elements(self):
        for elem in self.model_part.Elements:

            #  for RVELaw
            flag = elem.CalculateOnIntegrationPoints(
                sm.DAMAGE, self.model_part.ProcessInfo
            )
            if True in [x > 0.0 for x in flag]:
                return True

            #  for Damage CLs
            flag = elem.CalculateOnIntegrationPoints(
                km.DAMAGE_VARIABLE, self.model_part.ProcessInfo
            )
            if True in [x > 0.0 for x in flag]:
                return True

            #  for Plasticity CLs
            flag = elem.CalculateOnIntegrationPoints(
                sm.ACCUMULATED_PLASTIC_STRAIN, self.model_part.ProcessInfo,
            )
            if True in [x > 0.0 for x in flag]:
                return True

    def write_strain(self, group, filename, timestep):
        data_list = []
        strain_macro = self.model_part.ProcessInfo[km.INITIAL_STRAIN]
        for elem in self.model_part.Elements:
            strain_vectors = elem.CalculateOnIntegrationPoints(
                km.GREEN_LAGRANGE_STRAIN_VECTOR, self.model_part.ProcessInfo,
            )
            for strain_ip in strain_vectors:
                for i, strain_i in enumerate(strain_ip):
                    strain_fluctuant_i = strain_i - strain_macro[i]
                    data_list.append(strain_fluctuant_i)
        self.bases.write_field_to_hdf5(
            filename, group, "STRAIN_FLUCTUANT", timestep, data_list
        )

    def write_energy(self, group, filename, timestep):
        data_list = []
        for elem in self.model_part.Elements:
            strain_energy_values = elem.CalculateOnIntegrationPoints(
                km.STRAIN_ENERGY, self.model_part.ProcessInfo
            )
            for strain_energy_ip in strain_energy_values:
                data_list.append(strain_energy_ip)
        self.bases.write_field_to_hdf5(
            filename, group, "ENERGY_FREE", timestep, data_list
        )

    def write_rvalue(self, group, filename, timestep):
        data_list = []
        for elem in self.model_part.Elements:
            values = elem.CalculateOnIntegrationPoints(
                km.INTERNAL_VARIABLES, self.model_part.ProcessInfo
            )
            for value_ip in values:
                data_list.append(value_ip[0])
        self.bases.write_field_to_hdf5(filename, group, "R_VALUE", timestep, data_list)

    ###########################################################
    ###########################################################

    def ExecuteInitialize(self):
        # Create new file
        fname = self.bases.context["snapshots_fname"]
        h5py.File(fname, "w").close()

    def ExecuteFinalizeSolutionStep(self):
        if not self.inelastic_flag:
            self.inelastic_flag = self.has_damaged_elements()

        if not self.inelastic_flag:
            group = "ELASTIC"
        else:
            group = "INELASTIC"

        fname = self.bases.context["snapshots_fname"]
        self.write_strain(group, fname, self.timestep_counter)
        self.write_energy(group, fname, self.timestep_counter)
        self.write_rvalue(group, fname, self.timestep_counter)

        self.timestep_counter += 1

    def ExecuteFinalize(self):
        fname = self.bases.context["snapshots_fname"]
        field = "STRAIN_FLUCTUANT"
        self.bases.generate_local_bases(
            Path.cwd(),
            field,
            fname,
            self.bases.context["local_bases_fname_pattern"].format(field),
            self.bases.context["local_sv_fname_pattern"].format(field),
        )
        field = "ENERGY_FREE"
        self.bases.generate_local_bases(
            Path.cwd(),
            field,
            fname,
            self.bases.context["local_bases_fname_pattern"].format(field),
            self.bases.context["local_sv_fname_pattern"].format(field),
        )
        field = "R_VALUE"
        self.bases.generate_local_bases(
            Path.cwd(),
            field,
            fname,
            self.bases.context["local_bases_fname_pattern"].format(field),
            self.bases.context["local_sv_fname_pattern"].format(field),
        )
