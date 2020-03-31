import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import h5py
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
            "svd": true
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.filename = settings["filename"].GetString()
        self.svd = settings["svd"].GetBool()

    def has_damaged_elements(self):
        for elem in self.model_part.Elements:
            flag = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.DAMAGE_VARIABLE, self.model_part.ProcessInfo
            )
            if True in [x > 0.0 for x in flag]:
                return True

            flag = elem.GetValuesOnIntegrationPoints(
                KratosMultiphysics.StructuralMechanicsApplication.ACCUMULATED_PLASTIC_STRAIN,
                self.model_part.ProcessInfo,
            )
            if True in [x > 0.0 for y in flag for x in y]:
                return True

    def write_strain(self, group):
        data_list = []
        strain_macro = self.model_part.ProcessInfo[KratosMultiphysics.INITIAL_STRAIN]
        for elem in self.model_part.Elements:
            strain_vectors = elem.GetValuesOnIntegrationPoints(
                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR,
                self.model_part.ProcessInfo,
            )
            for strain_ip in strain_vectors:
                for i, strain_i in enumerate(strain_ip):
                    strain_fluctuant_i = strain_i - strain_macro[i]
                    data_list.append(strain_fluctuant_i)
        with h5py.File(self.filename, "a") as f:
            f.create_dataset(
                "{}/STRAIN_FLUCTUANT/{}".format(group, self.timestep_counter),
                data=data_list,
            )

    def write_energy(self, group):
        data_list = []
        for elem in self.model_part.Elements:
            strain_energy_values = elem.GetValuesOnIntegrationPoints(
                KratosMultiphysics.STRAIN_ENERGY, self.model_part.ProcessInfo
            )
            for strain_energy_ip in strain_energy_values:
                data_list.append(strain_energy_ip[0])
        with h5py.File(self.filename, "a") as f:
            f.create_dataset(
                "{}/ENERGY_FREE/{}".format(group, self.timestep_counter), data=data_list
            )

    def write_rvalue(self, group):
        data_list = []
        for elem in self.model_part.Elements:
            values = elem.GetValuesOnIntegrationPoints(
                KratosMultiphysics.INTERNAL_VARIABLES, self.model_part.ProcessInfo
            )
            for value_ip in values:
                data_list.append(value_ip[0])
        with h5py.File(self.filename, "a") as f:
            f.create_dataset(
                "{}/R_VALUE/{}".format(group, self.timestep_counter), data=data_list
            )

    ###########################################################
    ###########################################################

    def ExecuteInitialize(self):
        self.timestep_counter = 1
        self.inelastic_flag = False
        # Create new file
        h5py.File(self.filename, "w").close()

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
            self.write_strain(group)
            self.write_energy(group)
            self.write_rvalue(group)
        else:
            group = "INELASTIC"
            self.write_strain(group)
            self.write_energy(group)
            self.write_rvalue(group)

        self.timestep_counter += 1

    def ExecuteFinalize(self):
        if self.svd:
            field = "STRAIN_FLUCTUANT"
            offline_bases.generate_local_bases(".", field)
            field = "ENERGY_FREE"
            offline_bases.generate_local_bases(".", field)
            field = "R_VALUE"
            offline_bases.generate_local_bases(".", field)
