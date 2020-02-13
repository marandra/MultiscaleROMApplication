import KratosMultiphysics
import h5py


def Factory(settings, model):
    return WriteSnapshots(settings["Parameters"], model)


class WriteSnapshots(KratosMultiphysics.Process):
    def __init__(self, settings, model):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name": "unset_model_part_name",
            "filename": "snapshots.hdf5"
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()


    def write_strain(self, time):
        data_list = []
        strain_macro = self.model_part.ProcessInfo[KratosMultiphysics.INITIAL_STRAIN]
        for elem in self.model_part.Elements:
            strain_vectors = elem.GetValuesOnIntegrationPoints(
                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR,
                self.model_part.ProcessInfo)
            for strain_ip in strain_vectors:
                for i, strain_i in enumerate(strain_ip):
                    strain_fluctuant_i = strain_i - strain_macro[i]
                    data_list.append(strain_fluctuant_i)
        with h5py.File(self.filename, 'a') as f:
            f.create_dataset("strain/{}".format(time), data=data_list)

    def write_energy_strain(self, time):
        data_list = []
        for elem in self.model_part.Elements:
            strain_energy_values = elem.GetValuesOnIntegrationPoints(
                KratosMultiphysics.STRAIN_ENERGY, self.model_part.ProcessInfo)
            for strain_energy_ip in strain_energy_values:
                data_list.append(strain_energy_ip[0])
        with h5py.File(self.filename, 'a') as f:
           f.create_dataset("energy_strain/{}".format(time), data=data_list)


    def ExecuteInitialize(self):
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
        time = "{:.3f}".format(self.model_part.ProcessInfo[KratosMultiphysics.TIME])
        self.write_strain(time)
        self.write_energy_strain(time)

    def ExecuteFinalize(self):
        pass
