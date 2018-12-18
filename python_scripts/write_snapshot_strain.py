import KratosMultiphysics as Kratos
import os
import struct


def Factory(settings, model):
    return WriteSnapshotStrain(settings["Parameters"], model)


class WriteSnapshotStrain(Kratos.Process):
    def __init__(self, settings, model):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
        {
            "mesh_id": 0,
            "model_part_name": "unset_model_part_name",
            "filename": "unset_filename",
            "write_frequency": "every_timestep",
            "write_mode": "ascii"
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()
        self.write_frequency = settings['write_frequency'].GetString()
        self.write_mode = settings['write_mode'].GetString()
        self.time = "-{:.3f}".format(0.0)

    def write_results(self, filename):

        def write_value_binary(of, value):
            of.write(struct.pack('f', value))  # 'f'=float32

        def write_value_ascii(of, value):
            of.write("{:18.16f}\n".format(value))

        if self.write_mode == "binary":
            write_value = write_value_binary
            ofile = open(filename, 'wb')
        else:
            write_value = write_value_ascii
            ofile = open(filename, 'w')

        strain_macro = self.model_part.ProcessInfo[Kratos.INITIAL_STRAIN]
        for elem in self.model_part.Elements:
            strain_vectors = elem.GetValuesOnIntegrationPoints(
                Kratos.GREEN_LAGRANGE_STRAIN_VECTOR,
                self.model_part.ProcessInfo)
            for strain_ip in strain_vectors:
                for i, strain_i in enumerate(strain_ip):
                    strain_fluctuant_i = strain_i - strain_macro[i]
                    write_value(ofile, strain_fluctuant_i)

        if self.write_mode == "binary":
            ofile.write(b'\n')
        ofile.close()

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass

    def ExecuteInitializeSolutionStep(self):
        self.time = "-{:.3f}".format(self.model_part.ProcessInfo[Kratos.TIME])
        try:
            os.remove(self.filename + self.time)
        except OSError:
            pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        if self.write_frequency == "every_timestep":
            self.write_results(self.filename + self.time)

    def ExecuteFinalize(self):
        if self.write_frequency == "last_timestep":
            self.write_results(self.filename)
