import KratosMultiphysics as Kratos
import KratosMultiphysics.MultiscaleROMApplication as MultiscaleROM
import os
import struct


def Factory(settings, model):
    return WriteRveInterpolationParameters(settings["Parameters"], model)


class WriteRveInterpolationParameters(Kratos.Process):
    def __init__(self, settings, model):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
        {
            "model_part_name": "unset_model_part_name",
            "filename": "unset_filename",
            "write_mode": "ascii",
            "element": 1,
            "integration_point": 0
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()
        self.write_mode = settings['write_mode'].GetString()
        self.element = settings['element'].GetInt()
        self.ip = settings['integration_point'].GetInt()

    def append_results(self, filename):
        def write_value_binary(of, value):
            of.write(struct.pack('f', value))  # 'f'=float32

        def write_value_ascii(of, value):
            of.write("{:18.16f} ".format(value))

        if self.write_mode == "binary":
            write_value = write_value_binary
            ofile = open(filename, 'ab')
        else:
            write_value = write_value_ascii
            ofile = open(filename, 'a')

        for elem in self.model_part.Elements:
            if elem.Id == self.element:
                interpolation_params = elem.GetValuesOnIntegrationPoints(
                    MultiscaleROM.REDUCED_MODES_WEIGHTS, self.model_part.ProcessInfo)
                for v in interpolation_params[self.ip]:
                    write_value(ofile, v)

        if self.write_mode == "binary":
            ofile.write(b'\n')
        else:
            ofile.write("\n")
        ofile.close()

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        self.append_results(self.filename)

    def ExecuteFinalize(self):
        pass
