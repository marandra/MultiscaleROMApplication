import KratosMultiphysics as Kratos
import KratosMultiphysics.MultiscaleROMApplication as MultiscaleROM
import os
import KratosMultiphysics.MultiscaleROMApplication.io_utilities as io_utilities


def Factory(settings, model):
    return WriteRveReconstructionData(settings["Parameters"], model)


def append_to_json(filename, new_data):
   data = io_utilities.read_json(filename)
   data["interpolation_parameters"].append(new_data["interpolation_parameters"])
   data["macro_strain"].append(new_data["macro_strain"])
   data["stress"].append(new_data["stress"])
   io_utilities.write_json(filename, data)

class WriteRveReconstructionData(Kratos.Process):
    def __init__(self, settings, model):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters("""
        {
            "model_part_name": "unset_model_part_name",
            "filename": "unset_filename",
            "element": 1,
            "integration_point": 0
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()
        self.element = settings['element'].GetInt()
        self.ip = settings['integration_point'].GetInt()

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        self.data = {}
        self.data["interpolation_parameters"] = []
        self.data["macro_strain"] = []
        self.data["stress"] = []
        io_utilities.write_json(self.filename, self.data)

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        for elem in self.model_part.Elements:
            if elem.Id == self.element:
                # Get fluctuant displacement
                ip_data = elem.GetValuesOnIntegrationPoints(
                    MultiscaleROM.REDUCED_MODES_WEIGHTS, self.model_part.ProcessInfo)
                self.data["interpolation_parameters"] = ip_data[self.ip]
                # Get macro strain
                ip_data = elem.GetValuesOnIntegrationPoints(
                    Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, self.model_part.ProcessInfo)
                self.data["macro_strain"] = ip_data[self.ip]
                # Get stress vector list
                ip_data = elem.GetValuesOnIntegrationPoints(
                    MultiscaleROM.CAUCHY_STRESS_VECTOR_LIST, self.model_part.ProcessInfo)
                self.data["stress"] = ip_data[self.ip]
        append_to_json(self.filename, self.data)

    def ExecuteFinalize(self):
        pass
