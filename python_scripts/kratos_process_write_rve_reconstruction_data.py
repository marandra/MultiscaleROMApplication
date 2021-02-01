import KratosMultiphysics as Kratos
import KratosMultiphysics.MultiscaleROMApplication as MultiscaleROM
import os
import json


def write_json(filename, data_dict):
    with open(filename, "w") as fo:
        json.dump(data_dict, fo, indent=2)


def read_json(filename):
    with open(filename) as f:
        data_dict = json.load(f)
    return data_dict


def Factory(settings, model):
    return WriteRveReconstructionData(settings["Parameters"], model)


def append_to_json(filename, new_data):
    data = read_json(filename)
    data["interpolation_parameters"].append(new_data["interpolation_parameters"])
    data["macro_strain"].append(new_data["macro_strain"])
    data["strain_energy"].append(new_data["strain_energy"])
    data["r_value"].append(new_data["r_value"])
    write_json(filename, data)


class WriteRveReconstructionData(Kratos.Process):
    def __init__(self, settings, model):
        Kratos.Process.__init__(self)

        default_settings = Kratos.Parameters(
            """
        {
            "model_part_name": "unset_model_part_name",
            "filename": "unset_filename",
            "element": 1,
            "integration_point": 0
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[settings["model_part_name"].GetString()]
        self.filename = settings["filename"].GetString()
        self.element = settings["element"].GetInt()
        self.ip = settings["integration_point"].GetInt()

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        self.data = {}
        self.data["interpolation_parameters"] = []
        self.data["macro_strain"] = []
        self.data["strain_energy"] = []
        self.data["r_value"] = []
        write_json(self.filename, self.data)

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
                ip_data = elem.CalculateOnIntegrationPoints(
                    MultiscaleROM.REDUCED_MODES_WEIGHTS, self.model_part.ProcessInfo
                )
                data = [x for x in ip_data[self.ip]]
                self.data["interpolation_parameters"] = data

                # Get macro strain
                ip_data = elem.CalculateOnIntegrationPoints(
                    #Kratos.GREEN_LAGRANGE_STRAIN_VECTOR, self.model_part.ProcessInfo
                    Kratos.STRAIN, self.model_part.ProcessInfo
                )
                data = [x for x in ip_data[self.ip]]
                self.data["macro_strain"] = data

                # Get strain energy list
                ip_data = elem.CalculateOnIntegrationPoints(
                    MultiscaleROM.STRAIN_ENERGY_VECTOR, self.model_part.ProcessInfo
                )
                tmp_Vector = ip_data[self.ip]
                tmp_list = []
                for i in tmp_Vector:
                    tmp_list.append(i)
                self.data["strain_energy"] = tmp_list

                # Get r_value list
                ip_data = elem.CalculateOnIntegrationPoints(
                    Kratos.INTERNAL_VARIABLES, self.model_part.ProcessInfo
                )
                data = [x for x in ip_data[self.ip]]
                self.data["r_value"] = data

        append_to_json(self.filename, self.data)

    def ExecuteFinalize(self):
        pass
