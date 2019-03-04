import KratosMultiphysics as Kratos
import os


def Factory(settings, Model):
    return WriteIntegrationWeights(settings["Parameters"], Model)


class WriteIntegrationWeights(Kratos.Process):
    def __init__(self, settings, model):
        Kratos.Process.__init__(self)
        default_settings = Kratos.Parameters("""
        {
            "model_part_name": "unset_model_part_name",
            "filename": "integration_weights"
        }
        """)

        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()
        self.first_time_flag = True

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass

    def ExecuteInitializeSolutionStep(self):
        if self.first_time_flag:
            with open(self.filename, 'w') as ofile:
                for elem in self.model_part.Elements:
                    ip_weights = elem.GetValuesOnIntegrationPoints(
                        Kratos.INTEGRATION_WEIGHT, self.model_part.ProcessInfo)
                    for ip_weight in ip_weights:
                        ofile.write("{}\n".format(ip_weight[0]))
            self.first_time_flag = False

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteFinalize(self):
        pass
