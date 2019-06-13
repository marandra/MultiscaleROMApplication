import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as structural
import os
import operator


def Factory(settings, Model):
    return WriteElementsOutputScalar(settings["Parameters"], Model)


def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist


class WriteElementsOutputScalar(km.Process):
    def __init__(self, settings, Model):

        km.Process.__init__(self)

        default_settings = km.Parameters("""
        {
            "model_part_name": "unset_model_part_name",
            "filename": "unset_filename"
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()

    def write_results(self):
        with open(self.filename, 'w') as ofile:
            ofile.write("{}\n".format(self.timestep_counter))

    def ExecuteInitialize(self):
        self.timestep_counter = 0
        self.inelastic_flag = False
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
        if not self.inelastic_flag:
            for elem in self.model_part.Elements:
                flag = elem.CalculateOnIntegrationPoints(km.DAMAGE_VARIABLE, self.model_part.ProcessInfo)
                if True in [x > 0.0 for x in flag]:
                    self.inelastic_flag = True
                    self.write_results()
                    break
                flag = elem.GetValuesOnIntegrationPoints(structural.ACCUMULATED_PLASTIC_STRAIN, self.model_part.ProcessInfo)
                if True in [x for y in flag for x in y]:
                    print("DEBUG break in case a function enters here. We need to test it with plasticity CL")
                    error()
                    self.inelastic_flag = True
                    self.write_results()
                    break
        self.timestep_counter = self.timestep_counter + 1

    def ExecuteFinalize(self):
        pass
