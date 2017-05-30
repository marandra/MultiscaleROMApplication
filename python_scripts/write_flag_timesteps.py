import KratosMultiphysics as km
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

        default_settings = km.Parameters("""
        {
            "mesh_id": 0,
            "model_part_name": "unset_model_part_name",
            "filename": "unset_filename",
            "flag_location": "core",
            "flag_name": "unset_flag_name"
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()
        self.var_location = settings['flag_location'].GetString()
        self.var_name = settings['flag_name'].GetString()
        f = operator.attrgetter(self.var_name)
        if self.var_location == "core":
            self.var = f(km)
        else:
            f = operator.attrgetter(self.var_location)
            app = f(km)
            f = operator.attrgetter(self.var_name)
            self.var = f(app)

    def write_results(self):
        with open(self.filename, 'w') as ofile:
            ofile.write("{}\n".format(self.timestep_counter))

    def ExecuteInitialize(self):
        self.timestep_counter = -1
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
                flag = elem.GetValuesOnIntegrationPoints(self.var, self.model_part.ProcessInfo)
                if 1 in [x for y in flag for x in y]:
                    self.inelastic_flag = True
                    self.write_results()
        self.timestep_counter = self.timestep_counter + 1

    def ExecuteFinalize(self):
        pass
