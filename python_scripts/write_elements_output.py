import KratosMultiphysics as km
import os
import operator
import struct


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
            "write_frequency": "every_timestep",
            "write_mode": "ascii",
            "variable_location": "core",
            "variable_name": "unset_variable_name"
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings['model_part_name'].GetString()]
        self.filename = settings['filename'].GetString()
        self.write_frequency = settings['write_frequency'].GetString()
        self.write_mode = settings['write_mode'].GetString()
        self.var_location = settings['variable_location'].GetString()
        self.var_name = settings['variable_name'].GetString()
        f = operator.attrgetter(self.var_name)
        if self.var_location == "core":
            self.var = f(km)
        else:
            f = operator.attrgetter(self.var_location)
            app = f(km)
            f = operator.attrgetter(self.var_name)
            self.var = f(app)

    def write_results(self, filename):

        def write_results_binary():
            with open(filename, 'wb') as ofile:
                process_info = self.model_part.ProcessInfo
                for elem in self.model_part.Elements:
                    variables = elem.GetValuesOnIntegrationPoints(self.var, process_info)
                    for v in variables:
                        for comp in v:
                            ofile.write(struct.pack('f', comp)) # 'f'=float32
                ofile.write(b'\n')

        def write_results_ascii():
            with open(filename, 'w') as ofile:
                process_info = self.model_part.ProcessInfo
                for elem in self.model_part.Elements:
                    variables = elem.GetValuesOnIntegrationPoints(self.var, process_info)
                    for v in variables:
                        for comp in v:
                            ofile.write("{:18.16f}\n".format(comp))

        if self.write_mode == "binary":
            write_results_binary()
        else:
            write_results_ascii()

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass

    def ExecuteInitializeSolutionStep(self):
        self.timestep = "-{:.3f}".format(self.model_part.ProcessInfo[km.TIME])
        try:
            os.remove(self.filename + self.timestep)
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
            self.write_results(self.filename + self.timestep)

    def ExecuteFinalize(self):
        if self.write_frequency == "last_timestep":
            self.write_results(self.filename)
