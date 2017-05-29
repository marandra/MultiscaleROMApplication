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
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
        self.filename = param['filename'].GetString()
        self.write_frequency = param['write_frequency'].GetString()
        self.write_mode = param['write_mode'].GetString()
        self.var_reach = param['variable_reach'].GetString()
        self.var_name = param['variable_name'].GetString()
        f = operator.attrgetter(self.var_name)
        if self.var_reach == "core":
            self.Var = f(km)
        else:
            import KratosMultiphysics.MultiscaleROMApplication as msr
            self.Var = f(msr)

    def write_results(self, filename):

        def write_results_binary():
            with open(filename, 'wb') as ofile:
                process_info = self.model_part.ProcessInfo
                for elem in self.model_part.Elements:
                    variables = elem.GetValuesOnIntegrationPoints(self.Var, process_info)
                    for v in variables:
                        ofile.write(struct.pack('f', v[0])) # 'f'=float32
                ofile.write(b'\n')

        def write_results_ascii():
            with open(filename, 'w') as ofile:
                process_info = self.model_part.ProcessInfo
                for elem in self.model_part.Elements:
                    variables = elem.GetValuesOnIntegrationPoints(self.Var, process_info)
                    for v in variables:
                        ofile.write("{:18.16f}\n".format(v[0]))

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
