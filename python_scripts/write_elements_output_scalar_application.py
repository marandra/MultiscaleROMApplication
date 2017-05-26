import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as msr
import os
import operator
import struct


def Factory(settings, Model):
    return WriteElementsOutputScalarApplication(settings["Parameters"], Model)


def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist


class WriteElementsOutputScalarApplication(km.Process):
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
        self.filename = param['filename'].GetString()
        self.vname = param['variable_name'].GetString()
        f = operator.attrgetter(self.vname)
        self.Var = f(msr)

    def write_results(self, filename):
        #with open(filename, 'w') as ofile:
        with open(filename, 'wb') as ofile:
            process_info = self.model_part.ProcessInfo
            for elem in self.model_part.Elements:
                variables = elem.GetValuesOnIntegrationPoints(self.Var, process_info)
                for v in variables:
                    #ofile.write(" {:18.16f}\n".format(v[0]))
                    ofile.write(struct.pack('f', v[0])) #  'f'=float32
            ofile.write(b'\n')

    def ExecuteInitialize(self):
        pass

    def ExecuteInitializeSolutionStep(self):
        #self.timestep = "-{:.3f}".format(self.model_part.ProcessInfo[km.TIME])
        #try:
        #    os.remove(self.filename + self.timestep)
        #except OSError:
        #    pass
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        #self.write_results(self.filename + self.timestep)
        pass

    def ExecuteFinalize(self):
        self.write_results(self.filename)
        pass
