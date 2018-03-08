import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as msr
import os
import operator
import math

def Factory(settings, Model):
    return WriteElementsHomogenizedOutput(settings["Parameters"], Model)

def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist

def homogenization_function(self):
    var_ref = self.model_part.Elements[1].GetValuesOnIntegrationPoints(self.Var,self.model_part.ProcessInfo)
    nr_comp = len(var_ref[0])
    var_accum = [0.0] * nr_comp
    volume = 0.0

    for e, elem in enumerate(self.model_part.Elements):
        values = elem.GetValuesOnIntegrationPoints(self.Var, self.model_part.ProcessInfo)
        weights = elem.GetValuesOnIntegrationPoints(km.INTEGRATION_WEIGHT, self.model_part.ProcessInfo)
        weights = [x[0] for x in weights] # to unpack received list-inside-list
        for i, w in enumerate(weights):
            # used in HPROM case, to ignore GP
            if w == -1:
                continue
            for j in range(nr_comp):
                var_accum[j] += values[i][j] * w
            volume += w
    for i in range(nr_comp):
        var_accum[i] /= volume
    return var_accum

def compute_vonmisses_stress(hs):
    s = (hs[0] + hs[1] + hs[2]) / 3
    d = [hs[0] - s, hs[1] - s, hs[2] - s, hs[3], hs[4], hs[5]]
    dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + \
         d[3] * d[3] + d[4] * d[4] + d[5] * d[5] 
    vm = math.sqrt(dd)
    return vm

class WriteElementsHomogenizedOutput(km.Process):
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
        self.filename = param['filename'].GetString()
        self.vname = param['variable_name'].GetString()
        f = operator.attrgetter(self.vname)
        self.Var = f(km)

    def write_results(self, filename):
        with open(filename, 'a') as ofile:
            homog_value = homogenization_function(self)
            #von_misses = compute_vonmisses_stress(homog_value)
            for v in homog_value:
                ofile.write("{:.17f}   ".format(v))
            #ofile.write("{:.17f}".format(von_misses))
            ofile.write("\n")

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
        self.write_results(self.filename)

    def ExecuteFinalize(self):
        pass
