import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as msr
#import bisect
import os
import operator

def Factory(settings, Model):
    return WriteElementsOutputHomogenizedVector(settings["Parameters"], Model)

def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist

def homogenization_function(self):

    var_ref = self.model_part.Elements[1].GetValuesOnIntegrationPoints(self.Var,self.model_part.ProcessInfo)
    homog_comp = var_ref[0].__len__()

    var_acum =[0.0]*homog_comp
    volume = 0.0

    for elem in self.model_part.Elements:

        var_elem = elem.GetValuesOnIntegrationPoints(self.Var,self.model_part.ProcessInfo)
        weights= elem.GetValuesOnIntegrationPoints(msr.GAUSS_WEIGHTS,self.model_part.ProcessInfo)

        for iVar in range(weights.__len__()):
            for jVar in range(homog_comp):

                var_acum[jVar] = var_acum[jVar] + var_elem[iVar][jVar]*weights[iVar][0]

            volume += weights[iVar][0]

    for iComp in range(homog_comp):
        var_acum[iComp] /= volume

    #print(volume)
    #print(var_acum)
    return var_acum

class WriteElementsOutputHomogenizedVector(km.Process):
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
        self.filename = param['filename'].GetString()
        self.vname = param['variable_name'].GetString()
        f = operator.attrgetter(self.vname)
        self.Var = f(km)

    def write_results(self, filename):
        with open(filename, 'a') as ofile:
            homog_value = homogenization_function(self)
            for v in homog_value:
                ofile.write("{:.17f}   ".format(v))
            ofile.write("\n")

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass

    def ExecuteInitializeSolutionStep(self):
        #self.timestep = "-{:.3f}".format(self.model_part.ProcessInfo[km.TIME])
        #try:
        #    os.remove(self.filename + self.timestep)
        #except OSError:
        #    pass
        self.write_results(self.filename)

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        #t = self.Model.ProcessInfo[km.TIME]
        #if t == self.Model.ProcessInfo[km.END_TIME] or self.__check_write_freq(t):
        #self.write_results(self.filename + self.timestep)
        pass


    def ExecuteFinalize(self):
        pass