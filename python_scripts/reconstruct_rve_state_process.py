# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication
import os


def Factory(settings, Model):
    return ReconstructRVEStateProcess(settings["Parameters"], Model)

def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist


class ReconstructRVEStateProcess(km.Process):
    def __init__(self, params, model):

        default_settings = km.Parameters('''{
            "model_part_name": "unset model part name",
            "filename": "specify_filename",
            "element": -1,
            "point": 0,
            "variable_1": "UNSET_VARIABLE"
        }''')
        params.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[params["model_part_name"].GetString()]
        self.filename = params['filename'].GetString()
        self.vname1 = params['variable_1'].GetString()
        self.var1 = km.KratosGlobals.GetVariable(params['variable_1'].GetString())
        self.elem = self.model_part.Elements[params['element'].GetInt()]
        self.gp = params['point'].GetInt()

    def write_results(self):
        with open(self.filename, 'a') as ofile:
            process_info = self.model_part.ProcessInfo
            CX = self.elem.GetValuesOnIntegrationPoints(self.var1, process_info)[self.gp]
            print("DEBUG PROCESS - Reconstruct RVE State:")
            print(CX)
            #ofile.write(CX)
            #ofile.write("{: .3e} {: .3e} {: .3e}  {: .3e} {: .3e} {: .3e}"  #"  {: .3e} {: .3e} {: .3e} {: .3e}"
            #    .format(
            #        CX[0], CX[3], CX[1],  # CX[0],
            #        CY[0], CY[3], CY[1],  # CY[0],
            #        ))
            ofile.write("{}\n".format(" ".join(map(str, CX))))
    
    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        with open(self.filename, 'a') as ofile:
            ofile.write("# RVE {} {} Integration point: {}\n".format(self.vname1, self.elem, self.gp))
            ofile.write("#\n")

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        #t = self.Model.ProcessInfo[km.TIME]
        #if t == self.Model.ProcessInfo[km.END_TIME] or self.__check_write_freq(t):
        self.write_results()

    def ExecuteFinalize(self):
        pass
