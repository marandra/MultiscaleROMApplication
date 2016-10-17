import KratosMultiphysics as km
import KratosMultiphysics.MultiScaleApplication as mss # <- check is used
import os
import operator

def Factory(settings, Model):
    return WriteElementVariables(settings["Parameters"], Model)


def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist


class WriteElementVariables(km.Process):
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
    	self.filename = param['filename'].GetString()
        self.vname1 = param['variable_1'].GetString()
        f = operator.attrgetter(self.vname1)
        self.Var1 = f(km)
        self.vname2 = param['variable_2'].GetString()
        f = operator.attrgetter(self.vname2)
        self.Var2 = f(km)
    	self.elem = self.model_part.Elements[param['element'].GetInt()]
    	self.gp = param['gauss_point'].GetInt()
    
    def write_results(self):
    	with open(self.filename, 'a') as ofile:
            process_info = self.model_part.ProcessInfo
            CX = self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp]
            CY = self.elem.GetValuesOnIntegrationPoints(self.Var2, process_info)[self.gp]
            ofile.write("{: .3e} {: .3e} {: .3e}  {: .3e} {: .3e} {: .3e}" #"  {: .3e} {: .3e} {: .3e} {: .3e}"
                .format(
                CX[0], CX[3], CX[1], # CX[0],
                CY[0], CY[3], CY[1], # CY[0],
                ))
            ofile.write("\n")
    
    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        with open(self.filename, 'a') as ofile:
            ofile.write("#{:<32}  {:<32}\n".format(self.vname1, self.vname2))
            ofile.write("#{:<10} {:<10} {:<10}" #" {:<10} {:<10}"
                        "  {:<10} {:<10} {:<10}" #" {:<10} {:<10}"
                .format(
                "Comp XX", "Comp YY", "Comp XY", #"Comp 4", "Comp 5",
                "Comp XX", "Comp YY", "Comp XY", #"Comp 4", "Comp 5"
                ))
            ofile.write("\n")
                                                                            

        #self.write_results() 
    	#self.Tn = self.Model.ProcessInfo[km.TIME]
    
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
