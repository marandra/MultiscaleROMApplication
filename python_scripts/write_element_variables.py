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
            print("VAR")
            print(self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp])
            CompX0 = self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp][0]
            CompX1 = self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp][1]
            CompX2 = self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp][2]
            CompX3 = self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp][3]
            #CompX4 = self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp][4]
            #CompX5 = self.elem.GetValuesOnIntegrationPoints(self.Var1, process_info)[self.gp][5]
            CompY0 = self.elem.GetValuesOnIntegrationPoints(self.Var2, process_info)[self.gp][0]
            CompY1 = self.elem.GetValuesOnIntegrationPoints(self.Var2, process_info)[self.gp][1]
            CompY2 = self.elem.GetValuesOnIntegrationPoints(self.Var2, process_info)[self.gp][2]
            CompY3 = self.elem.GetValuesOnIntegrationPoints(self.Var2, process_info)[self.gp][3]
            #CompY4 = self.elem.GetValuesOnIntegrationPoints(self.Var2, process_info)[self.gp][4]
            #CompY5 = self.elem.GetValuesOnIntegrationPoints(self.Var2, process_info)[self.gp][5]
            ofile.write("{: .3e} {: .3e} {: .3e} {: .3e}  {: .3e} {: .3e} {: .3e} {: .3e}"# "  {: .3e} {: .3e} {: .3e} {: .3e}"
                .format(
                CompX0, CompX1, CompX2, CompX3, #CompX4, CompX5,
                CompY0, CompY1, CompY2, CompY3, #CompY4, CompY5
                ))
            ofile.write("\n")
    
    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        with open(self.filename, 'a') as ofile:
            ofile.write("#{:<43}  {:<43}\n".format(self.vname1, self.vname2))
            ofile.write("#{:<10} {:<10} {:<10} {:<10}"
                        #" {:<10} {:<10}"
                        "  {:<10} {:<10} {:<10} {:<10}"
                        #"  {:<10}  {:<10}"
                .format(
                "Comp 0", "Comp 1", "Comp 2", "Comp 3", #"Comp 4", "Comp 5",
                "Comp 0", "Comp 1", "Comp 2", "Comp 3", #"Comp 4", "Comp 5"
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
