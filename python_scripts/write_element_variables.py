import KratosMultiphysics as km
import KratosMultiphysics.MultiScaleApplication as mss # <- check is used
import os
import operator

def Factory(settings, Model):
    #if(type(settings) != km.Parameters):
    #    raise Exception("expected input is Parameters object, encapsulating a json string")
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
    	f = operator.attrgetter(param['variable_x'].GetString())
    	self.VarX = f(km)
    	f = operator.attrgetter(param['variable_y'].GetString())
    	self.VarY = f(km)
    	self.elem = self.model_part.Elements[param['element'].GetInt()]
    	self.gp = param['gauss_point'].GetInt()
    	self.compx = param['component_x'].GetInt()
    	self.compy = param['component_y'].GetInt()
    
    def write_results(self):
    	with open(self.filename, 'a') as ofile:
            process_info = self.model_part.ProcessInfo
            varx = self.elem.GetValuesOnIntegrationPoints(self.VarX, process_info)[self.gp]
            vary = self.elem.GetValuesOnIntegrationPoints(self.VarY, process_info)[self.gp]
    	    ofile.write("{} {}\n".format(varx[self.compx], vary[self.compy]))
    
    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
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
