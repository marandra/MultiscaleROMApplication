import KratosMultiphysics as km
import KratosMultiphysics.MultiScaleApplication as mss # <- check is used
import os

def Factory(settings, Model):
    if(type(settings) != km.Parameters):
        raise Exception("expected input is Parameters object, encapsulating a json string")
    return WriteGNUNodalMultiple(settings["Parameters"], Model)


def parameters_get_list_int(settings_list):
    olist = []
    for i in range(settings_list.size()):
        olist.append(settings_list[i].GetInt())
    return olist


class WriteGNUNodalMultiple(km.Process):
    def __init__(self, param, Model):
        self.Model = Model[param['modelpart_name'].GetString()]
    	self.BaseName = None
    	self.Name = param['filename'].GetString()
    	self.VarX = km.DISPLACEMENT_Y
    	#comm = "self.VarX = {}".format(param['variable_x'].GetString())
        #exec(comm)
    	self.VarY = km.REACTION_Y
    	#comm = "self.VarY = {}".format(param['variable_y'].GetString)
        #exec(comm)
    	self.NodesX = parameters_get_list_int(param['nodes_x'])
    	self.NodesY = parameters_get_list_int(param['nodes_y'])
    	self.FactorX = 1.
    	self.FactorY = 1.
    	self.XSumFactor = 1.
    	self.YSumFactor = 1.
    	self.Frequency = None
    	self.ofile = None
    	self.Tn = None
    
    def __check_write_freq(self,t):
    	r = True
    	f = self.Frequency
    	if(f is not None):
            if(self.Tn is None):
            	self.Tn = t
            else:
            	dt = t-self.Tn
            	if(dt >= f):
                    self.Tn = t
            	else:
            	    r = False
            	return r
            
    def write_results(self):
    	with open(self.Name, 'w+') as self.ofile:
    	    sum_x = 0.0
    	    sum_y = 0.0
    	    x_sum_fac = self.XSumFactor
    	    y_sum_fac = self.YSumFactor
    	    for i in self.NodesX:
    	    	sum_x += x_sum_fac*(self.Model.Nodes[i].GetSolutionStepValue(self.VarX))
    	    for i in self.NodesY:
    	    	sum_y += y_sum_fac*(self.Model.Nodes[i].GetSolutionStepValue(self.VarY))
    	    self.ofile.write("{}  {}\n".format(self.FactorX*sum_x, self.FactorY*sum_y))
    
    def ExecuteInitialize(self):
        try:
            os.remove(self.Name)
        except OSError:
            pass
        self.write_results() 
    	self.Tn = self.Model.ProcessInfo[km.TIME]
    
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
