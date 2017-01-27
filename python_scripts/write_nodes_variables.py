import KratosMultiphysics as km
import KratosMultiphysics.MultiScaleApplication as mss # <- check is used
import os
import operator

def Factory(settings, Model):
    return WriteNodesVariables(settings["Parameters"], Model)


def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist


class WriteNodesVariables(km.Process):
    def __init__(self, param, Model):
        self.Model = Model[param['model_part_name'].GetString()]
        self.BaseName = None
        self.filename = param['filename'].GetString()
        self.vname1 = param['variable_1'].GetString()
        f = operator.attrgetter(self.vname1)
        self.Var1 = f(km)
        self.vname2 = param['variable_2'].GetString()
        f = operator.attrgetter(self.vname2)
        self.Var2 = f(km)
        self.vname3 = param['variable_3'].GetString()
        f = operator.attrgetter(self.vname3)
        self.Var3 = f(km)
        self.vname4 = param['variable_4'].GetString()
        f = operator.attrgetter(self.vname4)
        self.Var4 = f(km)
        self.node = param['node'].GetInt()
    
    def write_results(self):
        with open(self.filename, 'a') as ofile:
                var1 = self.Model.Nodes[self.node].GetSolutionStepValue(self.Var1)
                ofile.write("  {}".format(var1))
                var2 = self.Model.Nodes[self.node].GetSolutionStepValue(self.Var2)
                ofile.write("  {}".format(var2))
                var3 = self.Model.Nodes[self.node].GetSolutionStepValue(self.Var3)
                ofile.write("  {}".format(var3))
                var4 = self.Model.Nodes[self.node].GetSolutionStepValue(self.Var4)
                ofile.write("  {}".format(var4))
                ofile.write("\n".format())
    
    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        with open(self.filename, 'a') as ofile:
            ofile.write("#  node {}\n#".format(self.node))
            ofile.write("  {}".format(self.vname1))
            ofile.write("  {}".format(self.vname2))
            ofile.write("  {}".format(self.vname3))
            ofile.write("  {}".format(self.vname4))
            ofile.write("\n".format())
    
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
