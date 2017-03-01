import KratosMultiphysics as km
import os
import operator

def Factory(settings, Model):
    return WriteNodeVariablesProcess(settings["Parameters"], Model)


def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist


class WriteNodeVariablesProcess(km.Process):
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
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
        for node in self.model_part.Nodes:
            if node.Id == param['node'].GetInt():
                self.node = node
                break

    def write_results(self):
        with open(self.filename, 'a') as ofile:
                var1 = self.node.GetSolutionStepValue(self.Var1)
                ofile.write("  {:.8e}".format(var1))
                var2 = self.node.GetSolutionStepValue(self.Var2)
                ofile.write("  {:.8e}".format(var2))
                var3 = self.node.GetSolutionStepValue(self.Var3)
                ofile.write("  {:.8e}".format(var3))
                var4 = self.node.GetSolutionStepValue(self.Var4)
                ofile.write("  {:.8e}".format(var4))
                ofile.write("\n".format())
    
    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        with open(self.filename, 'a') as ofile:
            ofile.write("#  node {}\n#".format(self.node.Id))
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
        self.write_results()
    
    def ExecuteFinalize(self):
        pass
