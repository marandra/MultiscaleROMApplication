import KratosMultiphysics as km
import bisect
import os

def Factory(settings, Model):
    if(type(settings) != km.Parameters):
        raise Exception("expected input is Parameters object, encapsulating a json string")
    return WriteGeneralStepInfo(Model, settings["Parameters"])


def get_strain(model_part):
    elem = model_part.Elements[1]
    process_info = model_part.ProcessInfo
    return elem.GetValuesOnIntegrationPoints(km.GREEN_LAGRANGE_STRAIN_TENSOR, process_info)[0]


def get_stress(model_part):
    elem = self.model_part.Elements[1]
    process_info = model_part.ProcessInfo
    return elem.GetValuesOnIntegrationPoints(km.CAUCHY_STRESS_TENSOR, process_info)[0]


def append_strain_stress_file(model_part, filename):
    with open(filename, 'a') as fo:
        fo.write("{} {}\n".format(get_strain(model_part)[0], get_stress(model_part)[0]))
        #fo.write("{} {}\n".format(get_strain(model_part), get_stress(model_part)))
    

def parameters_get_list_doubles(settings_list):
    olist = []
    for i in range(settings_list.size()):
        olist.append(settings_list[i].GetDouble())
    return olist


class WriteGeneralStepInfo(km.Process):

    def __init__(self, Model, settings):
        km.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.problem_name = "strain-stress"

    def ExecuteInitialize(self):
        filename = self.problem_name + ".out"
        if os.path.exists(filename):
            os.remove(filename)
        pass

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        append_strain_stress_file(self.model_part, self.problem_name + ".out")
        pass


    def ExecuteFinalize(self):
        pass


