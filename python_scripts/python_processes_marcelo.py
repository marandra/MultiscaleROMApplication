import KratosMultiphysics as km
import bisect
import os

def Factory(settings, Model):
    if(type(settings) != km.Parameters):
        raise Exception("expected input is Parameters object, encapsulating a json string")
    return ApplyCustomDisplacementProcess(Model, settings["Parameters"])


def get_strain(model_part):
    elem = model_part.Elements[1]
    process_info = model_part.ProcessInfo
    return elem.GetValuesOnIntegrationPoints(km.GREEN_LAGRANGE_STRAIN_TENSOR, process_info)[0]


def get_stress(model_part):
    print("ACA")
    elem = self.model_part.Elements[1]
    process_info = model_part.ProcessInfo
    return elem.GetValuesOnIntegrationPoints(km.CAUCHY_STRESS_TENSOR, process_info)[0]


def append_strain_stress_file(model_part, filename):
    with open(filename, 'a') as fo:
        #fo.write("{} {}\n".format(get_strain(model_part)[0], get_stress(model_part)[0]))
        fo.write("{} {}\n".format(get_strain(model_part), get_stress(model_part)))
    

def parameters_get_list_doubles(settings_list):
    olist = []
    for i in range(settings_list.size()):
        olist.append(settings_list[i].GetDouble())
    return olist


def get_multiplier(self):
    time = self.model_part.ProcessInfo[km.TIME]
    # interpolator is only valid within range defined in lookuptable
    # for times outside range, we force time to be the corresponding extreme.
    if(time <= self.lookuptable['time'][0]):
        time = self.lookuptable['time'][0]
    if(time >= self.lookuptable['time'][-1]):
        time = self.lookuptable['time'][-1]
    return self.time_interpolator[time]


class Interpolate(object):
  
    def __init__(self, x_list, y_list):
        if any([y - x <= 0 for x, y in zip(x_list, x_list[1:])]):
            raise ValueError("x_list must be in strictly ascending order!")
        x_list = self.x_list = list(map(float, x_list))
        y_list = self.y_list = list(map(float, y_list))
        intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
        self.slopes = [(y2 - y1)/(x2 - x1) for x1, x2, y1, y2 in intervals]

    def __getitem__(self, x):
        i = bisect.bisect(self.x_list, x) - 1
        if i >= len(self.slopes):
            i = len(self.slopes) - 1
        return self.y_list[i] + self.slopes[i] * (x - self.x_list[i])

        
class ApplyCustomDisplacementProcess(km.Process):

    def __init__(self, Model, settings):
        km.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        #self.variable = km.globals().get(settings["variable_name"].GetString)
        self.lookuptable = {'time': parameters_get_list_doubles(settings["lookuptable_time"]),
                            'mult': parameters_get_list_doubles(settings["lookuptable_mult"])}
        self.time_interpolator = Interpolate(self.lookuptable['time'],
                                             self.lookuptable['mult'])
        self.problem_name = "strain-stress"

    def ExecuteInitialize(self):
        filename = self.problem_name + ".out"
        if os.path.exists(filename):
            os.remove(filename)

        for node in self.model_part.Nodes:
             print(node)
             self.final_value = node.GetSolutionStepValue(km.DISPLACEMENT)

    def ExecuteInitializeSolutionStep(self):
        multiplier = get_multiplier(self)
        print("DEBUG Interpolate function. FACTOR {}".format(multiplier) )
        for node in self.model_part.Nodes:
             value = multiplier * self.final_value
             print(value)
             node.SetSolutionStepValue(km.DISPLACEMENT, 0, value)

    def ExecuteFinalizeSolutionStep(self):
        #append_strain_stress_file(self.model_part, self.problem_name + ".out")
        pass


    def ExecuteFinalize(self):
        print("DEBUG PROCESS ExecuteFinalize") 


