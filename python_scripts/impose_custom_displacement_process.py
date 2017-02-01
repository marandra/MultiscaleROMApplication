import KratosMultiphysics as km
import bisect
import os


def Factory(settings, Model):
    if(type(settings) != km.Parameters):
        raise Exception("expected input is Parameters object, encapsulating a json string")
    return ApplyCustomDisplacementProcess(Model, settings["Parameters"])


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
        self.lookuptable = {'time': parameters_get_list_doubles(settings["lookuptable_time"]),
                            'mult': parameters_get_list_doubles(settings["lookuptable_mult"])}
        self.time_interpolator = Interpolate(self.lookuptable['time'],
                                             self.lookuptable['mult'])

    def ExecuteInitialize(self):
        node = self.model_part.Nodes[1]
        self.final_value = node.GetSolutionStepValue(km.DISPLACEMENT)

    def ExecuteInitializeSolutionStep(self):
        multiplier = get_multiplier(self)
        for node in self.model_part.Nodes:
             value = multiplier * self.final_value
             node.SetSolutionStepValue(km.DISPLACEMENT, 0, value)

    def ExecuteFinalizeSolutionStep(self):
        pass


    def ExecuteFinalize(self):
        pass

