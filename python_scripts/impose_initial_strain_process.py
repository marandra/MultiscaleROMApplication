import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as msr
import bisect
import os


def Factory(settings, Model):
    if(type(settings) != km.Parameters):
        raise Exception("expected input is Parameters object, encapsulating a json string")
    return ImposeInitialStrainProcess(Model, settings["Parameters"])


def parameters_get_list_doubles(settings_list):
    olist = []
    for i in range(settings_list.size()):
        olist.append(settings_list[i].GetDouble())
    return olist


def get_scaling_factor(self):
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

        
class ImposeInitialStrainProcess(km.Process):

    def __init__(self, Model, settings):
        km.Process.__init__(self)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.lookuptable = {
            'time': parameters_get_list_doubles(settings["lookuptable_time"]),
            'mult': parameters_get_list_doubles(settings["lookuptable_mult"])}
        self.time_interpolator = Interpolate(
            self.lookuptable['time'], self.lookuptable['mult'])

        initial_strain_list = parameters_get_list_doubles(settings["initial_strain"])
        self.initial_strain = km.Vector(len(initial_strain_list))
        for i, s in enumerate(initial_strain_list):
            self.initial_strain[i] = s
        # TODO set it early so CL can check correct size. Not actually working.
        #self.model_part.ProcessInfo[msr.INITIAL_STRAIN_VECTOR] = self.initial_strain

    def ExecuteInitialize(self):
        self.model_part.ProcessInfo[km.INITIAL_STRAIN] = self.initial_strain

    def ExecuteInitializeSolutionStep(self):
        strain = get_scaling_factor(self) * self.initial_strain
        self.model_part.ProcessInfo[km.INITIAL_STRAIN] = strain

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteFinalize(self):
        pass

