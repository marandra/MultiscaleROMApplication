import KratosMultiphysics as km
import os
import operator
import math

def Factory(settings, Model):
    return WriteElementsHomogenizedOutput(settings["Parameters"], Model)


def homogenization_function(self):
    stress_ref = self.model_part.Elements[1].GetValuesOnIntegrationPoints(self.stress,self.model_part.ProcessInfo)
    nr_comp = len(stress_ref[0])
    stress_accum = [0.0] * nr_comp
    strain_accum = [0.0] * nr_comp
    tensor_accum = [0.0] * nr_comp * nr_comp
    volume = 0.0

    for e, elem in enumerate(self.model_part.Elements):
        stress = elem.GetValuesOnIntegrationPoints(self.stress, self.model_part.ProcessInfo)
        strain = elem.GetValuesOnIntegrationPoints(self.strain, self.model_part.ProcessInfo)
        tensor = elem.GetValuesOnIntegrationPoints(self.tensor, self.model_part.ProcessInfo)
        weights = elem.GetValuesOnIntegrationPoints(km.INTEGRATION_WEIGHT, self.model_part.ProcessInfo)
        weights = [x[0] for x in weights] # to unpack received list-inside-list
        for i, w in enumerate(weights):
            # used in HPROM case, to ignore GP
            if w == -1:
                continue
            for j in range(nr_comp):
                stress_accum[j] += stress[i][j] * w
                strain_accum[j] += strain[i][j] * w
            for j in range(nr_comp * nr_comp):
                tensor_accum[j] += tensor[i][j] * w
            volume += w
    for i in range(nr_comp):
        stress_accum[i] /= volume
        strain_accum[i] /= volume
        tensor_accum[i] /= volume
    return stress_accum, strain_accum, tensor_accum

def compute_vonmisses_stress(hs):
    s = (hs[0] + hs[1] + hs[2]) / 3
    d = [hs[0] - s, hs[1] - s, hs[2] - s, hs[3], hs[4], hs[5]]
    dd = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] + \
         d[3] * d[3] + d[4] * d[4] + d[5] * d[5]
    vm = math.sqrt(dd)
    return vm

class WriteElementsHomogenizedOutput(km.Process):
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
        self.filename = param['filename'].GetString()
        #self.vname = param['variable_name'].GetString()
        #f = operator.attrgetter(self.vname)
        #self.Var = f(km)
        self.stress = km.CAUCHY_STRESS_VECTOR
        self.strain = km.GREEN_LAGRANGE_STRAIN_VECTOR
        self.tensor = km.CONSTITUTIVE_MATRIX

    def write_results(self, filename):
        homog_stress, homog_strain, const_tensor = homogenization_function(self)
        line = "{:<5} "\
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
               "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  " \
               "{:<+1.4e}  {:<+1.4e}  " \
               "{:<+1.4e}\n".format(
            0,
            homog_strain[ 0], homog_strain[ 1], homog_strain[ 2], homog_strain[ 3], homog_strain[ 4], homog_strain[ 5],
            homog_stress[ 0], homog_stress[ 1], homog_stress[ 2], homog_stress[ 3], homog_stress[ 4], homog_stress[ 5],
            const_tensor[ 0], const_tensor[ 1], const_tensor[ 2], const_tensor[ 3], const_tensor[ 4], const_tensor[ 5],
                              const_tensor[ 7], const_tensor[ 8], const_tensor[ 9], const_tensor[10], const_tensor[11],
                                                const_tensor[14], const_tensor[15], const_tensor[16], const_tensor[17],
                                                                  const_tensor[21], const_tensor[22], const_tensor[23],
                                                                                    const_tensor[28], const_tensor[29],
                                                                                                      const_tensor[35]
        )
        with open(filename, 'a') as ofile:
            ofile.write(line)

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
        with open(self.filename, "w") as fo:
            fo.write("#col: "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12}\n".format(
                     "2", "3", "4", "5" , "6", "7",  # strain
                     "8", "9", "10", "11" , "12", "13", # stress
                     "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24",
                     "25", "26", "27", "28", "29", "30", "31", "32", "33", "34"))
            fo.write("#     "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
                     "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
                      "strain XX", "YY", "ZZ", "XY" , "YZ", "XZ",
                      "stress XX", "YY", "ZZ", "XY" , "YZ", "XZ",
                      "CT 11", "12", "13", "14" , "15", "16", "22", "23", "24", "25" , "26",
                      "33", "34", "35", "36" , "44", "45", "46", "55", "56", "66"))

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        self.write_results(self.filename)

    def ExecuteFinalize(self):
        pass
