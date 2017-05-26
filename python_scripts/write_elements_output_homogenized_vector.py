import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as msr
#import bisect
import os
import operator

def Factory(settings, Model):
    return WriteElementsOutputHomogenizedVector(settings["Parameters"], Model)

def parameters_get_list_int(ilist):
    olist = []
    for i in range(ilist.size()):
        olist.append(ilist[i].GetInt())
    return olist

def homogenization_function(self):

    var_ref = self.model_part.Elements[1].GetValuesOnIntegrationPoints(self.Var,self.model_part.ProcessInfo)
    homog_comp = var_ref[0].__len__()

    var_acum =[0.0]*homog_comp
    volume = 0.0

    for elem in self.model_part.Elements:

        var_elem = elem.GetValuesOnIntegrationPoints(self.Var,self.model_part.ProcessInfo)
        weights= elem.GetValuesOnIntegrationPoints(msr.GAUSS_WEIGHTS,self.model_part.ProcessInfo)

        for iVar in range(weights.__len__()):
            for jVar in range(homog_comp):

                var_acum[jVar] = var_acum[jVar] + var_elem[iVar][jVar]*weights[iVar][0]

            volume += weights[iVar][0]

    for iComp in range(homog_comp):
        var_acum[iComp] /= volume

    #print(volume)
    #print(var_acum)
    return var_acum

class WriteElementsOutputHomogenizedVector(km.Process):
    def __init__(self, param, Model):
        self.model_part = Model[param['model_part_name'].GetString()]
        self.filename = param['filename'].GetString()
        self.vname = param['variable_name'].GetString()
        f = operator.attrgetter(self.vname)
        self.Var = f(km)

    def write_results(self, filename):
        with open(filename, 'w') as ofile:

            homog_value = homogenization_function(self)
            #with open(filename, 'wb') as ofile:

            #process_info = self.model_part.ProcessInfo
            #for elem in self.model_part.Elements:
            #    variables = elem.GetValuesOnIntegrationPoints(self.Var, process_info)
            for v in homog_value:
                #print(v)
                ofile.write("{:18.16f}\n".format(v))
                #ofile.write("{:18.16f}\n".format(v[4]))
                #ofile.write("{:18.16f}\n".format(v[8]))
                #ofile.write("{:18.16f}\n".format(v[1]))
                #ofile.write("{:18.16f}\n".format(v[5]))
                #ofile.write("{:18.16f}\n".format(v[2]))
                    #ofile.write(struct.pack('f', v[0])) #  'f'=float32
                    #ofile.write(b'\n')
                    #print(variable)
                    #ofile.write("{: .3e} {: .3e} {: .3e}  {: .3e} {: .3e} {: .3e}"  #"  {: .3e} {: .3e} {: .3e} {: .3e}"
                    #    .format(
                    #        CX[0], CX[3], CX[1],  # CX[0],
                    #        CY[0], CY[3], CY[1],  # CY[0],
                    #        ))
                    #ofile.write("\n")

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass
            #with open(self.filename, 'a') as ofile:
            #    ofile.write("#{:<32}  {:<32}\n".format(self.vname1, self.vname2))
            #    ofile.write("#{:<10} {:<10} {:<10}" #" {:<10} {:<10}"
            #                "  {:<10} {:<10} {:<10}" #" {:<10} {:<10}"
            #        .format(
            #        "Comp XX", "Comp YY", "Comp XY", #"Comp 4", "Comp 5",
            #        "Comp XX", "Comp YY", "Comp XY", #"Comp 4", "Comp 5"
            #        ))
            #    ofile.write("\n")
            #

            #self.write_results()
            #self.Tn = self.Model.ProcessInfo[km.TIME]

    def ExecuteInitializeSolutionStep(self):
        self.timestep = "-{:.3f}".format(self.model_part.ProcessInfo[km.TIME])
        try:
            os.remove(self.filename + self.timestep)
        except OSError:
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
        self.write_results(self.filename + self.timestep)


    def ExecuteFinalize(self):
        pass