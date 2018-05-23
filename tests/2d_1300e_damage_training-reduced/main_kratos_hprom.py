from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as kmsr
import read_materials_process


class Kratos:
    def __init__(self, filename):
        self.model_part_rve = km.ModelPart("RVE")
        node = self.model_part_rve.CreateNewNode(1, 0.0, 0.0, 0.0)
        geom = km.Triangle2D3(node, node, node)
        Model = {"RVE" : self.model_part_rve}
        materials_rve  = km.Parameters("""
                   {
                       "Parameters": {
                               "materials_filename": "materials_rve.json"
                       }
               }
               """)
        read_materials_process.Factory(materials_rve, Model)
        # import rve_data json string
        with open (filename, "r") as myfile:
            rve_data = km.Parameters(myfile.read())
        print(self.model_part_rve)
        self.cl = kmsr.RVELaw(self.model_part_rve, rve_data)
        #cl_clone = cl.Clone()
        self.cl.Check(km.ModelPart("dummy").Properties[1], geom, self.model_part_rve.ProcessInfo)
        self.cl.InitializeMaterial(km.ModelPart("dummy").Properties[1], geom, km.Vector(3))

    def init_case(self, strain, path_out="."):
        self.path_out = path_out
        #nr_comp = self.cl.GetStrainSize()
        nr_comp = 3
        self.init_strain_macro = km.Vector(nr_comp)
        for i, e in enumerate(strain):
            self.init_strain_macro[i] = e
        print(self.init_strain_macro)
    
    def run(self):
        node1 = self.model_part_rve.CreateNewNode(1,0.0,0.0,0.0)
        geom = km.Triangle2D3(node1, node1, node1) # create point geom
        #nr_comp = self.cl.GetStrainSize()
        nr_comp = 3
        homog_stress = km.Vector(nr_comp)
        homog_constit = km.Matrix(nr_comp, nr_comp)
        cl_params = km.ConstitutiveLawParameters()
        cl_options = km.Flags()
        cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
        cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
        cl_params.SetOptions(cl_options)
        cl_params.SetStressVector(homog_stress)
        cl_params.SetConstitutiveMatrix(homog_constit)
        cl_params.SetMaterialProperties(self.model_part_rve.Properties[1])
    
        nr_timesteps = 250
        t = dt = 1. / nr_timesteps
        fo = open("{}/hs_hprom.dat".format(self.path_out) ,'w')
        while (t <= 1. + dt / 10.):
            print("time {:.3f}".format(t))
            self.model_part_rve.CloneTimeStep(t)
            strain_macro = t * self.init_strain_macro
            cl_params.SetStrainVector(strain_macro)
            self.cl.CalculateMaterialResponseCauchy(cl_params)
            self.cl.FinalizeSolutionStep(km.ModelPart("dummy").Properties[1], geom,
                                    km.Vector(3), self.model_part_rve.ProcessInfo)
            cl_params.GetStressVector(homog_stress)
            print("{}: {}".format(t, homog_stress))
            cl_params.GetConstitutiveMatrix(homog_constit)
            print("{}: {}".format(t, homog_constit))
            for ih in homog_stress:
                 fo.write("{:.16e} ".format(ih))
            fo.write("\n")
            fo.flush()
            t += dt
            #print("BREAK"); exit()
        fo.close()
