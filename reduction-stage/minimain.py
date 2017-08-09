from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as kmsr
import read_materials_process
import configparser

if __name__ == "__main__":

    # parse configuration file
    conf = configparser.ConfigParser() 
    conf.read("reduced_bases.cfg")
    nr_modes = int(conf['Parameters']['nr_active_modes']) 


    model_part_rve = km.ModelPart("RVE")
    node1 = model_part_rve.CreateNewNode(1,0.0,0.0,0.0)
    geom = km.Triangle2D3(node1, node1, node1) # create point geom
    Model = {"RVE" : model_part_rve}
    materials_rve  = km.Parameters("""
               {
                   "Parameters": {
                           "materials_filename": "materials_rve.json"
                   }
           }
           """)
    read_materials_process.Factory(materials_rve, Model)
    #read_materials_process.ReadMaterialsProcess(Model, materials_rve)

    # import rve_data json string
    with open ("rve.json", "r") as myfile:
        rve_data = km.Parameters(myfile.read())

    cl = kmsr.RVELaw(model_part_rve, rve_data)
    cl_clone = cl.Clone()

    cl.Check(km.ModelPart("dummy").Properties[1], geom, model_part_rve.ProcessInfo)
    cl.InitializeMaterial(km.ModelPart("dummy").Properties[1], geom, km.Vector(3))

    nr_comp = cl.GetStrainSize()
    # creation and init
    init_strain_macro = km.Vector(nr_comp)
    homog_stress = km.Vector(nr_comp)
    homog_constit = km.Matrix(nr_comp, nr_comp)
    # trajectory 31:
    init_strain_macro[0] = 0.001
    init_strain_macro[1] = 0.0
    init_strain_macro[2] = 0.001
    init_strain_macro[3] = 0.0
    init_strain_macro[4] = 0.0
    init_strain_macro[5] = 0.001

    cl_params = km.ConstitutiveLawParameters()
    cl_options = km.Flags()
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
    cl_params.SetOptions(cl_options)
    cl_params.SetStressVector(homog_stress)
    cl_params.SetConstitutiveMatrix(homog_constit)
    cl_params.SetMaterialProperties(model_part_rve.Properties[1])

    nr_timesteps = 250
    t = dt = 1. / nr_timesteps
    fo = open("homog_stress.dat",'w')
    while (t <= 1. + dt / 10.):
        print("time {:.3f}".format(t))
        model_part_rve.CloneTimeStep(t)
        strain_macro = t * init_strain_macro
        cl_params.SetStrainVector(strain_macro)
        cl.CalculateMaterialResponseCauchy(cl_params)
        cl.FinalizeSolutionStep(km.ModelPart("dummy").Properties[1], geom,
                                km.Vector(3), model_part_rve.ProcessInfo)
        # Print output
        #modes_weights = km.Vector(nr_modes)
        #print(cl.GetValue(kmsr.REDUCED_MODES_WEIGHTS, modes_weights))
        cl_params.GetStressVector(homog_stress)
        #homog_stress = cl_params.GetStressVector()
        #print("{}: {}".format(t, homog_stress))
        #cl_params.GetConstitutiveMatrix(homog_constit)
        #print("{}: {}".format(t, homog_constit))
        t += dt
        fo.write("{}\n".format(homog_stress[0]))
    fo.close()
