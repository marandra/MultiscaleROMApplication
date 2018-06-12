from __future__ import print_function, absolute_import, division
import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as kmsr
import read_materials_process

if __name__ == "__main__":
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
        mystr=myfile.read()
        rve_data = km.Parameters(mystr)
    #from load_rve_data import rve_data

    cl = kmsr.RVELaw(model_part_rve, rve_data)
    cl_clone = cl.Clone()
    print(cl)
    print(cl_clone)

    cl.Check(km.ModelPart("dummy").Properties[1], geom, model_part_rve.ProcessInfo)
    cl.InitializeMaterial(km.ModelPart("dummy").Properties[1], geom, km.Vector(3))

    nr_comp = cl.GetStrainSize()
    # creation and init
    init_strain_macro = km.Vector(nr_comp)
    homog_stress = km.Vector(nr_comp)
    homog_constit = km.Matrix(nr_comp, nr_comp)
    for i in range(nr_comp):
        init_strain_macro[i] = 0.
        # homog_stress[i] = 0.
        # for j in range(nr_comp):
        #     homog_constit[i, j] = 0.
    init_strain_macro[0] = 0.001

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
    fo=open("homog_stress.dat",'w')
    while (t <= 1. + dt / 10.):
        model_part_rve.CloneTimeStep(t)
        strain_macro = t * init_strain_macro
        cl_params.SetStrainVector(strain_macro)
        cl.CalculateMaterialResponseCauchy(cl_params)
        cl.FinalizeSolutionStep(km.ModelPart("dummy").Properties[1], geom,
                                km.Vector(3), model_part_rve.ProcessInfo)
        # Print output
        modes_weights = km.Vector(10)
        print(cl.Has(kmsr.REDUCED_MODES_WEIGHTS))
        print(cl.GetValue(kmsr.REDUCED_MODES_WEIGHTS, modes_weights))
        cl_params.GetStressVector(homog_stress)
        print("{}: {}".format(t, homog_stress))
        # This does not work:
        #homog_stress = cl_params.GetStressVector()
        #print("{}: {}".format(t, homog_stress))
        cl_params.GetConstitutiveMatrix(homog_constit)
        print("{}: {}".format(t, homog_constit))
        t += dt
        fo.write("{}\n".format(homog_stress[0]))
    fo.close()

