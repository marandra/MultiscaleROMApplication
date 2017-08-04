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
    rve_data = km.Parameters(open("rve.json", 'r').read())
    cl = kmsr.RVELaw(model_part_rve, rve_data)
    cl_clone = cl.Clone()
    print(cl)
    print(cl_clone)

    cl.Check(km.ModelPart("dummy").Properties[1], geom, model_part_rve.ProcessInfo)
    cl.InitializeMaterial(km.ModelPart("dummy").Properties[1], geom, km.Vector(3))

    nr_comp = cl.GetStrainSize()
    init_strain_macro = km.Vector(6)
    init_strain_macro[0] = 0.001
    init_strain_macro[1] = 0.
    init_strain_macro[2] = 0.
    init_strain_macro[3] = 0.
    init_strain_macro[4] = 0.
    init_strain_macro[5] = 0.
    cl_params = km.ConstitutiveLawParameters()
    cl_options = km.Flags()
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_STRESS, True)
    cl_options.Set(km.ConstitutiveLaw.COMPUTE_CONSTITUTIVE_TENSOR, True)
    cl_params.SetOptions(cl_options)
    cl_params.SetStressVector(km.Vector(nr_comp))
    cl_params.SetConstitutiveMatrix(km.Matrix(nr_comp, nr_comp))
    cl_params.SetMaterialProperties(model_part_rve.Properties[1])

    nr_timesteps = 250
    t = dt = 1. / nr_timesteps
    while (t <= 1. + dt / 10.):
        model_part_rve.CloneTimeStep(t)
        strain = t * init_strain_macro
        cl_params.SetStrainVector(strain)
        cl.CalculateMaterialResponseCauchy(cl_params)
        CM = cl_params.GetConstitutiveMatrix()
        stress = cl_params.GetStressVector()
        cl.FinalizeSolutionStep(km.ModelPart("dummy").Properties[1], geom,
                                km.Vector(3), model_part_rve.ProcessInfo)
        print("{}: {}".format(t, stress))
        t += dt
