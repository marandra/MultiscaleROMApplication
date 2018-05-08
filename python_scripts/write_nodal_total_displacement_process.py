import KratosMultiphysics as km
#import KratosMultiphysics.MultiscaleROMApplication as msr
#import os
#import struct

def Factory(settings, Model):
    if(type(settings) != km.Parameters):
        raise Exception("expected input is Parameters object, encapsulating a json string")
    return ComputeTotalDisplacementProcess(Model, settings["Parameters"])

def TotalDisplacement(node, strain):
    s_xx = strain[0]
    s_yy = strain[1]
    s_zz = strain[2]
    s_xy = 0.5 * strain[3]
    s_yz = 0.5 * strain[4]
    s_xz = 0.5 * strain[5]
    comp_x = s_xx * node.X0 + s_xy * node.Y0 + s_xz * node.Z0
    comp_y = s_xy * node.X0 + s_yy * node.Y0 + s_yz * node.Z0
    comp_z = s_xz * node.X0 + s_yz * node.Y0 + s_zz * node.Z0
    displ = node.GetSolutionStepValue(km.DISPLACEMENT)
    total_displ_X = comp_x + displ[0]
    total_displ_Y = comp_y + displ[1]
    total_displ_Z = comp_z + displ[2]
    return [total_displ_X, total_displ_Y, total_displ_Z]


class ComputeTotalDisplacementProcess(km.Process):

    def __init__(self, Model, settings):
        km.Process.__init__(self)

        default_settings = km.Parameters("""
        {
            "mesh_id": 0,
            "model_part_name": "unset_model_part_name",
        }
        """)
        #"filename": "unset_filename",
        #"write_mode": "ascii"
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        #self.filename = settings['filename'].GetString()
        #self.write_mode = settings['write_mode'].GetString()

    #def write_results(self, filename):
    #    def write_results_binary():
    #        with open(filename, 'wb') as ofile:
    #            initial_strain = self.model_part.ProcessInfo[km.INITIAL_STRAIN]
    #            for node in self.model_part.Nodes:
    #                nodal_displ = TotalDisplacement(node, initial_strain)
    #                #for v in total_displ:
    #                for v in nodal_displ:
    #                    ofile.write(struct.pack('f', v)) # 'f'=float32
    #                ofile.write(b'\n')

    #    def write_results_ascii():
    #        with open(filename, 'w') as ofile:
    #            initial_strain = self.model_part.ProcessInfo[km.INITIAL_STRAIN]
    #            for node in self.model_part.Nodes:
    #                nodal_displ = TotalDisplacement(node, initial_strain)
    #                for v in nodal_displ:
    #                    ofile.write("{:18.16f} ".format(v))
    #                ofile.write("\n")

    #    if self.write_mode == "binary":
    #        write_results_binary()
    #    else:
    #        write_results_ascii()

    #def ExecuteInitialize(self):
    #    try:
    #        os.remove(self.filename)
    #    except OSError:
    #        pass

    #def ExecuteInitializeSolutionStep(self):
    #    #self.timestep = "-{:.3f}".format(self.model_part.ProcessInfo[km.TIME])
    #    #try:
    #    #    os.remove(self.filename + self.timestep)
    #    #except OSError:
    #    #    pass
    #    pass

    def ExecuteFinalizeSolutionStep(self):
        #self.write_results(self.filename + self.timestep)
        initial_strain = self.model_part.ProcessInfo[km.INITIAL_STRAIN]
        for node in self.model_part.Nodes:
            total_displ = TotalDisplacement(node, initial_strain)
            node.SetSolutionStepValue(km.LAGRANGE_DISPLACEMENT, 0, total_displ)

    #def ExecuteFinalize(self):
    #    pass