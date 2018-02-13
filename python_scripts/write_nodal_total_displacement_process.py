import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication as msr
import os
import struct

def Factory(settings, Model):
    if(type(settings) != km.Parameters):
        raise Exception("expected input is Parameters object, encapsulating a json string")
    return ComputeTotalDisplacementProcess(Model, settings["Parameters"])

def TotalDisplacement(node, initial_strain):
    # compute regular displacements
    comp_x =  initial_strain[0] * node.X0 + 0.5 * initial_strain[3] * node.Y0 + 0.5 * initial_strain[5] * node.Z0
    comp_y = 0.5 * initial_strain[3] * node.X0 + initial_strain[1] * node.Y0 + initial_strain[4] * node.Z0
    comp_z = 0.5 * initial_strain[5] * node.X0 + 0.5 * initial_strain[4] * node.Y0 + initial_strain[3] * node.Z0

    displ = node.GetSolutionStepValue(km.DISPLACEMENT)
    total_dis_X = comp_x + displ[0]
    total_dis_Y = comp_y + displ[1]
    total_dis_Z = comp_z + displ[2]

    # Total displacement
    total_displ = [node.Id, total_dis_X, total_dis_Y, total_dis_Z]

    return total_displ


class ComputeTotalDisplacementProcess(km.Process):

    def __init__(self, Model, settings):
        km.Process.__init__(self)

        default_settings = km.Parameters("""
        {
            "mesh_id": 0,
            "model_part_name": "unset_model_part_name",
            "filename": "unset_filename",
            "write_mode": "ascii"
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.filename = settings['filename'].GetString()
        self.write_mode = settings['write_mode'].GetString()

    def write_results(self, filename):
        def write_results_binary():
            with open(filename, 'wb') as ofile:
                initial_strain = self.model_part.ProcessInfo[km.INITIAL_STRAIN]
                for node in self.model_part.Nodes:
                    nodal_displ = TotalDisplacement(node, initial_strain)
                    #for v in total_displ:
                    for v in nodal_displ:
                        ofile.write(struct.pack('f', v)) # 'f'=float32
                    ofile.write(b'\n')

        def write_results_ascii():
            with open(filename, 'w') as ofile:
                initial_strain = self.model_part.ProcessInfo[km.INITIAL_STRAIN]
                for node in self.model_part.Nodes:
                    nodal_displ = TotalDisplacement(node, initial_strain)
                    for v in nodal_displ:
                        ofile.write("{:18.16f} ".format(v))
                    ofile.write("\n")

        if self.write_mode == "binary":
            write_results_binary()
        else:
            write_results_ascii()

    def ExecuteInitialize(self):
        try:
            os.remove(self.filename)
        except OSError:
            pass

    def ExecuteInitializeSolutionStep(self):
        self.timestep = "-{:.3f}".format(self.model_part.ProcessInfo[km.TIME])
        try:
            os.remove(self.filename + self.timestep)
        except OSError:
            pass
        #pass

    def ExecuteFinalizeSolutionStep(self):
        self.write_results(self.filename + self.timestep)

    def ExecuteFinalize(self):
        pass