import KratosMultiphysics as km


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
            "model_part_name": "unset_model_part_name"
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]

    def ExecuteFinalizeSolutionStep(self):
        initial_strain = self.model_part.ProcessInfo[km.INITIAL_STRAIN]
        for node in self.model_part.Nodes:
            total_displ = TotalDisplacement(node, initial_strain)
            node.SetSolutionStepValue(km.LAGRANGE_DISPLACEMENT, 0, total_displ)
