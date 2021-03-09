"""
Kratos process to write xdmf time series of model.
IP values are averaged in cell (one value per cell)
"""
import numpy
import KratosMultiphysics as km
import meshio


def kratos_vector_to_numpy_array(_vector):
    """Convert Kratos Vector to numpy array.

    Arguments:
        _vector {KratosMultiphysics.Vector} -- Vector to be converted

    Returns:
        numpy.array -- the array with the values of _vector
    """
    _list = []
    for value in _vector:
        _list.append(value)
    return numpy.array(_list)


def Factory(settings, Model):
    return WriteXdmf(settings["Parameters"], Model)


class WriteXdmf(km.Process):
    """Write timeseries xdmf

    Arguments:
        km {Process} -- Kratos process
    """

    def __init__(self, settings, Model):
        km.Process.__init__(self)

        default_settings = km.Parameters(
            """
        {
            "model_part_name": "unset_model_part_name",
            "xdmf_filename": "model_field_output.xdmf",
            "mdpa_filename": "model.mdpa"
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.xdmf_fname = settings["xdmf_filename"].GetString()
        self.mdpa_fname = settings["mdpa_filename"].GetString()
        self.cell_data_list = []
        self.point_data_list = []

    def get_averaged_cell_data_scalar(self, kratos_label):
        """Average cell values and return an array.

        Returns:
            numpy.array -- One value per cell. Format compatible with meshio
        """
        # damage:  KratosMultiphysics.DAMAGE_VARIABLE
        # energy:  KratosMultiphysics.STRAIN_ENERGY
        _list = []
        for elem in self.model_part.Elements:
            values = elem.CalculateOnIntegrationPoints(
                kratos_label, self.model_part.ProcessInfo
            )
            _list.append(numpy.mean(values))
            arr = numpy.array(_list).reshape((-1, 1))
        return arr

    def get_averaged_cell_data_tensor(self, kratos_label):
        """Average cell values and return an array.

        Returns:
            numpy.array -- One value per cell. Format compatible with meshio
        """
        # KratosMultiphysics.CAUCHY_STRESS_VECTOR
        # KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR
        _list = []
        for elem in self.model_part.Elements:
            values = elem.CalculateOnIntegrationPoints(
                kratos_label, self.model_part.ProcessInfo
            )
            values = [kratos_vector_to_numpy_array(v) for v in values]
            _list.append(numpy.mean(values, axis=0))
        return numpy.array(_list)

    def get_point_data_reaction(self):
        """Create array of nodal fields.
        Arguments:
            field {str} -- "FLUCTUANT_DISPLACEMENT" or "TOTAL_DISPLACEMENT"
        Returns:
            array -- nodal values,  formatted for meshio [[x0, y0, z0], [x1, y1, z1], [x2. y2, z2]]
        """
        _list = []
        for node in self.model_part.Nodes:
            displ = node.GetSolutionStepValue(km.REACTION)
            _list.append([displ[0], displ[1], displ[2]])
        return numpy.array(_list)

    def get_point_data_displacement_fluct(self):
        """Create array of nodal fields.
        Arguments:
            field {str} -- "FLUCTUANT_DISPLACEMENT" or "TOTAL_DISPLACEMENT"
        Returns:
            array -- nodal values,  formatted for meshio [[x0, y0, z0], [x1, y1, z1], [x2. y2, z2]]
        """
        _list = []
        for node in self.model_part.Nodes:
            displ = node.GetSolutionStepValue(km.DISPLACEMENT)
            _list.append([displ[0], displ[1], displ[2]])
        return numpy.array(_list)

    def get_point_data_displacement_total(self):
        """Create array of nodal fields.
        Arguments:
            field {str} -- "FLUCTUANT_DISPLACEMENT" or "TOTAL_DISPLACEMENT"
        Returns:
            array -- nodal values,  formatted for meshio [[x0, y0, z0], [x1, y1, z1], [x2. y2, z2]]
        """
        _list = []

        for elem in self.model_part.Elements:
            values = elem.CalculateOnIntegrationPoints(
                km.INITIAL_STRAIN_VECTOR, self.model_part.ProcessInfo
            )
            break
        strain = values[0]
        for node in self.model_part.Nodes:
            displ = node.GetSolutionStepValue(km.DISPLACEMENT)
            s_xx = strain[0]
            s_yy = strain[1]
            s_zz = strain[2]
            s_xy = 0.5 * strain[3]
            s_yz = 0.5 * strain[4]
            s_xz = 0.5 * strain[5]
            comp_x = s_xx * node.X0 + s_xy * node.Y0 + s_xz * node.Z0
            comp_y = s_xy * node.X0 + s_yy * node.Y0 + s_yz * node.Z0
            comp_z = s_xz * node.X0 + s_yz * node.Y0 + s_zz * node.Z0
            total_displ_x = comp_x + displ[0]
            total_displ_y = comp_y + displ[1]
            total_displ_z = comp_z + displ[2]
            _list.append([total_displ_x, total_displ_y, total_displ_z])
        return numpy.array(_list)

    def ExecuteInitialize(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteFinalize(self):
        """
        Write point and cell data (specific for meshio)
        """
        mesh = meshio.read(self.mdpa_fname)
        cells = []
        for cell_block in mesh.cells:
            element_type = cell_block[0]
            # if "hexa" in element_type or "wedge" in element_type:
            if "line8" in element_type:
                cells.append(meshio.CellBlock("hexahedron", cell_block[1]))
            if "line6" in element_type:
                cells.append(meshio.CellBlock("wedge", cell_block[1]))
        meshio.write_points_cells(self.xdmf_fname, mesh.points, cells)
        with meshio.xdmf.TimeSeriesWriter(self.xdmf_fname) as writer:
            writer.write_points_cells(mesh.points, cells)
            displ_fluct = self.get_point_data_displacement_fluct()
            displ_total = self.get_point_data_displacement_total()
            reaction = self.get_point_data_reaction()
            damage = self.get_averaged_cell_data_scalar(km.DAMAGE_VARIABLE)
            energy = self.get_averaged_cell_data_scalar(km.STRAIN_ENERGY)
            stress = self.get_averaged_cell_data_tensor(km.CAUCHY_STRESS_VECTOR)
            #strain = self.get_averaged_cell_data_tensor(km.GREEN_LAGRANGE_STRAIN_VECTOR)
            strain = self.get_averaged_cell_data_tensor(km.STRAIN)
            for elem in self.model_part.Elements:
                values = elem.CalculateOnIntegrationPoints(
                    km.INITIAL_STRAIN_VECTOR, self.model_part.ProcessInfo
                )
                break
            initial_strain = values[0]
            strain_fluctuant = strain - initial_strain

            point_data = {
                "DISPLACEMENT_FLUCT": displ_fluct,
                "DISPLACEMENT": displ_total,
                "REACTION": reaction,
            }
            cell_data = {
                "DAMAGE": damage,
                "ENERGY": energy,
                "STRAIN_FLUCT": strain_fluctuant,
                "STRESS": stress,
            }
            writer.write_data(
                self.model_part.ProcessInfo[km.STEP],
                point_data=point_data,
                cell_data=cell_data,
            )
