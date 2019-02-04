import KratosMultiphysics as Kratos
import KratosMultiphysics.MultiscaleROMApplication as MSA
import os
import struct
import numpy


def Factory(settings, model):
    return LoadModesToProperties(settings["Parameters"], model)


def numpy_to_kratos(np_matrix):
    nr_rows = numpy.shape(np_matrix)[0]
    nr_cols = numpy.shape(np_matrix)[1]
    k_matrix = Kratos.Matrix(nr_rows, nr_cols)
    for r in range(nr_rows):
        for c in range(nr_cols):
            k_matrix[r, c] = np_matrix[r, c]
    return k_matrix


def kratos_to_numpy(k_matrix):
    nr_rows = k_matrix.Size1()
    nr_cols = k_matrix.Size2()
    np_matrix = numpy.empty((nr_rows, nr_cols))
    for r in range(nr_rows):
        for c in range(nr_cols):
           np_matrix[r, c] = k_matrix[r, c]
    return np_matrix


class LoadModesToProperties(Kratos.Process):
    def __init__(self, settings, model):
        Kratos.Process.__init__(self)
        default_settings = Kratos.Parameters("""
        {
            "mesh_id": 0,
            "model_part_name": "unset_model_part_name",
            "modes_filename": "unset_filename",
            "modes_file_format": "binary",
            "number_modes_to_load": 0
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[settings['model_part_name'].GetString()]
        self.modes_filename = settings['modes_filename'].GetString()
        self.modes_file_format = settings['modes_file_format'].GetString()
        self.nr_modes = settings['number_modes_to_load'].GetInt()

    def ExecuteInitialize(self):
        def read_modes(filename, file_format, nr_modes):
            if file_format == 'binary':
                modes = numpy.load(filename)[:, :nr_modes]
            else:
                modes = numpy.loadtxt(filename)[:, :nr_modes]
            return modes
        # Create global modes matrix
        modes_numpy = read_modes(self.modes_filename, self.modes_file_format, self.nr_modes)
        modes_matrix = numpy_to_kratos(modes_numpy)
        self.model_part.ProcessInfo[MSA.GLOBAL_MODES_MATRIX] = modes_matrix

        # Create LHS and RHS matrices
        print("WARNING: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx")
        print("WARNING: harcoded nr of DOFs in process")
        print("WARNING: XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxx")
        nr_nodes = self.model_part.NumberOfNodes(0)
        nr_dofs = 3 * nr_nodes
        self.model_part.ProcessInfo[MSA.RHS_MATRIX] = Kratos.Matrix(nr_dofs, self.nr_modes)
        self.model_part.ProcessInfo[MSA.LHS_MATRIX] = Kratos.Matrix(nr_dofs, nr_dofs)

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        pass

    def ExecuteFinalize(self):
        def save_matrix(filename, file_format, matrix):
            if file_format == 'binary':
                modes = numpy.save(filename, matrix)
            else:
                modes = numpy.savetxt(filename, matrix)
            return modes

        rhs_matrix = self.model_part.ProcessInfo[MSA.RHS_MATRIX]
        filename = self.modes_filename.rsplit(".", 1)[0]+ "_rhs"
        save_matrix(filename, self.modes_file_format, kratos_to_numpy(rhs_matrix))

        lhs_matrix = self.model_part.ProcessInfo[MSA.LHS_MATRIX]
        filename = self.modes_filename.rsplit(".", 1)[0]+ "_lhs"
        save_matrix(filename, self.modes_file_format, kratos_to_numpy(lhs_matrix))
