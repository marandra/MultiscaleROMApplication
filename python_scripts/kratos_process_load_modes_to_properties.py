import numpy
from pathlib import Path
import KratosMultiphysics as Kratos
import KratosMultiphysics.MultiscaleROMApplication as MSA
from common import Common


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


def numpy_to_kratos_vector(np_vector):
    nr_rows = numpy.shape(np_vector)[0]
    k_vector = Kratos.Vector(nr_rows)
    for r in range(nr_rows):
        k_vector[r] = np_vector[r]
    return k_vector


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
        default_settings = Kratos.Parameters(
            """
        {
            "mesh_id": 0,
            "model_part_name": "unset_model_part_name",
            "global_index_filename": "unset_global_index_filename",
            "number_modes_to_load": 0,
            "root_path": "unset_path"
        }
        """
        )
        settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[settings["model_part_name"].GetString()]
        self.global_index_filename = settings["global_index_filename"].GetString()
        self.nr_modes = settings["number_modes_to_load"].GetInt()
        self.common = Common(Path(settings["root_path"].GetString()))

    def ExecuteInitialize(self):
        # Create global modes matrix
        modes_numpy = self.common.get_dataset("BASES", "STRAIN")[:, : self.nr_modes]
        modes_matrix = numpy_to_kratos(modes_numpy)
        self.model_part.ProcessInfo[MSA.GLOBAL_MODES_MATRIX] = modes_matrix

        # Create output matrix
        nr_dofs = 3 * self.model_part.NumberOfNodes(0)
        self.model_part.ProcessInfo[MSA.RHS_MATRIX] = Kratos.Matrix(
            nr_dofs, self.nr_modes
        )

        # Load global starting elements index vector
        global_index_numpy = numpy.loadtxt(self.global_index_filename)
        global_index_vector = numpy_to_kratos_vector(global_index_numpy)
        self.model_part.ProcessInfo[MSA.GLOBAL_INDEX_VECTOR] = global_index_vector

    def ExecuteFinalize(self):
        rhs_matrix = self.model_part.ProcessInfo[MSA.RHS_MATRIX]
        self.common.set_dataset(
            kratos_to_numpy(rhs_matrix), "CORRELATION", "STRAIN", self.nr_modes
        )
