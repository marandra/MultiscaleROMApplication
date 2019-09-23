import KratosMultiphysics as km
import KratosMultiphysics.MultiscaleROMApplication.io_utilities as io_utilities


def Factory(settings, Model):
    return WriteElementsHomogenizedOutput(settings["Parameters"], Model)


def homogenization_function(self):
    stress_ref = self.model_part.Elements[1].GetValuesOnIntegrationPoints(self.stress, self.model_part.ProcessInfo)
    nr_comp = len(stress_ref[0])
    stress_accum = [0.0] * nr_comp
    strain_accum = [0.0] * nr_comp
    tensor_accum = [0.0] * nr_comp * nr_comp
    volume = 0.0

    for e, elem in enumerate(self.model_part.Elements):
        stress = elem.GetValuesOnIntegrationPoints(self.stress, self.model_part.ProcessInfo)
        strain = elem.GetValuesOnIntegrationPoints(self.strain, self.model_part.ProcessInfo)
        tensor = elem.GetValuesOnIntegrationPoints(self.tensor, self.model_part.ProcessInfo)
        weights = elem.GetValuesOnIntegrationPoints(km.INTEGRATION_WEIGHT, self.model_part.ProcessInfo)
        weights = [x[0] for x in weights]  # to unpack received list-inside-list
        for i, w in enumerate(weights):
            # used in HPROM case, to ignore GP
            if w == -1:
                continue
            for j in range(nr_comp):
                stress_accum[j] += stress[i][j] * w
                strain_accum[j] += strain[i][j] * w
            for j in range(nr_comp * nr_comp):
                tensor_accum[j] += tensor[i][j] * w
            volume += w
    for i in range(nr_comp):
        stress_accum[i] /= volume
        strain_accum[i] /= volume
        tensor_accum[i] /= volume
    return stress_accum, strain_accum, tensor_accum


class WriteElementsHomogenizedOutput(km.Process):
    def __init__(self, param, Model):
        km.Process.__init__(self)

        self.model_part = Model[param['model_part_name'].GetString()]
        self.filename = param['filename'].GetString()
        self.stress = km.CAUCHY_STRESS_VECTOR
        self.strain = km.GREEN_LAGRANGE_STRAIN_VECTOR
        self.tensor = km.CONSTITUTIVE_MATRIX

    def ExecuteInitialize(self):
        io_utilities.write_strain_stress_header(self.filename)
        # WORKAROUND: added 0,0 row for consistency
        nr_comp = 6
        stress = [0.0] * nr_comp
        strain = [0.0] * nr_comp
        io_utilities.write_strain_stress(self.filename, strain, stress)

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteAfterOutputStep(self):
        pass

    def ExecuteBeforeOutputStep(self):
        pass

    def ExecuteBeforeSolutionLoop(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        stress, strain, const_tensor = homogenization_function(self)
        io_utilities.write_strain_stress(self.filename, strain, stress)

    def ExecuteFinalize(self):
        pass
