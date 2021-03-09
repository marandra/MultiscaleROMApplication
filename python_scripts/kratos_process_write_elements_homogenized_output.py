import os
import KratosMultiphysics as km


def Factory(settings, Model):
    return WriteElementsHomogenizedOutput(settings["Parameters"], Model)


def homogenization_function(self):
    stress_ref = self.model_part.Elements[1].CalculateOnIntegrationPoints(
        km.CAUCHY_STRESS_VECTOR, self.model_part.ProcessInfo
        #km.STRESSES, self.model_part.ProcessInfo
    )
    nr_comp = len(stress_ref[0])
    stress_accum = [0.0] * nr_comp
    strain_accum = [0.0] * nr_comp
    volume = 0.0

    for e, elem in enumerate(self.model_part.Elements):
        stress = elem.CalculateOnIntegrationPoints(
            #km.STRESSES, self.model_part.ProcessInfo
            km.CAUCHY_STRESS_VECTOR, self.model_part.ProcessInfo
        )
        strain = elem.CalculateOnIntegrationPoints(
            km.STRAIN, self.model_part.ProcessInfo
        )
        weights = elem.CalculateOnIntegrationPoints(
            km.INTEGRATION_WEIGHT, self.model_part.ProcessInfo
        )
        for i, w in enumerate(weights):
            # used in HPROM case, to ignore GP
            if w == -1:
                continue
            for j in range(nr_comp):
                stress_accum[j] += stress[i][j] * w
                strain_accum[j] += strain[i][j] * w
            volume += w
    for i in range(nr_comp):
        stress_accum[i] /= volume
        strain_accum[i] /= volume
    return stress_accum, strain_accum


def write_strain_stress_header(filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    with open(filename, "w") as fo:
        fo.write(
            "#{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
            "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
                "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"  # strain
            )
        )  # stress
        fo.write(
            "#{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} "
            "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n".format(
                "strain XX",
                "YY",
                "ZZ",
                "XY",
                "YZ",
                "XZ",
                "stress XX",
                "YY",
                "ZZ",
                "XY",
                "YZ",
                "XZ",
            )
        )


def write_strain_stress(filename, strain, stress):
    line = (
        "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  "
        "{:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}  {:<+1.4e}\n".format(
            strain[0],
            strain[1],
            strain[2],
            strain[3],
            strain[4],
            strain[5],
            stress[0],
            stress[1],
            stress[2],
            stress[3],
            stress[4],
            stress[5],
        )
    )
    with open(filename, "a") as ofile:
        ofile.write(line)


class WriteElementsHomogenizedOutput(km.Process):
    def __init__(self, param, Model):
        km.Process.__init__(self)

        self.model_part = Model[param["model_part_name"].GetString()]
        self.filename = param["filename"].GetString()


    def ExecuteInitialize(self):
        write_strain_stress_header(self.filename)
        # WORKAROUND: added 0,0 row for consistency
        nr_comp = 6
        stress = [0.0] * nr_comp
        strain = [0.0] * nr_comp
        write_strain_stress(self.filename, strain, stress)

    def ExecuteFinalizeSolutionStep(self):
        stress, strain = homogenization_function(self)
        write_strain_stress(self.filename, strain, stress)
