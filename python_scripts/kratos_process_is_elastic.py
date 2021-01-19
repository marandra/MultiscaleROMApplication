from pathlib import Path
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication


def Factory(settings, model):
    return IsElastic(settings["Parameters"], model)


class IsElastic(km.Process):
    def __init__(self, settings, model):
        km.Process.__init__(self)

        default_settings = km.Parameters(
            """{"model_part_name": "unset_model_part_name",
                "filename": "unset_filename"}"""
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.inelastic_flag = False
        self.model_part = model[settings["model_part_name"].GetString()]
        self.pf = Path(settings["filename"].GetString())

    def has_damaged_elements(self):
        for elem in self.model_part.Elements:
            flag = elem.CalculateOnIntegrationPoints(
                km.DAMAGE_VARIABLE, self.model_part.ProcessInfo
            )
            if True in [x > 0.0 for x in flag]:
                return True

            flag = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.StructuralMechanicsApplication.ACCUMULATED_PLASTIC_STRAIN,
                self.model_part.ProcessInfo,
            )
            if True in [x > 0.0 for x in flag]:
                return True

###########################################################
###########################################################

    def ExecuteFinalizeSolutionStep(self):
        self.time = self.model_part.ProcessInfo[km.TIME]
        if not self.inelastic_flag:
            self.inelastic_flag = self.has_damaged_elements()

    def ExecuteFinalize(self):
        if not self.inelastic_flag:
            # ELASTIC
            line = "1"
        else:
            # INELASTIC
            line = "0"
        self.pf.write_text(line)
