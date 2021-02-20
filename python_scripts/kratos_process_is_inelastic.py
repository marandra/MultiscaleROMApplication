from pathlib import Path
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication as sm


def Factory(settings, model):
    return IsInelastic(settings["Parameters"], model)


class IsInelastic(km.Process):
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

    def has_damage(self):
        for elem in self.model_part.Elements:
            flag = elem.CalculateOnIntegrationPoints(
                km.DAMAGE_VARIABLE, self.model_part.ProcessInfo
            )
            if True in [x > 0.0 for x in flag]:
                return True
        return False

###########################################################
###########################################################

    def ExecuteFinalizeSolutionStep(self):
        self.time = self.model_part.ProcessInfo[km.TIME]
        self.inelastic_flag = self.has_damage()

    def ExecuteFinalize(self):
        if not self.inelastic_flag:
            # ELASTIC
            line = "1"
        else:
            # INELASTIC
            line = "0"
        self.pf.write_text(line)
