from pathlib import Path
import KratosMultiphysics as km
import KratosMultiphysics.StructuralMechanicsApplication


def Factory(settings, model):
    return IsElastic(settings["Parameters"], model)


class IsElastic(km.Process):
    def __init__(self, settings, model):
        km.Process.__init__(self)

        default_settings = km.Parameters(
            """{"model_part_name": "unset_model_part_name"}"""
        )
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.inelastic_flag = False

    def has_damaged_elements(self):
        for elem in self.model_part.Elements:
            flag = elem.CalculateOnIntegrationPoints(
                km.DAMAGE_VARIABLE, self.model_part.ProcessInfo
            )
            if True in [x > 0.0 for x in flag]:
                return True

            flag = elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.StructuralMechanicsApplication.ACCUMULATED_PLASTIC_STRAIN, self.model_part.ProcessInfo,
            )
            if True in [x > 0.0 for x in flag]:
                return True


    ###########################################################
    ###########################################################

    def ExecuteInitialize(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        self.time = self.model_part.ProcessInfo[km.TIME]
        if not self.inelastic_flag:
            self.inelastic_flag = self.has_damaged_elements()

    def ExecuteFinalize(self):
        #time = self.model_part.ProcessInfo[km.TIME]
        pf = Path("is_elastic.dat")
        #pf.touch(exist_ok=True)
        if not self.inelastic_flag:
            # ELASTIC
            line = "1"
        else:
            # INELASTIC
            line = "0"
        pf.write_text(line)
       
