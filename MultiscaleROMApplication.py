import KratosMultiphysics as KM
from KratosMultiscaleROMApplication import *
application = KratosMultiscaleROMApplication()
application_name = "KratosMultiscaleROMApplication"
application_folder = "MultiscaleROMApplication"

# The following lines are common for all applications
KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
