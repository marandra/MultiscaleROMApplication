set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultiscaleROMApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MultiscaleROMApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## generate variables with the sources
set( KRATOS_MULTISCALE_ROM_APPLICATION_SOURCES
    ## MAIN FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application_variables.cpp

    ## CUSTOM PYTHON
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multiscale_rom_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp

    ## UTILITIES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/lagrange_multiplier_utility.h

    ## CONSTITUTIVE LAW
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/rve_law.cpp

    ## CONDITIONS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/minimal_kinetic_2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/minimal_kinetic_3D.cpp

    ## ELEMENTS
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_custom_element.cpp
)

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMultiscaleROMApplication MODULE ${KRATOS_MULTISCALE_ROM_APPLICATION_SOURCES})

target_link_libraries(KratosMultiscaleROMApplication PRIVATE KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosMultiscaleROMApplication PROPERTIES COMPILE_DEFINITIONS "MULTISCALE_ROM_APPLICATION=EXPORT,API")
set_target_properties(KratosMultiscaleROMApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMultiscaleROMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMultiscaleROMApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN "*.json" PATTERN "*.mdpa" PATTERN "*.ref")
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN "*.csv")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

if(USE_COTIRE MATCHES ON)
    cotire(KratosMultiscaleROMApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosMultiscaleROMApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultiscaleROMApplication.py" DESTINATION KratosMultiphysics )
