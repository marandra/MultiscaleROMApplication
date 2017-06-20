set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultiscaleROMApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/SolidMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## generate variables with the sources
set( KRATOS_MULTISCALE_ROM_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multiscale_rom_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_conditions_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_isotropic_damage_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/exponential_isotropic_damage_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_J2_plasticity_plane_strain_2D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_J2_plasticity_3D_law.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/linear_isotropic_damage_3D_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/rve_law.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/minimal_kinetic_2D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/minimal_kinetic_3D.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_str_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_std_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_bbar_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_hprom_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_hprom_J2_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_fbar_element.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/total_lagrangian_fbar_element.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/lagrange_multiplier_utility.h )

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosMultiscaleROMApplication SHARED ${KRATOS_MULTISCALE_ROM_APPLICATION_SOURCES})
target_link_libraries(KratosMultiscaleROMApplication KratosCore KratosSolidMechanicsApplication KratosStructuralMechanicsApplication)
set_target_properties(KratosMultiscaleROMApplication PROPERTIES PREFIX "")
set_target_properties(KratosMultiscaleROMApplication PROPERTIES COMPILE_DEFINITIONS "MULTISCALE_ROM_APPLICATION=EXPORT,API")
install(TARGETS KratosMultiscaleROMApplication DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

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
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultiscaleROMApplication.py" DESTINATION KratosMultiphysics )
