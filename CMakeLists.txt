set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultiscaleROMApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/MultiscaleROMApplication )

## generate variables with the sources
set( KRATOS_MULTISCALE_ROM_APPLICATION_CORE
  ## MAIN FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application_variables.cpp

  ## UTILITIES
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/lagrange_multiplier_utility.h

  ## CONSTITUTIVE LAW
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/rve_law.cpp

  ## CONDITIONS
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/minimal_kinetic_vol_3D.cpp

  ## ELEMENTS
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/small_displacement_custom.cpp
)

## generate variables with the testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_MULTISCALE_ROM_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## generate variables with the sources
set( KRATOS_MULTISCALE_ROM_APPLICATION_PYTHON_INTERFACE
  ## CUSTOM PYTHON
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/multiscale_rom_python_application.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_constitutive_laws_to_python.cpp
)

add_library(KratosMultiscaleROMCore SHARED ${KRATOS_MULTISCALE_ROM_APPLICATION_CORE} ${KRATOS_MULTISCALE_ROM_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMultiscaleROMCore PUBLIC KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosMultiscaleROMCore PROPERTIES COMPILE_DEFINITIONS "MULTISCALE_ROM_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMultiscaleROMApplication MODULE THIN_LTO ${KRATOS_MULTISCALE_ROM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMultiscaleROMApplication PUBLIC KratosMultiscaleROMCore)
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
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN "*.csv")
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

if(USE_COTIRE MATCHES ON)
    cotire(KratosMultiscaleROMCore)
    cotire(KratosMultiscaleROMApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosMultiscaleROMCore DESTINATION libs )
install(TARGETS KratosMultiscaleROMApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/MultiscaleROMApplication.py" DESTINATION "KratosMultiphysics/MultiscaleROMApplication" RENAME "__init__.py")
