set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMultiscaleROMApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${KRATOS_SOURCE_DIR}/applications/StructuralMechanicsApplication )

## MultiscaleROMApplications sources
file(GLOB_RECURSE KRATOS_MULTISCALE_ROM_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/multiscale_rom_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_constitutive/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/builders_and_solvers/*.h*
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.h*
)

## MultiscaleROMApplications Python interface sources
file(GLOB_RECURSE KRATOS_MULTISCALE_ROM_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

## MultiscaleROMApplication testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE KRATOS_MULTISCALE_ROM_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosMultiscaleROMCore SHARED ${KRATOS_MULTISCALE_ROM_APPLICATION_CORE} ${KRATOS_MULTISCALE_ROM_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosMultiscaleROMCore PUBLIC KratosCore KratosStructuralMechanicsCore)
set_target_properties(KratosMultiscaleROMCore PROPERTIES COMPILE_DEFINITIONS "MULTISCALE_ROM_APPLICATION=EXPORT,API")

# Define library Kratos which defines the basic python interface
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

# Set batch size in the unity build
if(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosMultiscaleROMCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosMultiscaleROMApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
endif(CMAKE_UNITY_BUILD MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MultiscaleROMApplication.py KratosMultiphysics/MultiscaleROMApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
    PATTERN "*.git" EXCLUDE
    PATTERN "*.c" EXCLUDE
    PATTERN "*.h" EXCLUDE
    PATTERN "*.cpp" EXCLUDE
    PATTERN "*.hpp" EXCLUDE
)

# Install targets
install(TARGETS KratosMultiscaleROMCore DESTINATION libs )
install(TARGETS KratosMultiscaleROMApplication DESTINATION libs )
