KRATOS Multiphysics Multiscale Reduce Order Application

To add application to Kratos, add to kratos/CMakeLists.txt
set(MULTISCALE_ROM_APPLICATION "ON")


in the proper locations of kratos/applications/CMakeLists.txt:
message("MULTISCALE_ROM_APPLICATION ......... ${MULTISCALE_ROM_APPLICATION}")

if(${MULTISCALE_ROM_APPLICATION} MATCHES ON)
  add_subdirectory(MultiScaleROMApplication)
endif(${MULTISCALE_ROM_APPLICATION} MATCHES ON)

