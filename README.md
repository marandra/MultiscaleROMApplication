KRATOS Multiphysics Multiscale Reduce Order Application




in the proper locations of kratos/applications/CMakeLists.txt:
message("MULTISCALE_ROM_APPLICATION ......... ${MULTISCALE_ROM_APPLICATION}")

if(${MULTISCALE_ROM_APPLICATION} MATCHES ON)
  add_subdirectory(MultiScaleROMApplication)
endif(${MULTISCALE_ROM_APPLICATION} MATCHES ON)

