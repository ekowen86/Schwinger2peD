set (SCHWINGER_2PED_OBJS utils.cpp hmc.cpp inverters.cpp dirac_op.cpp measurements.cpp blas.cpp io.cpp)

add_library(schwinger2peD STATIC ${SCHWINGER_2PED_OBJS})
target_include_directories(schwinger2peD PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include/)

if(ENABLE_OPENMP)
  target_link_libraries(schwinger2peD PUBLIC OpenMP::OpenMP_CXX)
endif()
