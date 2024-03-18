add_library(PM_Simulation STATIC Simulation.cpp Utils.cpp)
target_include_directories(PM_Simulation PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(PM_Simulation PUBLIC fftw3 OpenMP::OpenMP_CXX)