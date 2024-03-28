add_executable(NBody_Visualiser NBody_Visualiser.cpp)
target_link_libraries(NBody_Visualiser PUBLIC PM_Simulation)
target_compile_options(NBody_Visualiser PRIVATE -O2)


add_executable(NBody_Comparison NBody_Comparison.cpp)
target_link_libraries(NBody_Comparison PUBLIC PM_Simulation MPI::MPI_CXX)