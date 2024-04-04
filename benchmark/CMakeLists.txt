find_package(OpenMP REQUIRED)
add_executable(BenchmarkSimulation benchmarks.cpp)
target_link_libraries(BenchmarkSimulation PUBLIC PM_Simulation OpenMP::OpenMP_CXX)

