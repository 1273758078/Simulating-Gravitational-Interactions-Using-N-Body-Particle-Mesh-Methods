# add_executable(BenchmarkSimulation benchmarks.cpp)
# target_link_libraries(BenchmarkSimulation PUBLIC PM_Simulation)

# cmake_minimum_required(VERSION 3.1)
# project(BenchmarkProject) # 更改为你的项目名

# 寻找OpenMP包
find_package(OpenMP REQUIRED)

# 创建名为BenchmarkSimulation的可执行文件
add_executable(BenchmarkSimulation benchmarks.cpp)

# 假设PM_Simulation是你项目中的一个库，我们将它链接到BenchmarkSimulation可执行文件
# 同时确保OpenMP也被链接

target_link_libraries(BenchmarkSimulation PUBLIC PM_Simulation OpenMP::OpenMP_CXX)

