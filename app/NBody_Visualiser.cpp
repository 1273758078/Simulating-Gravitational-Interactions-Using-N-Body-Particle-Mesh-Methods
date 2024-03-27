// NBody_Visualiser.cpp
#include "Simulation.hpp"
#include "Utils.hpp"
#include <optional>
#include <iostream>
#include <string>
#include <cstdlib>
#include <filesystem>

int main(int argc, char* argv[])
{
    // 解析命令行参数
    if (argc == 2 && std::string(argv[1]) == "-h") {
        std::cout << "Help message and usage instructions" << std::endl;
        return 0;
    }

    std::optional<std::string> output_folder;
    int nc, np_per_cell, total_time, random_seed;
    double dt, expansion_factor;
    std::string output;

    for (int i = 1; i < argc; i += 2) {
        std::string arg(argv[i]);
        if (arg == "-nc") {
            nc = std::atoi(argv[i + 1]);
        } else if (arg == "-np") {
            np_per_cell = std::atoi(argv[i + 1]);
        } else if (arg == "-t") {
            total_time = std::atoi(argv[i + 1]);
        } else if (arg == "-dt") {
            dt = std::atof(argv[i + 1]);
        } else if (arg == "-F") {
            expansion_factor = std::atof(argv[i + 1]);
        } else if (arg == "-o") {
            output_folder = std::string(argv[i + 1]);
            if (!std::filesystem::exists(output_folder.value())) {
                std::filesystem::create_directory(output_folder.value());
            }
        } else if (arg == "-s") {
            random_seed = std::atoi(argv[i + 1]);
        }
    }

    // 根据提供的参数设置模拟
    int total_particles = nc * nc * nc * np_per_cell;
    double box_width = 100;
    double particle_mass = 1e5 / total_particles;
    
    // 初始化模拟
    Simulation sim(total_time, dt, box_width, expansion_factor, nc, particle_mass);
    sim.initializeParticles(total_particles, random_seed);

    // 运行模拟，并每10个时间步长保存一次密度分布图像
    sim.run(output_folder, 10);

    return 0;
}
