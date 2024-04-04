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
    std::optional<std::string> output_folder;
    int nc = 10, np_per_cell = 10, random_seed = 93170929;
    double dt = 0.1, expansion_factor = 1.0, total_time = 1.0;
    std::string output;

    
    // Parsing command-line parameters
    if (argc == 2 && std::string(argv[1]) == "-h") {
        std::cout << "Help message and usage instructions" << std::endl;
        return 0;
    }


    for (int i = 1; i < argc; i += 2) {
        std::string arg(argv[i]);
        if (arg == "-nc") {
            nc = std::atoi(argv[i + 1]);
        } else if (arg == "-np") {
            np_per_cell = std::atoi(argv[i + 1]);
        } else if (arg == "-t") {
            total_time = std::atof(argv[i + 1]);
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

    // Simulate based on the provided parameter settings
    int total_particles = nc * nc * nc * np_per_cell;
    double box_width = 100;
    double particle_mass = 1e5 / total_particles;
    
    // Initialize simulation
    Simulation sim(total_time, dt, box_width, expansion_factor, nc, particle_mass);
    sim.initializeParticles(total_particles, random_seed);

    // Run the simulation and save the density distribution image every 10 time steps
    sim.run(output_folder);

    return 0;
}
