// NBody_Comparison.cpp
#include <mpi.h>
#include "Simulation.hpp"
#include "Utils.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

// Entry point of the program
int main(int argc, char** argv) {
    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Determine the size of the world, i.e., how many processes are running
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Determine the rank of the current process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // Ensure the correct number of arguments are passed
    if (argc < 4 || argc > 5) { // Allow for 4 or 5 arguments
        // Only the master process (rank 0) should print the error message
        if (world_rank == 0) {
            std::cerr << "Usage: " << argv[0] << " output_folder min_expansion_factor max_expansion_factor [num_particles]" << std::endl;
        }
        // Abort the MPI execution if the argument count is incorrect
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Parse command line arguments
    std::string output_folder(argv[1]);
    double min_expansion_factor = atof(argv[2]);
    double max_expansion_factor = atof(argv[3]);
    int num_particles = 101 * 101 * 101 * 10; // Set a default value for num_particles
    if (argc == 5) { // If num_particles is provided, overwrite the default value
        num_particles = atoi(argv[4]);
    }

    // Calculate the expansion factor for this process based on its rank
    double expansion_factor = min_expansion_factor + 
                              (max_expansion_factor - min_expansion_factor) / (world_size - 1) * world_rank;

    // Print out the role of the process and the expansion factor it's working with
    if (world_rank == 0) {
        std::cout << "Master process, running simulations with different expansion factors" << std::endl;
    } else {
        std::cout << "Process " << world_rank << " running simulation with expansion factor " << expansion_factor << std::endl;
    }

    // Simulate based on the provided parameter settings
    int total_particles = 101 * 101 * 101 * 10;
    double particle_mass = 1e5 / total_particles;

    // Initialize and run the simulation
    Simulation simulation(1.5, 0.01, 100.0, expansion_factor, 101, particle_mass);
    simulation.initializeParticles(num_particles, 93170929);
    simulation.run(std::nullopt);

    // Retrieve particle positions from the simulation
    std::vector<std::array<double, 3>> positions = simulation.getParticlesPositions();
    // Calculate the correlation function based on the positions
    std::vector<double> correlation = correlationFunction(positions, 100);

    // Send the calculated correlation data to the master process
    if (world_rank != 0) {
        MPI_Send(correlation.data(), correlation.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
        // The master process collects correlation data from all processes
        std::vector<double> all_correlations(correlation.size() * world_size);
        std::copy(correlation.begin(), correlation.end(), all_correlations.begin());

        for (int i = 1; i < world_size; i++) {
            MPI_Recv(&all_correlations[i * correlation.size()], correlation.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Save the collected correlation data to a CSV file
        std::ofstream csv_file(output_folder + "/correlation_functions.csv");
        if (!csv_file.is_open()) {
            std::cerr << "Failed to open file for writing." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        // Write the header of the CSV file
        csv_file << "ExpansionFactor";
        for (int i = 0; i < world_size; i++) {
            csv_file << "," << min_expansion_factor + (max_expansion_factor - min_expansion_factor) / (world_size - 1) * i;
        }
        csv_file << "\n";

        // Write the correlation data for each radius
        for (size_t i = 0; i < correlation.size(); i++) {
            csv_file << "Radius" << i;
            for (int j = 0; j < world_size; j++) {
                csv_file << "," << all_correlations[j * correlation.size() + i];
            }
            csv_file << "\n";
        }
        // Close the file to save the data
        csv_file.close();
    }

    // Finalize the MPI environment and exit the program
    MPI_Finalize();
    return 0;
}

