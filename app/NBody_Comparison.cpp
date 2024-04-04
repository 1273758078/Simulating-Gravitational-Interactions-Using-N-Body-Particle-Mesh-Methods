// // NBody_Comparison.cpp
// #include <mpi.h>
// #include "Simulation.hpp"
// #include "Utils.hpp"
// #include <iostream>
// #include <vector>
// #include <fstream>
// #include <sstream>
// #include <string>

// int main(int argc, char** argv) {
//     MPI_Init(&argc, &argv);

//     int world_size;
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);

//     int world_rank;
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//     if (argc != 5) {
//         if (world_rank == 0) {
//             std::cerr << "Usage: " << argv[0] << " output_folder min_expansion_factor max_expansion_factor num_particles" << std::endl;
//         }
//         MPI_Abort(MPI_COMM_WORLD, 1);
//     }

//     std::string output_folder(argv[1]);
//     double min_expansion_factor = atof(argv[2]);
//     double max_expansion_factor = atof(argv[3]);
//     int num_particles = atoi(argv[4]);

//     // 分布式运算，每个进程计算一个扩张因子
//     double expansion_factor = min_expansion_factor + 
//                               (max_expansion_factor - min_expansion_factor) / (world_size - 1) * world_rank;

//     // 进程0是主进程
//     if (world_rank == 0) {
//         std::cout << "Master process, running simulations with different expansion factors" << std::endl;
//     } else {
//         std::cout << "Process " << world_rank << " running simulation with expansion factor " << expansion_factor << std::endl;
//     }

//     // 初始化并运行模拟
//     Simulation simulation(100.0, 0.1, 100.0, expansion_factor, 100, 1e5);
//     simulation.initializeParticles(num_particles, 12345);  // Make sure to use the same seed for reproducibility
//     simulation.run(std::nullopt);

//     // 计算径向相关函数
//     auto positions = simulation.getParticlesPositions();
//     std::vector<double> correlation = correlationFunction(positions, 100);

//     // 非主进程发送其数据到主进程
//     if (world_rank != 0) {
//         MPI_Send(correlation.data(), correlation.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
//     } else {
//         std::vector<double> all_correlations(correlation.size() * world_size);
//         std::copy(correlation.begin(), correlation.end(), all_correlations.begin());

//         // 主进程收集数据
//         for (int i = 1; i < world_size; i++) {
//             MPI_Recv(&all_correlations[i * correlation.size()], correlation.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         }

//         // 将所有径向相关函数输出到csv文件
//         std::ofstream csv_file(output_folder + "/correlation_functions.csv");
//         csv_file << "ExpansionFactor";
//         for (int i = 0; i < world_size; i++) {
//             csv_file << "," << min_expansion_factor + (max_expansion_factor - min_expansion_factor) / (world_size - 1) * i;
//         }
//         csv_file << "\n";

//         for (size_t i = 0; i < correlation.size(); i++) {
//             csv_file << "Radius" << i;
//             for (int j = 0; j < world_size; j++) {
//                 csv_file << "," << all_correlations[j * correlation.size() + i];
//             }
//             csv_file << "\n";
//         }
//         csv_file.close();
//     }

//     MPI_Finalize();
//     return 0;
// }


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

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (argc != 5) {
        if (world_rank == 0) {
            std::cerr << "Usage: " << argv[0] << " output_folder min_expansion_factor max_expansion_factor num_particles" << std::endl;
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    std::string output_folder(argv[1]);
    double min_expansion_factor = atof(argv[2]);
    double max_expansion_factor = atof(argv[3]);
    int num_particles = atoi(argv[4]);

    double expansion_factor = min_expansion_factor + 
                              (max_expansion_factor - min_expansion_factor) / (world_size - 1) * world_rank;

    if (world_rank == 0) {
        std::cout << "Master process, running simulations with different expansion factors" << std::endl;
    } else {
        std::cout << "Process " << world_rank << " running simulation with expansion factor " << expansion_factor << std::endl;
    }

    Simulation simulation(1.0, 0.01, 100.0, expansion_factor, 10, 100);
    // Simulation simulation(1.0, 0.01, 100.0, expansion_factor, 100, 1e5);
    simulation.initializeParticles(num_particles, 93170929);
    simulation.run(std::nullopt);

    std::vector<std::array<double, 3>> positions = simulation.getParticlesPositions();
    std::vector<double> correlation = correlationFunction(positions, 100);

    // // 打印前五个粒子的位置（如果它们存在）
    // for (size_t i = 0; i < positions.size() && i < 5; ++i) {
    //     std::cout << "Particle " << i << ": (" 
    //             << positions[i][0] << ", " 
    //             << positions[i][1] << ", " 
    //             << positions[i][2] << ")" << std::endl;
    // }

    if (world_rank != 0) {
        MPI_Send(correlation.data(), correlation.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    } else {
        std::vector<double> all_correlations(correlation.size() * world_size);
        std::copy(correlation.begin(), correlation.end(), all_correlations.begin());

        for (int i = 1; i < world_size; i++) {
            MPI_Recv(&all_correlations[i * correlation.size()], correlation.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        std::ofstream csv_file(output_folder + "/correlation_functions.csv");
        if (!csv_file.is_open()) {
            std::cerr << "Failed to open file for writing." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        csv_file << "ExpansionFactor";
        for (int i = 0; i < world_size; i++) {
            csv_file << "," << min_expansion_factor + (max_expansion_factor - min_expansion_factor) / (world_size - 1) * i;
        }
        csv_file << "\n";

        for (size_t i = 0; i < correlation.size(); i++) {
            csv_file << "Radius" << i;
            for (int j = 0; j < world_size; j++) {
                csv_file << "," << all_correlations[j * correlation.size() + i];
            }
            csv_file << "\n";
        }
        csv_file.close();
    }

    MPI_Finalize();
    return 0;
}
