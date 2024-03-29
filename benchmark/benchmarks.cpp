#include <iostream>
#include <string>
#include <omp.h>
#include <chrono>
#include "Simulation.hpp"

class BenchmarkData
{
public:
    BenchmarkData(std::string benchmark_name, int threads) : name(benchmark_name), num_threads(threads) {}
    double time;
    std::string name;
    int num_threads;
    //use for additional context such as number of particles / cells 
    std::string info;

    void start()
    {
        t1 = std::chrono::high_resolution_clock::now();
    }

    void finish()
    {
        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        time = (t2 - t1).count() / 1e9;
    }

    std::chrono::high_resolution_clock::time_point t1;
};

std::ostream& operator<<(std::ostream &os, const BenchmarkData& b)
{
    std::cout << "Benchmarking " << b.name << " with " << b.num_threads << " threads." << std::endl;
    std::cout << "Time = " << b.time << std::endl;
    std::cout << "Info: " << b.info << std::endl;
    return os;
}

int main() {
    // 初始化模拟参数
    double time_max = 10.0, delta_t = 0.01, box_width = 100.0, expansion_factor = 1.0;
    int nc = 64, num_particles = 1024;
    double particle_mass = 1.0;
    std::optional<std::string> output_folder = std::nullopt;

    // 创建模拟实例
    Simulation sim(time_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.initializeParticles(num_particles, 42); // 使用随机种子初始化粒子

    // 线程数数组
    int thread_nums[] = {1, 2, 4, 8};

    for (int num_threads : thread_nums) {
        omp_set_num_threads(num_threads); // 设置OpenMP的线程数

        // 重置模拟状态，以保证每次测试的一致性
        sim.initializeParticles(num_particles, 42);

        // calculateDensity 基准测试
        BenchmarkData bdDensity("calculateDensity", num_threads);
        bdDensity.start();
        sim.calculateDensity();
        bdDensity.finish();
        std::cout << bdDensity << std::endl;

        // calculatePotential 基准测试
        BenchmarkData bdPotential("calculatePotential", num_threads);
        bdPotential.start();
        sim.calculatePotential();
        bdPotential.finish();
        std::cout << bdPotential << std::endl;

        // calculateGradient 基准测试
        BenchmarkData bdGradient("calculateGradient", num_threads);
        bdGradient.start();
        auto gradients = sim.calculateGradient(sim.getDensityBuffer()); // 假设有获取密度缓冲区的函数
        bdGradient.finish();
        std::cout << bdGradient << std::endl;

        // updateParticles 基准测试
        BenchmarkData bdUpdate("updateParticles", num_threads);
        bdUpdate.start();
        sim.updateParticles(gradients, delta_t);
        bdUpdate.finish();
        std::cout << bdUpdate << std::endl;
    }

    return 0;
}