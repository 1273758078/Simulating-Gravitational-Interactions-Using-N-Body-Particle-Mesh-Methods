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

int main()
{

    return 0;
}