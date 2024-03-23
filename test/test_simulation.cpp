#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <cmath>
#include "Simulation.hpp"

using namespace Catch::Matchers;

/**
 * @brief Fill out this test function by filling in the TODOs
 * Tests the calculation of the gravitational potential due to a single particle
 * Potential is examined along the x axis throught the centre of the cube
 */
TEST_CASE("Test potential function for single particle", "Potential Tests")
{
    //One particle in the centre of the box with mass 0.01
    //TODO: Set up single particle with position {0.5, 0.5, 0.5}
    double mass = 0.01;

    // These properties should apply to your box
    double width = 100;
    int ncells = 101;

    //Declare simulation with particle setup, width, and number of cells. 
    //TODO: Simulation sim(...);
    //TODO: calculate density
    //TODO: calculate potential

    double w_c = width / ncells;
    // Look at the potential function along the x axis
    for (int i = 0; i < 101; i++)
    {
        int j = 50;
        int k = 50;
        if (i == 50 && j == 50 && k == 50)
        {
            //ignore the centre where point potential would be singular
            continue;
        }
        else
        {
            // ignore imaginary component
            // approximate potential for multiple sources due to periodicity
            double dx1 = std::abs(i - 50) * w_c;
            double dx2 = std::abs(i - 151) * w_c;
            double dx3 = std::abs(i + 51) * w_c;
            double pot; //TODO = your potential function at indices (i,j,k)
            double expected_pot =  -mass *(1/ dx1 + 1/dx2 + 1/dx3);
            double diff = pot - expected_pot;
            REQUIRE_THAT(pot, WithinRel(expected_pot, 0.3));
        }
    }
}


//////////////////////////////////////////////////////////////////

#define CATCH_CONFIG_MAIN  // 该宏指令会让Catch自动生成main函数
#include "catch.hpp"
#include "Simulation.hpp"  // 你的Simulation类的头文件

TEST_CASE("Imaginary components are zero", "[density]") {
    // 设置模拟的初始参数
    const double t_max = 10.0;  // 模拟的最大时间
    const double delta_t = 0.1;  // 时间步长
    const double box_width = 100.0; // 盒子宽度
    const double expansion_factor = 1.0; // 扩展因子
    const int nc = 10; // 单边网格数
    const double particle_mass = 1.0; // 粒子质量

    // 创建一个不包含粒子的模拟实例
    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.initializeDensityBuffer(); // 初始化密度缓冲区
    sim.calculateDensity(); // 计算密度

    for (int i = 0; i < sim.getTotalCells(); ++i) {
        REQUIRE(sim.getDensityBuffer()[i][1] == 0.0);
    }
}

TEST_CASE("Density is zero without particles", "[density]") {
    // 设置模拟的初始参数
    const double t_max = 10.0;  // 模拟的最大时间
    const double delta_t = 0.1;  // 时间步长
    const double box_width = 100.0; // 盒子宽度
    const double expansion_factor = 1.0; // 扩展因子
    const int nc = 10; // 单边网格数
    const double particle_mass = 1.0; // 粒子质量

    // 创建一个不包含粒子的模拟实例
    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.initializeDensityBuffer(); // 初始化密度缓冲区
    sim.calculateDensity(); // 计算密度

    for (int i = 0; i < sim.getTotalCells(); ++i) {
        REQUIRE(sim.getDensityBuffer()[i][0] == 0.0);
    }
}

TEST_CASE("Density with a single particle", "[density]") {
    // 设置模拟的初始参数
    const double t_max = 10.0;  // 模拟的最大时间
    const double delta_t = 0.1;  // 时间步长
    const double box_width = 100.0; // 盒子宽度
    const double expansion_factor = 1.0; // 扩展因子
    const int nc = 10; // 单边网格数
    const double particle_mass = 1.0; // 粒子质量

    std::vector<Particle> particles = {Particle({0.5, 0.5, 0.5}, {0.0, 0.0, 0.0})};
    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass, particles);
    sim.initializeDensityBuffer();
    sim.calculateDensity();

    int expected_index = sim.getCellIndex(0.5, 0.5, 0.5);
    REQUIRE(sim.getDensityBuffer()[expected_index][0] == sim.getParticleMass() / sim.getCellVolume());

    for (int i = 0; i < sim.getTotalCells(); ++i) {
        if (i != expected_index) {
            REQUIRE(sim.getDensityBuffer()[i][0] == 0.0);
        }
    }
}

TEST_CASE("Density with multiple particles", "[density]") {
    // 设置模拟的初始参数
    const double t_max = 10.0;  // 模拟的最大时间
    const double delta_t = 0.1;  // 时间步长
    const double box_width = 100.0; // 盒子宽度
    const double expansion_factor = 1.0; // 扩展因子
    const int nc = 10; // 单边网格数
    const double particle_mass = 1.0; // 粒子质量

    std::vector<Particle> particles = {
        Particle({0.1, 0.1, 0.1}, {0.0, 0.0, 0.0}),
        Particle({0.1, 0.1, 0.1}, {0.0, 0.0, 0.0}),  // Same cell as first particle
        Particle({0.9, 0.9, 0.9}, {0.0, 0.0, 0.0})   // Different cell
    };
    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass, particles);
    sim.initializeDensityBuffer();
    sim.calculateDensity();

    int cellIndex1 = sim.getCellIndex(0.1, 0.1, 0.1);
    int cellIndex2 = sim.getCellIndex(0.9, 0.9, 0.9);

    // Check the density in the cell with two particles
    REQUIRE(sim.getDensityBuffer()[cellIndex1][0] == 2 * sim.getParticleMass() / sim.getCellVolume());

    // Check the density in the cell with one particle
    REQUIRE(sim.getDensityBuffer()[cellIndex2][0] == sim.getParticleMass() / sim.getCellVolume());

    // Check that all other cells have zero density
    for (int i = 0; i < sim.getTotalCells(); ++i) {
        if (i != cellIndex1 && i != cellIndex2) {
            REQUIRE(sim.getDensityBuffer()[i][0] == 0.0);
        }
    }
}

// ... 其他测试案例 ...
