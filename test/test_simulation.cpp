// test_simulation.cpp
#define CATCH_CONFIG_MAIN  // 该宏指令让Catch自动生成main函数
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <cmath>
#include "Simulation.hpp"  // 确保这是你的Simulation类的正确路径

using namespace Catch; //::Matchers

TEST_CASE("Test potential function for single particle", "[Potential Tests]") {
    // 一颗粒子位于盒子中心，质量为0.01
    double mass = 0.01;
    double time_max = 1.0; // 仅为初始化Simulation对象所用
    double delta_t = 0.1;  // 同上，时间步长
    double width = 100;
    double expansion_factor = 1.0; // 在这个测试中不会用到膨胀因子，但需要为Simulation构造函数提供
    int ncells = 101;
    double particle_mass = 0.01; // 粒子质量

    // 使用初始参数声明Simulation对象
    Simulation sim(time_max, delta_t, width, expansion_factor, ncells, particle_mass);
    sim.addParticle({0.5, 0.5, 0.5}); // 在盒子中心添加一个粒子
    sim.calculateDensity(); // 计算密度
    sim.calculatePotential(); // 计算势能

    double w_c = width / ncells;
    // 沿x轴查看势能函数
    for (int i = 0; i < 101; i++) {
        int j = 50;
        int k = 50;
        if (i == 50 && j == 50 && k == 50) {
            // 忽略中心点，因为点势能在此会是奇异的
            continue;
        } else {
            double pot = sim.getPotentialAtGridIndex(i, j, k); // 获取势能
            // 由于周期性，考虑多源势能的近似值
            double dx1 = std::abs(i - 50) * w_c;
            double dx2 = std::abs(i - 151) * w_c;
            double dx3 = std::abs(i + 51) * w_c;
            double expected_pot =  -mass * (1/dx1 + 1/dx2 + 1/dx3);
            REQUIRE_THAT(pot, Matchers::WithinRel(expected_pot, 0.3)); // 验证势能是否在预期的相对误差范围内
        }
    }
}

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


    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.addParticle({0.5, 0.5, 0.5});
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


    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.addParticle({0.1, 0.1, 0.1});
    sim.addParticle({0.1, 0.1, 0.1}); // Same cell as first particle
    sim.addParticle({0.9, 0.9, 0.9}); // Different cell
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

// ...包含其他测试案例...

TEST_CASE("Gradient of potential is calculated correctly", "[Simulation]") {
    // 初始化测试用的Simulation对象和势能缓冲区
    int nc = 10;
    double box_width = 100.0;
    Simulation sim(1.0, 0.1, box_width, 1.0, nc, 1.0);
    fftw_complex* potential = new fftw_complex[nc * nc * nc]();
    
    // 设置一个简单的测试情况：在中心处势能非零，在其他地方为零
    int centerIndex = (nc / 2) * nc * nc + (nc / 2) * nc + (nc / 2);
    potential[centerIndex][0] = 100.0;  // 只设置实部
    
    // 计算梯度
    auto gradient = sim.calculateGradient(potential);
    
    // 进行测试：中心点周围的梯度值应该是已知的
    REQUIRE(gradient[centerIndex][0] == Approx(0.0).margin(1e-5));
    REQUIRE(gradient[centerIndex][1] == Approx(0.0).margin(1e-5));
    REQUIRE(gradient[centerIndex][2] == Approx(0.0).margin(1e-5));
    
    // 清理
    delete[] potential;
}

TEST_CASE("Gradient calculation with periodic boundaries", "[Simulation]") {
    // 测试边界条件
    int nc = 10;
    double box_width = 100.0;
    Simulation sim(1.0, 0.1, box_width, 1.0, nc, 1.0);
    fftw_complex* potential = new fftw_complex[nc * nc * nc]();
    
    // 设置潜在的线性增加势能
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            for (int k = 0; k < nc; ++k) {
                int index = (k + nc * (j + nc * i));
                potential[index][0] = static_cast<double>(i);
            }
        }
    }
    
    // 计算梯度并检查边缘处是否正确包裹
    auto gradient = sim.calculateGradient(potential);
    
    int edgeIndex = nc - 1;
    int wrappedIndex = sim.wrapIndex(edgeIndex + 1, nc); // 应该为0
    
    REQUIRE(wrappedIndex == 0);
    REQUIRE(gradient[edgeIndex][0] == Approx(-1.0).margin(1e-5));  // 边缘处梯度应该为-1
    
    // 清理
    delete[] potential;
}

TEST_CASE("Particles are updated correctly", "[Simulation]") {
    // 设置初始条件
    Particle particle({0.5, 0.5, 0.5}, {0.0, 0.0, 0.0});
    std::array<double, 3> acceleration = {1.0, 0.0, 0.0}; // 假设沿x轴有一个单位加速度
    double delta_t = 1.0; // 时间步长为1秒

    // 执行更新
    particle.update(delta_t, acceleration);

    // 检查结果
    auto position = particle.getPosition();
    auto velocity = particle.getVelocity();

    REQUIRE(position[0] == Approx(0.5 + 1.0 * delta_t)); // 检查x位置是否正确更新
    REQUIRE(velocity[0] == Approx(1.0 * delta_t)); // 检查x速度是否正确更新
    // 由于沿y和z方向没有加速度，所以它们的位置和速度不应该改变
    REQUIRE(position[1] == Approx(0.5));
    REQUIRE(position[2] == Approx(0.5));
    REQUIRE(velocity[1] == Approx(0.0));
    REQUIRE(velocity[2] == Approx(0.0));
}

// 可以继续添加其他测试案例...
