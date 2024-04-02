// Simulation.hpp
#pragma once
#include <fftw3.h>
#ifndef SIMULATION_HPP
#define SIMULATION_HPP
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include <optional>
#include <string>

// 粒子类
class Particle {
public:
    Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity);
    std::array<double, 3> getPosition() const;
    std::array<double, 3> getVelocity() const;
    void setPosition(const std::array<double, 3>& position);
    void setVelocity(const std::array<double, 3>& velocity);
    void update(double delta_t, const std::array<double, 3>& acceleration);


private:
    std::array<double, 3> position_; // 粒子位置
    std::array<double, 3> velocity_; // 粒子速度
};

// 模拟类
class Simulation {
public:
    Simulation(double time_max, double delta_t, double box_width, double expansion_factor, int nc, double particle_mass);
    ~Simulation();
    void addParticle(const std::array<double, 3>& position);
    void initializeParticles(int num_particles, unsigned seed);
    void run(const std::optional<std::string>& output_folder);
    void initializeDensityBuffer();
    void calculateDensity();
    //得到某一个cell的索引
    int getCellIndex(double x, double y, double z);
    void calculatePotential();
    // 获取给定网格索引(i, j, k)处的势能
    double getPotentialAtGridIndex(int i, int j, int k);
    // 梯度计算函数
    std::vector<std::array<double, 3>> calculateGradient(const fftw_complex* potential);
    // 更新所有粒子的函数
    void updateParticles(const std::vector<std::array<double, 3>>& gradients, double delta_t);
    // 新增扩展盒子的函数
    void expandBox(double expansion_factor);
    // 得到单元格的总数
    int getTotalCells() const;
    // 得到密度函数的数据
    fftw_complex* getDensityBuffer() const;
    // 将两个边界接在一起
    int wrapIndex(int index, int max);
    // 得到粒子的质量
    double getParticleMass() const;
    // 得到单元格的体积
    double getCellVolume() const;
    // 提取所有粒子的位置，返回一个向量
    std::vector<std::array<double, 3>> getParticlesPositions() const;

private:
    std::vector<Particle> particles_;  // 模拟中的粒子数组
    double time_max_;                  // 模拟的最大时间
    double delta_t_;                   // 模拟的时间步长
    double box_width_;                 // 模拟的空间盒子宽度
    double expansion_factor_;          // 空间膨胀因子
    int nc_;                           // 网格单元数
    double particle_mass_;             // 粒子质量
    fftw_complex* density_buffer_;     // 密度缓冲区
    fftw_complex* k_buffer_; // 过渡区
    fftw_complex* potential_buffer_;   // 势能缓冲区
    double cell_width_; // 单元格的宽度，如果有的话
};

#endif // SIMULATION_HPP
