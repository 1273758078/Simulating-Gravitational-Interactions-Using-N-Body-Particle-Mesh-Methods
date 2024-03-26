// Simulation.cpp
#include "Simulation.hpp"
#include "Utils.hpp"
#include <fftw3.h>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <algorithm>

// Particle 类的实现
Particle::Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity)
    : position_(position), velocity_(velocity) {}

// 获取粒子的位置
std::array<double, 3> Particle::getPosition() const {
    return position_;
}

// 获取粒子的速度
std::array<double, 3> Particle::getVelocity() const {
    return velocity_;
}

// 设置粒子的位置
void Particle::setPosition(const std::array<double, 3>& position) {
    position_ = position;
}

// 根据加速度和时间步长更新粒子的位置和速度
void Particle::update(double delta_t, const std::array<double, 3>& acceleration) {
    for (int i = 0; i < 3; ++i) {
        velocity_[i] += acceleration[i] * delta_t;
        position_[i] += velocity_[i] * delta_t;
    }
}

// Simulation 类的实现
Simulation::Simulation(double time_max, double delta_t, double box_width, double expansion_factor, int nc, double particle_mass)
    : time_max_(time_max), delta_t_(delta_t), box_width_(box_width), expansion_factor_(expansion_factor), nc_(nc), particle_mass_(particle_mass) {
    // 分配内存给密度和势能缓冲区
    density_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
    potential_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
    // 初始化密度缓冲区
    initializeDensityBuffer();
}

// 析构函数，释放分配的内存
Simulation::~Simulation() {
    if (density_buffer_ != nullptr) {
        fftw_free(density_buffer_);
    }
    if (potential_buffer_ != nullptr) {
        fftw_free(potential_buffer_);
    }
}

// 向模拟中添加粒子
void Simulation::addParticle(const std::array<double, 3>& position) {
    particles_.emplace_back(position, std::array<double, 3>{0.0, 0.0, 0.0});
}

// 初始化粒子，根据给定的数量和种子，随机分布粒子
void Simulation::initializeParticles(int num_particles, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0.0, box_width_); // 确保粒子位置在箱子范围内
    for (int i = 0; i < num_particles; ++i) {
        addParticle({dis(gen), dis(gen), dis(gen)});
    }
}

// 运行模拟
void Simulation::run() {
    for (double t = 0; t < time_max_; t += delta_t_) {
        for (auto& particle : particles_) {
            std::array<double, 3> acceleration = {0, 0, 0}; // 加速度计算的占位符
            particle.update(delta_t_, acceleration);
        }
        box_width_ *= expansion_factor_; // 模拟箱子的膨胀
    }
}

// 初始化密度缓冲区
void Simulation::initializeDensityBuffer() {
    std::fill_n(density_buffer_, nc_ * nc_ * nc_, fftw_complex{0.0, 0.0});
}

// 根据粒子的位置计算密度
void Simulation::calculateDensity() {
    initializeDensityBuffer();
    double cell_volume = std::pow(box_width_ / nc_, 3);
    for (const auto& particle : particles_) {
        std::array<double, 3> position = particle.getPosition();
        int i = static_cast<int>(position[0] * nc_ / box_width_);
        int j = static_cast<int>(position[1] * nc_ / box_width_);
        int k = static_cast<int>(position[2] * nc_ / box_width_);
        int index = (k + nc_ * (j + nc_ * i)) % (nc_ * nc_ * nc_);
        density_buffer_[index][0] += particle_mass_ / cell_volume;
    }
}

// 计算势能
void Simulation::calculatePotential() {
    // 创建并执行正向FFT计划，将密度数据转换到频率空间
    fftw_plan forward_plan = fftw_plan_dft_3d(nc_, nc_, nc_, density_buffer_, potential_buffer_, FFTW_FORWARD, FFTW_MEASURE);
    fftw_execute(forward_plan);
    // 完成后立即释放FFT计划资源
    fftw_destroy_plan(forward_plan);

    // 在频率空间处理势能数据
    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            for (int k = 0; k < nc_; ++k) {
                int index = k + nc_ * (j + nc_ * i);
                // 计算每个点的k向量
                std::array<double, 3> k_vec = getKVector(i, j, k);
                // 计算k向量的平方
                double k_squared = k_vec[0] * k_vec[0] + k_vec[1] * k_vec[1] + k_vec[2] * k_vec[2];
                // 避免除以零
                if (k_squared == 0) {
                    potential_buffer_[index][0] = 0;
                    potential_buffer_[index][1] = 0;
                } else {
                    // 根据k向量的平方调整势能值
                    double scaling_factor = -4 * M_PI / k_squared;
                    potential_buffer_[index][0] *= scaling_factor;
                    potential_buffer_[index][1] *= scaling_factor;
                }
            }
        }
    }

    // 创建并执行逆向FFT计划，将势能数据从频率空间转换回实空间
    fftw_plan inverse_plan = fftw_plan_dft_3d(nc_, nc_, nc_, potential_buffer_, density_buffer_, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_execute(inverse_plan);
    // 完成后立即释放FFT计划资源
    fftw_destroy_plan(inverse_plan);

    // 应用归一化因子，以确保逆FFT后的结果正确
    double norm_factor = 1.0 / 8 * (nc_ * nc_ * nc_);
    for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
        density_buffer_[i][0] *= norm_factor;
        // 对于势能，我们只关心实部
        density_buffer_[i][1] = 0;
    }
}

// 根据网格索引计算对应的k向量，考虑Nyquist频率和频率折叠
std::array<double, 3> Simulation::getKVector(int i, int j, int k) {
    double kx = (i < nc_ / 2) ? i : (i - nc_);
    double ky = (j < nc_ / 2) ? j : (j - nc_);
    double kz = (k < nc_ / 2) ? k : (k - nc_);
    return {2 * M_PI * kx / box_width_, 2 * M_PI * ky / box_width_, 2 * M_PI * kz / box_width_};
}

// 在给定网格索引处获取势能值
double Simulation::getPotentialAtGridIndex(int i, int j, int k) {
    int index = k + nc_ * (j + nc_ * i);
    // 注意：这里使用density_buffer_来保存逆FFT后的结果，实际上应该是势能值
    return density_buffer_[index][0]; // 我们只关心实部
}

std::vector<std::array<double, 3>> Simulation::calculateGradient(const fftw_complex* potential) {
    std::vector<std::array<double, 3>> gradient(nc_ * nc_ * nc_);
    
    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            for (int k = 0; k < nc_; ++k) {
                int index = (k + nc_ * (j + nc_ * i));
                gradient[index][0] = -(potential[wrapIndex(i+1, nc_) * nc_ * nc_ + j * nc_ + k][0] - potential[wrapIndex(i-1, nc_) * nc_ * nc_ + j * nc_ + k][0]) / (2 * box_width_ / nc_);
                gradient[index][1] = -(potential[i * nc_ * nc_ + wrapIndex(j+1, nc_) * nc_ + k][0] - potential[i * nc_ * nc_ + wrapIndex(j-1, nc_) * nc_ + k][0]) / (2 * box_width_ / nc_);
                gradient[index][2] = -(potential[i * nc_ * nc_ + j * nc_ + wrapIndex(k+1, nc_)][0] - potential[i * nc_ * nc_ + j * nc_ + wrapIndex(k-1, nc_)][0]) / (2 * box_width_ / nc_);
            }
        }
    }
    return gradient;
}

int Simulation::wrapIndex(int index, int max) {
    return (index + max) % max; // 模运算
}
