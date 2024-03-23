#pragma once

#include <fftw3.h>

// class Simulation
// {
// public:

//     /**
//      * @brief Run a particle mesh simulation from t=0 to t_max
//      * 
//      */
//     void run();

// };


#include <vector>
#include <array>
#include <random>
#include <cmath>

// 粒子类定义
class Particle {
public:
    // 构造函数用于初始化位置和速度
    Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity)
        : position_(position), velocity_(velocity) {}

    // 获取位置和速度的函数
    std::array<double, 3> getPosition() const { return position_; }
    std::array<double, 3> getVelocity() const { return velocity_; }

    // 设置位置的函数
    void setPosition(const std::array<double, 3>& position) { position_ = position; }

    // 更新粒子位置和速度的方法
    void update(double delta_t, const std::array<double, 3>& acceleration) {
        for (int i = 0; i < 3; ++i) {
            velocity_[i] += acceleration[i] * delta_t;
            position_[i] += velocity_[i] * delta_t;
        }
    }

private:
    std::array<double, 3> position_; // 位置
    std::array<double, 3> velocity_; // 速度
};

// 模拟类定义
class Simulation {
public:
    // 构造函数
    Simulation(double time_max, double delta_t, double box_width, double expansion_factor, int nc, double particle_mass)
        : time_max_(time_max), delta_t_(delta_t), box_width_(box_width), expansion_factor_(expansion_factor, nc_(nc), particle_mass_(particle_mass)) {}


    // 添加单个粒子到自定义位置的方法
    void addParticle(const std::array<double, 3>& position) {
        std::array<double, 3> velocity = {0, 0, 0}; // 假设初始速度为零
        particles_.emplace_back(position, velocity);
    }

    // 生成随机分布的粒子的方法
    void initializeParticles(int num_particles, unsigned seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> dis(0.0, 1.0);
        
        for (int i = 0; i < num_particles; ++i) {
            std::array<double, 3> position = {dis(gen), dis(gen), dis(gen)};
            addParticle(position);
        }
    }

    // 运行模拟的方法
    void run() {
        for (double t = 0; t < time_max_; t += delta_t_) {
            for (auto& particle : particles_) {
                std::array<double, 3> acceleration = {0, 0, 0}; // 占位符，实际应计算加速度
                particle.update(delta_t_, acceleration);
            }
            box_width_ *= expansion_factor_; // 扩展盒子
        }
    }

    fftw_complex* density_buffer_;

    // 初始化密度缓冲区的函数
    void initializeDensityBuffer() {
        int n_total_cells = nc_ * nc_ * nc_;
        density_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n_total_cells);
        for (int i = 0; i < n_total_cells; ++i) {
            density_buffer_[i][0] = 0.0; // 实部设为0
            density_buffer_[i][1] = 0.0; // 虚部设为0（因为密度是实函数）
        }
    }

    // 计算密度函数的函数
    void calculateDensity() {
        double cell_width = box_width_ / nc_; // 网格单元的宽度
        double cell_volume = std::pow(cell_width, 3);  // 网格单元的体积

        // 清除之前的密度值
        int n_total_cells = nc_ * nc_ * nc_;
        for (int i = 0; i < n_total_cells; ++i) {
            density_buffer_[i][0] = 0.0; // 只需要清理实部
        }

        // 对每个粒子进行循环，更新它们在密度缓冲区中的贡献
        for (const auto& particle : particles_) {
            std::array<double, 3> position = particle.getPosition();
            int i = static_cast<int>(position[0] * nc_);
            int j = static_cast<int>(position[1] * nc_);
            int k = static_cast<int>(position[2] * nc_);
            int index = k + nc_ * (j + nc_ * i);
            density_buffer_[index][0] += particle_mass_ / cell_volume; // 更新密度值
        }
    }

    // 析构函数来释放密度缓冲区
    ~Simulation() {
        if (density_buffer_ != nullptr) {
            fftw_free(density_buffer_);
        }
    }


private:
    std::vector<Particle> particles_; // 粒子数组
    double time_max_; // 模拟的最大时间
    double delta_t_; // 时间步长
    double box_width_; // 盒子宽度
    double expansion_factor_; // 扩展因子
    int nc_; //number of cells wide
    double particle_mass_; //粒子质量
};
