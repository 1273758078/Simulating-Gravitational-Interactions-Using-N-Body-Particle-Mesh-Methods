#include "Simulation.hpp"
#include "Utils.hpp"

// void Simulation::run()
// {

// }

#include <fftw3.h>
#include <vector>
#include <array>
#include <random>
#include <cmath>

// 粒子类
class Particle {
public:
    // 构造函数：初始化位置和速度
    Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity)
        : position_(position), velocity_(velocity) {}

    // 获取粒子位置
    std::array<double, 3> getPosition() const { return position_; }
    // 获取粒子速度
    std::array<double, 3> getVelocity() const { return velocity_; }

    // 设置粒子位置
    void setPosition(const std::array<double, 3>& position) { position_ = position; }

    // 根据给定的加速度和时间步长更新粒子的位置和速度
    void update(double delta_t, const std::array<double, 3>& acceleration) {
        for (int i = 0; i < 3; ++i) {
            velocity_[i] += acceleration[i] * delta_t;
            position_[i] += velocity_[i] * delta_t;
        }
    }

private:
    std::array<double, 3> position_; // 粒子位置
    std::array<double, 3> velocity_; // 粒子速度
};

// 模拟类
class Simulation {
public:
    // 构造函数：设置模拟的参数
    Simulation(double time_max, double delta_t, double box_width, double expansion_factor, int nc, double particle_mass)
        : time_max_(time_max), delta_t_(delta_t), box_width_(box_width), expansion_factor_(expansion_factor), nc_(nc), particle_mass_(particle_mass) {
        // 分配内存给密度缓冲区和势能缓冲区
        density_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
        potential_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
        // 初始化密度缓冲区
        initializeDensityBuffer();
    }

    // 析构函数：释放分配的内存
    ~Simulation() {
        if (density_buffer_ != nullptr) {
            fftw_free(density_buffer_);
        }
        if (potential_buffer_ != nullptr) {
            fftw_free(potential_buffer_);
        }
    }

    // 添加粒子到模拟中
    void addParticle(const std::array<double, 3>& position) {
        particles_.emplace_back(position, std::array<double, 3>{0.0, 0.0, 0.0});
    }

    // 初始化粒子，创建随机分布的粒子
    void initializeParticles(int num_particles, unsigned seed) {
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < num_particles; ++i) {
            addParticle({dis(gen), dis(gen), dis(gen)});
        }
    }

    // 运行模拟：更新所有粒子的位置和速度
    void run() {
        for (double t = 0; t < time_max_; t += delta_t_) {
            for (auto& particle : particles_) {
                std::array<double, 3> acceleration = {0, 0, 0}; // 加速度计算的占位符
                particle.update(delta_t_, acceleration);
            }
            box_width_ *= expansion_factor_; // 模拟箱子的膨胀
        }
    }

    // 初始化密度缓冲区，将所有值设置为0
    void initializeDensityBuffer() {
        std::fill_n(density_buffer_, nc_ * nc_ * nc_, fftw_complex{0.0, 0.0});
    }

    // 根据粒子的位置计算密度
    void calculateDensity() {
        initializeDensityBuffer(); // 首先重置密度缓冲区
        double cell_volume = std::pow(box_width_ / nc_, 3); // 计算每个网格单元的体积
        for (const auto& particle : particles_) {
            std::array<double, 3> position = particle.getPosition();
            int i = static_cast<int>(position[0] * nc_);
            int j = static_cast<int>(position[1] * nc_);
            int k = static_cast<int>(position[2] * nc_);
            int index = (k + nc_ * (j + nc_ * i));
            // 更新密度值，密度等于粒子质量除以单元体积
            density_buffer_[index][0] += particle_mass_ / cell_volume;
        }
    }

    // 根据密度缓冲区计算势能
    void calculatePotential() {
        // 创建FFT的正向计划
        fftw_plan forward_plan = fftw_plan_dft_3d(nc_, nc_, nc_, density_buffer_, potential_buffer_, FFTW_FORWARD, FFTW_MEASURE);
        // 执行正向FFT变换
        fftw_execute(forward_plan);
        // 释放正向计划资源
        fftw_destroy_plan(forward_plan);
        
        // 在频率空间处理势能数据
        for (int i = 0; i < nc_; ++i) {
            for (int j = 0; j < nc_; ++j) {
                for (int k = 0; k < nc_; ++k) {
                    int index = k + nc_ * (j + nc_ * i);
                    std::array<double, 3> k_vec = getKVector(i, j, k);
                    double k_squared = k_vec[0] * k_vec[0] + k_vec[1] * k_vec[1] + k_vec[2] * k_vec[2];
                    // 根据k的平方值调整势能缓冲区的值
                    if (k_squared == 0) { // 避免除以0
                        potential_buffer_[index][0] = 0;
                        potential_buffer_[index][1] = 0;
                    } else {
                        double scaling_factor = -4 * M_PI / k_squared;
                        potential_buffer_[index][0] *= scaling_factor;
                        potential_buffer_[index][1] *= scaling_factor;
                    }
                }
            }
        }

        // 创建FFT的逆向计划
        fftw_plan inverse_plan = fftw_plan_dft_3d(nc_, nc_, nc_, potential_buffer_, density_buffer_, FFTW_BACKWARD, FFTW_MEASURE);
        // 执行逆向FFT变换
        fftw_execute(inverse_plan);
        // 释放逆向计划资源
        fftw_destroy_plan(inverse_plan);
        
        // 对整个缓冲区应用归一化因子
        double norm_factor = 1.0 / (nc_ * nc_ * nc_);
        for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
            density_buffer_[i][0] *= norm_factor;
            // 虚部不需要，因为我们只关心实部
            density_buffer_[i][1] = 0;
        }
    }

    // 获取给定网格索引对应的k向量，考虑Nyquist频率和频率折叠
    std::array<double, 3> getKVector(int i, int j, int k) {
        double kx = (i < nc_ / 2) ? i : (i - nc_);
        double ky = (j < nc_ / 2) ? j : (j - nc_);
        double kz = (k < nc_ / 2) ? k : (k - nc_);
        return {2 * M_PI * kx / box_width_, 2 * M_PI * ky / box_width_, 2 * M_PI * kz / box_width_};
    }

private:
    std::vector<Particle> particles_;  // 模拟中的粒子数组
    double time_max_;                  // 模拟的最大时间
    double delta_t_;                   // 模拟的时间步长
    double box_width_;                 // 模拟的空间盒子宽度
    double expansion_factor_;          // 空间膨胀因子
    int nc_;                           // 网格单元数
    double particle_mass_;             // 粒子质量
    fftw_complex* density_buffer_;     // 密度缓冲区
    fftw_complex* potential_buffer_;   // 势能缓冲区
};

