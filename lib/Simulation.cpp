// Simulation.cpp
#include "Simulation.hpp"
#include "Utils.hpp"
#include <fftw3.h>
#include <vector>
#include <array>
#include <random>
#include <cmath>
#include <algorithm>
#include <optional>
#include <string>
#include <iostream>

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

// 设置粒子的速度
void Particle::setVelocity(const std::array<double, 3>& velocity) {
    velocity_ = velocity;
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
    k_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
    potential_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
    // 初始化密度缓冲区
    initializeDensityBuffer();
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

    // // 打印最多前五个粒子的位置
    // for (size_t i = 0; i < particles_.size() && i < 5; i++) {
    //     std::array<double, 3> pos = particles_[i].getPosition(); // 使用getPosition()获取位置
    //     std::cout << "Particle " << i << ": (" 
    //               << pos[0] << ", " // 使用正确的索引访问数组元素
    //               << pos[1] << ", " 
    //               << pos[2] << ")" << std::endl;
    // }


}

// 运行模拟
void Simulation::run(const std::optional<std::string>& output_folder) {
    int steps = 0;
    int save_interval = 10;
    for (double t = 0; t < time_max_; t += delta_t_) {
        calculateDensity();
        calculatePotential();
        auto gradients = calculateGradient(potential_buffer_);
        
        
        std::vector<std::array<double, 3>> positions = getParticlesPositions();
        // // 打印前五个粒子的位置（如果它们存在）
        // for (size_t i = 0; i < positions.size() && i < 5; ++i) {
        //     std::cout << "Particle " << i << ": (" 
        //             << positions[i][0] << ", " 
        //             << positions[i][1] << ", " 
        //             << positions[i][2] << ")" << std::endl;
        // }





        updateParticles(gradients, delta_t_);
        expandBox(expansion_factor_);

        // 如果提供了输出文件夹并且当前步骤是保存间隔的倍数，则保存密度图像
        if (output_folder && steps % save_interval == 0) {
            std::string filename = *output_folder + "/density_" + std::to_string(steps) + ".ppm";
            SaveToFile(density_buffer_, nc_, filename);
        }
        steps++;
    }
}


// 初始化密度缓冲区
void Simulation::initializeDensityBuffer() {
    for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
        density_buffer_[i][0] = 0.0;
        density_buffer_[i][1] = 0.0;
    }
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

// OpenMP版
// void Simulation::calculateDensity() {
//     initializeDensityBuffer();
//     double cell_volume = std::pow(box_width_ / nc_, 3);

//     #pragma omp parallel for
//     for (const auto& particle : particles_) {
//         std::array<double, 3> position = particle.getPosition();
//         int i = static_cast<int>(position[0] * nc_ / box_width_);
//         int j = static_cast<int>(position[1] * nc_ / box_width_);
//         int k = static_cast<int>(position[2] * nc_ / box_width_);
//         int index = (k + nc_ * (j + nc_ * i)) % (nc_ * nc_ * nc_);
//         density_buffer_[index][0] += particle_mass_ / cell_volume;
//     }
// }


//得到某一个单元格的索引
int Simulation::getCellIndex(double x, double y, double z){
    int i = static_cast<int>(x * nc_ / box_width_);
    int j = static_cast<int>(y * nc_ / box_width_);
    int k = static_cast<int>(z * nc_ / box_width_);
    int index = (k + nc_ * (j + nc_ * i)) % (nc_ * nc_ * nc_);
    return index;
}

// 计算势能
// void Simulation::calculatePotential() {
//     // 创建并执行正向FFT计划，将密度数据转换到频率空间
//     fftw_plan forward_plan = fftw_plan_dft_3d(nc_, nc_, nc_, density_buffer_, k_buffer_, FFTW_FORWARD, FFTW_MEASURE);
//     fftw_execute(forward_plan);
//     // 完成后立即释放FFT计划资源
//     fftw_destroy_plan(forward_plan);

//     // 在频率空间处理势能数据
//     for (int i = 0; i < nc_; ++i) {
//         for (int j = 0; j < nc_; ++j) {
//             for (int k = 0; k < nc_; ++k) {
//                 int index = k + nc_ * (j + nc_ * i);
//                 // 计算每个点的k向量
//                 // std::array<double, 3> k_vec = getKVector(i, j, k);没用
//                 // 计算k向量的平方
//                 // double k_squared = k_vec[0] * k_vec[0] + k_vec[1] * k_vec[1] + k_vec[2] * k_vec[2];
//                 double k_squared = (i * i + j * j + k * k) / (box_width_ * box_width_);
//                 // 避免除以零
//                 if (k_squared == 0) {
//                     k_buffer_[index][0] = 0;
//                     k_buffer_[index][1] = 0;
//                 } else {
//                     // 根据k向量的平方调整势能值
//                     double scaling_factor = -4 * M_PI / k_squared;
//                     k_buffer_[index][0] *= scaling_factor;
//                     k_buffer_[index][1] *= scaling_factor;
//                 }
//             }
//         }
//     }

//     // 创建并执行逆向FFT计划，将势能数据从频率空间转换回实空间
//     fftw_plan inverse_plan = fftw_plan_dft_3d(nc_, nc_, nc_, k_buffer_, potential_buffer_, FFTW_BACKWARD, FFTW_MEASURE);
//     fftw_execute(inverse_plan);
//     // 完成后立即释放FFT计划资源
//     fftw_destroy_plan(inverse_plan);

//     // 应用归一化因子，以确保逆FFT后的结果正确
//     double norm_factor = 1.0 / 8 * (nc_ * nc_ * nc_);
//     for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
//         potential_buffer_[i][0] *= norm_factor;
//         // 对于势能，我们只关心实部
//         potential_buffer_[i][1] = 0;
//     }
// }

// OpenMP版
void Simulation::calculatePotential() {
    // 创建并执行正向FFT计划，将密度数据转换到频率空间
    fftw_plan forward_plan = fftw_plan_dft_3d(nc_, nc_, nc_, density_buffer_, k_buffer_, FFTW_FORWARD, FFTW_MEASURE);
    fftw_execute(forward_plan);
    // 完成后立即释放FFT计划资源
    fftw_destroy_plan(forward_plan);

    // 在频率空间处理势能数据
    #pragma omp parallel for collapse(3)
    for (int i = 0; i < nc_; ++i) {
        for (int j = 0; j < nc_; ++j) {
            for (int k = 0; k < nc_; ++k) {
                int index = k + nc_ * (j + nc_ * i);
                double k_squared = (i * i + j * j + k * k) / (box_width_ * box_width_);
                if (k_squared == 0) {
                    k_buffer_[index][0] = 0;
                    k_buffer_[index][1] = 0;
                } else {
                    double scaling_factor = -4 * M_PI / k_squared;
                    k_buffer_[index][0] *= scaling_factor;
                    k_buffer_[index][1] *= scaling_factor;
                }
            }
        }
    }

    // 创建并执行逆向FFT计划，将势能数据从频率空间转换回实空间
    fftw_plan inverse_plan = fftw_plan_dft_3d(nc_, nc_, nc_, k_buffer_, potential_buffer_, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_execute(inverse_plan);
    // 完成后立即释放FFT计划资源
    fftw_destroy_plan(inverse_plan);

    // 应用归一化因子，以确保逆FFT后的结果正确
    double norm_factor = 1.0 / 8 * (nc_ * nc_ * nc_);
    for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
        potential_buffer_[i][0] *= norm_factor;
        // 对于势能，我们只关心实部
        potential_buffer_[i][1] = 0;
    }
}


// 在给定网格索引处获取势能值
double Simulation::getPotentialAtGridIndex(int i, int j, int k) {
    int index = k + nc_ * (j + nc_ * i);
    // 注意：这里使用potential_buffer_来保存逆FFT后的结果
    return potential_buffer_[index][0]; // 我们只关心实部
}

// 计算梯度
// std::vector<std::array<double, 3>> Simulation::calculateGradient(const fftw_complex* potential) {
//     std::vector<std::array<double, 3>> gradient(nc_ * nc_ * nc_);
    
//     for (int i = 0; i < nc_; ++i) {
//         for (int j = 0; j < nc_; ++j) {
//             for (int k = 0; k < nc_; ++k) {
//                 int index = (k + nc_ * (j + nc_ * i));
//                 gradient[index][0] = -(potential[wrapIndex(i+1, nc_) * nc_ * nc_ + j * nc_ + k][0] - potential[wrapIndex(i-1, nc_) * nc_ * nc_ + j * nc_ + k][0]) / (2 * box_width_ / nc_);
//                 gradient[index][1] = -(potential[i * nc_ * nc_ + wrapIndex(j+1, nc_) * nc_ + k][0] - potential[i * nc_ * nc_ + wrapIndex(j-1, nc_) * nc_ + k][0]) / (2 * box_width_ / nc_);
//                 gradient[index][2] = -(potential[i * nc_ * nc_ + j * nc_ + wrapIndex(k+1, nc_)][0] - potential[i * nc_ * nc_ + j * nc_ + wrapIndex(k-1, nc_)][0]) / (2 * box_width_ / nc_);
//             }
//         }
//     }
//     return gradient;
// }

// OpenMP版
std::vector<std::array<double, 3>> Simulation::calculateGradient(const fftw_complex* potential) {
    std::vector<std::array<double, 3>> gradient(nc_ * nc_ * nc_);
    
    #pragma omp parallel for collapse(3)
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


// 将两个边界接在一起
int Simulation::wrapIndex(int index, int max) {
    return (index + max) % max; // 模运算
}

// 更新粒子群状态
void Simulation::updateParticles(const std::vector<std::array<double, 3>>& gradients, double delta_t) {
    for (Particle& particle : particles_) {
        std::array<double, 3> pos = particle.getPosition();
        
        // 计算粒子位置所在的格点索引，同时确保索引不会超出边界
        int i = std::min(static_cast<int>(pos[0] * nc_ / box_width_), nc_ - 1);
        int j = std::min(static_cast<int>(pos[1] * nc_ / box_width_), nc_ - 1);
        int k = std::min(static_cast<int>(pos[2] * nc_ / box_width_), nc_ - 1);

        // 确保索引不会为负
        i = std::max(i, 0);
        j = std::max(j, 0);
        k = std::max(k, 0);

        // 获得格点上的加速度
        std::array<double, 3> acceleration = gradients[i * nc_ * nc_ + j * nc_ + k];

        // 更新粒子状态
        particle.update(delta_t, acceleration);
    }
}

// OpenMP版
// void Simulation::updateParticles(const std::vector<std::array<double, 3>>& gradients, double delta_t) {
//     #pragma omp parallel for
//     for (size_t idx = 0; idx < particles_.size(); ++idx) {
//         Particle& particle = particles_[idx];
//         std::array<double, 3> pos = particle.getPosition();
        
//         int i = static_cast<int>(pos[0] * nc_ / box_width_);
//         int j = static_cast<int>(pos[1] * nc_ / box_width_);
//         int k = static_cast<int>(pos[2] * nc_ / box_width_);

//         int index = i * nc_ * nc_ + j * nc_ + k;

//         std::array<double, 3> acceleration = gradients[index];

//         particle.update(delta_t, acceleration);
//     }
// }




void Simulation::expandBox(double expansion_factor) {
    // 放大盒子的宽度
    box_width_ *= expansion_factor;
    
    // 如果有存储单元格宽度，也需要更新它
    cell_width_ = box_width_ / nc_;
    
    // 缩小所有粒子的速度
    for (auto& particle : particles_) {
        std::array<double, 3> velocity = particle.getVelocity();
        velocity[0] /= expansion_factor;
        velocity[1] /= expansion_factor;
        velocity[2] /= expansion_factor;
        particle.setVelocity(velocity);
    }
}

// 得到单元格的总数
int Simulation::getTotalCells() const {
    return nc_ * nc_ * nc_;
}

// 得到密度函数的数据
fftw_complex* Simulation::getDensityBuffer() const {
    return density_buffer_;
}

// 得到粒子的质量
double Simulation::getParticleMass() const {
    return particle_mass_;
}

// 得到单元格的体积
double Simulation::getCellVolume() const {
    // 假设你已经计算了单元格体积
    return std::pow(box_width_ / nc_, 3);
}

// 提取所有粒子的位置，返回一个向量
std::vector<std::array<double, 3>> Simulation::getParticlesPositions() const {
    std::vector<std::array<double, 3>> positions;
    positions.reserve(particles_.size()); // 预留足够的空间以避免多次内存分配

    for (const auto& particle : particles_) {
        positions.push_back(particle.getPosition());
        std::cout << "Particle " << ": (" 
                    << particle.getPosition()[0] << ", " // 使用正确的索引访问数组元素
                    << particle.getPosition()[1] << ", " 
                    << particle.getPosition()[2] << ")" << std::endl;    
    }

    return positions;
}