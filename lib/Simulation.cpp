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

// Implementation of the Particle class
Particle::Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity)
    : position_(position), velocity_(velocity) {}

// Get the position of the particle
std::array<double, 3> Particle::getPosition() const {
    return position_;
}

// Get the velocity of the particle
std::array<double, 3> Particle::getVelocity() const {
    return velocity_;
}

// Set the position of the particle
void Particle::setPosition(const std::array<double, 3>& position) {
    position_ = position;
}

// Set the velocity of the particle
void Particle::setVelocity(const std::array<double, 3>& velocity) {
    velocity_ = velocity;
}

// Update the position and velocity of the particle based on acceleration and time step
void Particle::update(double delta_t, const std::array<double, 3>& acceleration) {
    for (int i = 0; i < 3; ++i) {
        velocity_[i] += acceleration[i] * delta_t;
        position_[i] += velocity_[i] * delta_t;
        
        // Ensure the position wraps around the simulation box bounds
        if (position_[i] < 0) position_[i] += 1.0;  // If position is negative, wrap around to the positive side
        else if (position_[i] >= 1.0) position_[i] -= 1.0;  // If position exceeds 1.0, wrap back to the start
    }
}

// Implementation of the Simulation class
Simulation::Simulation(double time_max, double delta_t, double box_width, double expansion_factor, int nc, double particle_mass)
    : time_max_(time_max), delta_t_(delta_t), box_width_(box_width), expansion_factor_(expansion_factor), nc_(nc), particle_mass_(particle_mass) {
    // Allocate memory for the density and potential energy buffers
    density_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
    k_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
    potential_buffer_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nc_ * nc_ * nc_);
    // Initialize the density buffer
    initializeDensityBuffer();
}

// Destructor, free allocated memory
Simulation::~Simulation() {
    if (density_buffer_ != nullptr) {
        fftw_free(density_buffer_);
    }
    if (k_buffer_ != nullptr) {
        fftw_free(k_buffer_);
    }
    if (potential_buffer_ != nullptr) {
        fftw_free(potential_buffer_);
    }
}

// Add a particle to the simulation
void Simulation::addParticle(const std::array<double, 3>& position) {
    particles_.emplace_back(position, std::array<double, 3>{0.0, 0.0, 0.0});
}

// Initialize particles with a given number and seed, randomly distributing them
void Simulation::initializeParticles(int num_particles, unsigned seed) {
    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0); // Ensure particle positions are within the box width
    for (int i = 0; i < num_particles; ++i) {
        addParticle({dis(gen), dis(gen), dis(gen)});
    }
}

// Run the simulation
void Simulation::run(const std::optional<std::string>& output_folder) {
    int steps = 0;
    int save_interval = 10;
    for (double t = 0; t < time_max_; t += delta_t_) {
        calculateDensity();
        calculatePotential();
        auto gradients = calculateGradient(potential_buffer_);
        updateParticles(gradients, delta_t_);
        expandBox(expansion_factor_);

        // If an output folder is provided and the current step is a multiple of the save interval, save the density image
        if (output_folder && steps % save_interval == 0) {
            std::string filename = *output_folder + "/density_" + std::to_string(steps) + ".pbm";
            SaveToFile(density_buffer_, nc_, filename);
        }
        steps++;
    }
}

// Initialize the density buffer
void Simulation::initializeDensityBuffer() {
    for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
        density_buffer_[i][0] = 0.0;
        density_buffer_[i][1] = 0.0;
    }
}

// Calculate the density based on the positions of the particles
void Simulation::calculateDensity() {
    initializeDensityBuffer();
    double cell_volume = std::pow(box_width_ / nc_, 3);
    for (const auto& particle : particles_) {
        std::array<double, 3> position = particle.getPosition();

        // Calculate the grid index for the particle's position, ensuring the index does not exceed the boundaries
        int i = std::min(static_cast<int>(position[0] * nc_), nc_ - 1);
        int j = std::min(static_cast<int>(position[1] * nc_), nc_ - 1);
        int k = std::min(static_cast<int>(position[2] * nc_), nc_ - 1);

        // Ensure the index is not negative
        i = std::max(i, 0);
        j = std::max(j, 0);
        k = std::max(k, 0);

        int index = (k + nc_ * (j + nc_ * i));
        density_buffer_[index][0] += particle_mass_ / cell_volume;
    }
}

// Calculate the density based on the positions of the particles(OpenMP version)
// void Simulation::calculateDensity() {
//     initializeDensityBuffer();
//     double cell_volume = std::pow(box_width_ / nc_, 3);

//     #pragma omp parallel for
//     for (const auto& particle : particles_) {
//         std::array<double, 3> position = particle.getPosition();

//         // Calculate the grid index for the particle's position, ensuring the index does not exceed the boundaries
//         int i = std::min(static_cast<int>(position[0] * nc_), nc_ - 1);
//         int j = std::min(static_cast<int>(position[1] * nc_), nc_ - 1);
//         int k = std::min(static_cast<int>(position[2] * nc_), nc_ - 1);

//         // Ensure the index is not negative
//         i = std::max(i, 0);
//         j = std::max(j, 0);
//         k = std::max(k, 0);

//         int index = (k + nc_ * (j + nc_ * i)) % (nc_ * nc_ * nc_);
//         density_buffer_[index][0] += particle_mass_ / cell_volume;
//     }
// }

// Get the index of a particular cell
int Simulation::getCellIndex(double x, double y, double z){

    // Calculate the grid index for the particle's position, ensuring the index does not exceed the boundaries
    int i = std::min(static_cast<int>(x * nc_), nc_ - 1);
    int j = std::min(static_cast<int>(y * nc_), nc_ - 1);
    int k = std::min(static_cast<int>(z * nc_), nc_ - 1);

    // Ensure the index is not negative
    i = std::max(i, 0);
    j = std::max(j, 0);
    k = std::max(k, 0);

    int index = (k + nc_ * (j + nc_ * i));
    return index;
}

// Calculate potential energy
// void Simulation::calculatePotential() {
//     // Create and execute the forward FFT plan, transforming density data to frequency space
//     fftw_plan forward_plan = fftw_plan_dft_3d(nc_, nc_, nc_, density_buffer_, k_buffer_, FFTW_FORWARD, FFTW_MEASURE);
//     fftw_execute(forward_plan);
//     // Immediately release FFT plan resources after execution
//     fftw_destroy_plan(forward_plan);

//     // Process potential energy data in frequency space
//     for (int i = 0; i < nc_; ++i) {
//         for (int j = 0; j < nc_; ++j) {
//             for (int k = 0; k < nc_; ++k) {
//                 int index = k + nc_ * (j + nc_ * i);
//                 // Calculate the k-vector for each point
//                 // std::array<double, 3> k_vec = getKVector(i, j, k);没用
//                 // Calculate the square of the k-vector
//                 // double k_squared = k_vec[0] * k_vec[0] + k_vec[1] * k_vec[1] + k_vec[2] * k_vec[2];
//                 double k_squared = (i * i + j * j + k * k) / (box_width_ * box_width_);
//                 // Avoid dividing by zero
//                 if (k_squared == 0) {
//                     k_buffer_[index][0] = 0;
//                     k_buffer_[index][1] = 0;
//                 } else {
//                     // Adjusting the potential energy value based on the square of the k-vector
//                     double scaling_factor = -4 * M_PI / k_squared;
//                     k_buffer_[index][0] *= scaling_factor;
//                     k_buffer_[index][1] *= scaling_factor;
//                 }
//             }
//         }
//     }

//     // Create and execute the inverse FFT plan, transforming potential energy data back to real space
//     fftw_plan inverse_plan = fftw_plan_dft_3d(nc_, nc_, nc_, k_buffer_, potential_buffer_, FFTW_BACKWARD, FFTW_MEASURE);
//     fftw_execute(inverse_plan);
//     // Immediately release FFT plan resources after execution
//     fftw_destroy_plan(inverse_plan);

//     // Apply a normalization factor to ensure correct results after the inverse FFT
//     double norm_factor = 1.0 / 8 * (nc_ * nc_ * nc_);
//     for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
//         potential_buffer_[i][0] *= norm_factor;
//         // For potential energy, we only care about the real part
//         potential_buffer_[i][1] = 0;
//     }
// }

// Calculate potential energy(OpenMP version)
void Simulation::calculatePotential() {
    // Create and execute the forward FFT plan, transforming density data to frequency space
    fftw_plan forward_plan = fftw_plan_dft_3d(nc_, nc_, nc_, density_buffer_, k_buffer_, FFTW_FORWARD, FFTW_MEASURE);
    fftw_execute(forward_plan);
    // Immediately release FFT plan resources after execution
    fftw_destroy_plan(forward_plan);

    // Process potential energy data in frequency space
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

    // Create and execute the inverse FFT plan, transforming potential energy data back to real space
    fftw_plan inverse_plan = fftw_plan_dft_3d(nc_, nc_, nc_, k_buffer_, potential_buffer_, FFTW_BACKWARD, FFTW_MEASURE);
    fftw_execute(inverse_plan);
    // Immediately release FFT plan resources after execution
    fftw_destroy_plan(inverse_plan);

    // Apply a normalization factor to ensure correct results after the inverse FFT
    double norm_factor = 1.0 / (8 * (nc_ * nc_ * nc_));
    for (int i = 0; i < nc_ * nc_ * nc_; ++i) {
        potential_buffer_[i][0] *= norm_factor;
        // For potential energy, we only care about the real part
        potential_buffer_[i][1] = 0;
    }
}

// Get the potential energy value at a given grid index
double Simulation::getPotentialAtGridIndex(int i, int j, int k) {
    int index = k + nc_ * (j + nc_ * i);
    // Note: We use potential_buffer_ to store the results of the inverse FFT
    return potential_buffer_[index][0]; // We are only interested in the real part
}

// Calculate the gradient
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

// Calculate the gradient(OpenMP version)
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

// Wrap two boundaries together
int Simulation::wrapIndex(int index, int max) {
    return (index + max) % max; // Modular operation
}

// Update the state of the particle swarm
void Simulation::updateParticles(const std::vector<std::array<double, 3>>& gradients, double delta_t) {
    for (Particle& particle : particles_) {
        std::array<double, 3> pos = particle.getPosition();
        
        // Calculate the grid index for the particle's position, ensuring the index does not exceed the boundaries
        int i = std::min(static_cast<int>(pos[0] * nc_), nc_ - 1);
        int j = std::min(static_cast<int>(pos[1] * nc_), nc_ - 1);
        int k = std::min(static_cast<int>(pos[2] * nc_), nc_ - 1);

        // Ensure the index is not negative
        i = std::max(i, 0);
        j = std::max(j, 0);
        k = std::max(k, 0);

        // Get the acceleration at the grid point
        std::array<double, 3> acceleration = gradients[i * nc_ * nc_ + j * nc_ + k];

        // Update the particle state
        particle.update(delta_t, acceleration);

        // Re-obtain the updated position
        pos = particle.getPosition();

        // Apply periodic boundary conditions
        for (int dim = 0; dim < 3; ++dim) {
            if (pos[dim] < 0) pos[dim] += 1.0;  // If the position is less than 0, wrap to the positive side
            else if (pos[dim] >= 1.0) pos[dim] -= 1.0;  // If the position is greater than or equal to 1, wrap back to the starting side
        }

        // Update the particle position with the updated, boundary condition-considered position
        particle.setPosition(pos);
    }
}

// Update the state of the particle swarm(OpenMP version)
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
    // Enlarge the box width
    box_width_ *= expansion_factor;
    
    // If there's a stored cell width, it also needs to be updated
    cell_width_ = box_width_ / nc_;
    
    // Reduce the velocity of all particles
    for (auto& particle : particles_) {
        std::array<double, 3> velocity = particle.getVelocity();
        velocity[0] /= expansion_factor;
        velocity[1] /= expansion_factor;
        velocity[2] /= expansion_factor;
        particle.setVelocity(velocity);
    }
}

// Get the total number of cells
int Simulation::getTotalCells() const {
    return nc_ * nc_ * nc_;
}

// Get the data of the density function
fftw_complex* Simulation::getDensityBuffer() const {
    return density_buffer_;
}

// Get the mass of a particle
double Simulation::getParticleMass() const {
    return particle_mass_;
}

// Get the volume of a cell
double Simulation::getCellVolume() const {
    // Assuming you've already calculated the cell volume
    return std::pow(box_width_ / nc_, 3);
}

// Extract the positions of all particles, returning a vector
std::vector<std::array<double, 3>> Simulation::getParticlesPositions() const {
    std::vector<std::array<double, 3>> positions;
    positions.reserve(particles_.size()); // Reserve sufficient space to avoid multiple memory allocations

    for (const auto& particle : particles_) {
        positions.push_back(particle.getPosition());    
    }
    return positions;
}
