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

// Particle class
class Particle {
public:
    Particle(const std::array<double, 3>& position, const std::array<double, 3>& velocity);
    std::array<double, 3> getPosition() const;
    std::array<double, 3> getVelocity() const;
    void setPosition(const std::array<double, 3>& position);
    void setVelocity(const std::array<double, 3>& velocity);
    void update(double delta_t, const std::array<double, 3>& acceleration);

private:
    std::array<double, 3> position_; // Particle position
    std::array<double, 3> velocity_; // Particle velocity
};

// Simulation class
class Simulation {
public:
    Simulation(double time_max, double delta_t, double box_width, double expansion_factor, int nc, double particle_mass);
    ~Simulation();
    void addParticle(const std::array<double, 3>& position);
    void initializeParticles(int num_particles, unsigned seed);
    void run(const std::optional<std::string>& output_folder);
    void initializeDensityBuffer();
    void calculateDensity();
    // Get the index of a cell
    int getCellIndex(double x, double y, double z);
    void calculatePotential();
    // Get the potential energy at a given grid index (i, j, k)
    double getPotentialAtGridIndex(int i, int j, int k);
    // Function for gradient calculation
    std::vector<std::array<double, 3>> calculateGradient(const fftw_complex* potential);
    // Function to update all particles
    void updateParticles(const std::vector<std::array<double, 3>>& gradients, double delta_t);
    // Function to expand the simulation box
    void expandBox(double expansion_factor);
    // Get the total number of cells
    int getTotalCells() const;
    // Get the data of the density function
    fftw_complex* getDensityBuffer() const;
    // Wrap two boundaries together
    int wrapIndex(int index, int max);
    // Get the mass of a particle
    double getParticleMass() const;
    // Get the volume of a cell
    double getCellVolume() const;
    // Extract the positions of all particles, returning a vector
    std::vector<std::array<double, 3>> getParticlesPositions() const;

private:
    std::vector<Particle> particles_;  // Array of particles in the simulation
    double time_max_;                  // Maximum simulation time
    double delta_t_;                   // Simulation time step
    double box_width_;                 // Width of the simulation space box
    double expansion_factor_;          // Space expansion factor
    int nc_;                           // Number of grid cells
    double particle_mass_;             // Particle mass
    fftw_complex* density_buffer_;     // Density buffer
    fftw_complex* k_buffer_; // Transition area
    fftw_complex* potential_buffer_;   // Potential energy buffer
    double cell_width_; // Width of a cell, if applicable
};

#endif // SIMULATION_HPP
