// test_simulation.cpp
#define CATCH_CONFIG_MAIN  // This macro directive tells Catch to generate a main function
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <cmath>
#include "Simulation.hpp"  // Make sure this is the correct path to your Simulation class

using namespace Catch; //::Matchers

TEST_CASE("Test potential function for single particle", "[Potential Tests]") {
    // A single particle is located at the center of the box with a mass of 0.01
    double mass = 0.01;
    double time_max = 1.0; // Used only to initialize the Simulation object
    double delta_t = 0.1;  // Same as above, time step
    double width = 100;
    double expansion_factor = 1.0; // The expansion factor is not used in this test, but is required for the Simulation constructor
    int ncells = 101;
    double particle_mass = 0.01; // Particle mass

    // Declare a Simulation object with initial parameters
    Simulation sim(time_max, delta_t, width, expansion_factor, ncells, particle_mass);
    sim.addParticle({0.5, 0.5, 0.5}); // Add a single particle at the center of the box
    sim.calculateDensity(); // Calculate density
    sim.calculatePotential(); // Calculate potential energy

    double w_c = width / ncells;
    // Examine the potential function along the x-axis
    for (int i = 0; i < 101; i++) {
        int j = 50;
        int k = 50;
        if (i == 50 && j == 50 && k == 50) {
            // Ignore the center point, as the point potential here would be singular
            continue;
        } else {
            double pot = sim.getPotentialAtGridIndex(i, j, k); // Get potential energy
            // Considering periodicity, approximate the value of multiple source potential
            double dx1 = std::abs(i - 50) * w_c;
            double dx2 = std::abs(i - 151) * w_c;
            double dx3 = std::abs(i + 51) * w_c;
            double expected_pot =  -mass * (1/dx1 + 1/dx2 + 1/dx3);
            REQUIRE_THAT(pot, Matchers::WithinRel(expected_pot, 0.3)); // Verify that the potential is within the expected relative error range
        }
    }
}

TEST_CASE("Imaginary components are zero", "[density]") {
    // Set initial parameters for the simulation
    const double t_max = 10.0;  // The maximum time of the simulation
    const double delta_t = 0.1;  // Time step
    const double box_width = 100.0; // Box width
    const double expansion_factor = 1.0; // Expansion factor
    const int nc = 10; // Number of grid cells per side
    const double particle_mass = 1.0; // Particle mass

    // Create a simulation instance without particles
    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.initializeDensityBuffer(); // Initialize the density buffer
    sim.calculateDensity(); // Calculate density

    for (int i = 0; i < sim.getTotalCells(); ++i) {
        REQUIRE(sim.getDensityBuffer()[i][1] == 0.0);
    }
}

TEST_CASE("Density is zero without particles", "[density]") {
    // Set initial parameters for the simulation
    const double t_max = 10.0;  // The maximum time of the simulation
    const double delta_t = 0.1;  // Time step
    const double box_width = 100.0; // Box width
    const double expansion_factor = 1.0; // Expansion factor
    const int nc = 10; // Number of grid cells per side
    const double particle_mass = 1.0; // Particle mass

    // Create a simulation instance without particles
    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.initializeDensityBuffer(); // Initialize the density buffer
    sim.calculateDensity(); // Calculate density

    for (int i = 0; i < sim.getTotalCells(); ++i) {
        REQUIRE(sim.getDensityBuffer()[i][0] == 0.0);
    }
}

TEST_CASE("Density with a single particle", "[density]") {
    // Set initial parameters for the simulation
    const double t_max = 10.0;  // The maximum time of the simulation
    const double delta_t = 0.1;  // Time step
    const double box_width = 100.0; // Box width
    const double expansion_factor = 1.0; // Expansion factor
    const int nc = 10; // Number of grid cells per side
    const double particle_mass = 1.0; // Particle mass

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
    // Set initial parameters for the simulation
    const double t_max = 10.0;  // The maximum time of the simulation
    const double delta_t = 0.1;  // Time step
    const double box_width = 100.0; // Box width
    const double expansion_factor = 1.0; // Expansion factor
    const int nc = 10; // Number of grid cells per side
    const double particle_mass = 1.0; // Particle mass

    Simulation sim(t_max, delta_t, box_width, expansion_factor, nc, particle_mass);
    sim.addParticle({0.1, 0.1, 0.1});
    sim.addParticle({0.1, 0.1, 0.1}); // Same cell as the first particle
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

// ...include other test cases...

TEST_CASE("Gradient of potential is calculated correctly", "[Simulation]") {
    // Initialize a Simulation object and potential buffer for the test
    int nc = 10;
    double box_width = 100.0;
    Simulation sim(1.0, 0.1, box_width, 1.0, nc, 1.0);
    fftw_complex* potential = new fftw_complex[nc * nc * nc]();
    
    // Set up a simple test case: non-zero potential at the center, zero elsewhere
    int centerIndex = (nc / 2) * nc * nc + (nc / 2) * nc + (nc / 2);
    potential[centerIndex][0] = 100.0;  // Set the real part only
    
    // Calculate the gradient
    auto gradient = sim.calculateGradient(potential);
    
    // Test: the gradient values around the center should be known
    REQUIRE(gradient[centerIndex][0] == Approx(0.0).margin(1e-5));
    REQUIRE(gradient[centerIndex][1] == Approx(0.0).margin(1e-5));
    REQUIRE(gradient[centerIndex][2] == Approx(0.0).margin(1e-5));
    
    // Cleanup
    delete[] potential;
}

TEST_CASE("Gradient calculation with periodic boundaries", "[Simulation]") {
    // Test boundary conditions
    int nc = 10;
    double box_width = 100.0;
    Simulation sim(1.0, 0.1, box_width, 1.0, nc, 1.0);
    fftw_complex* potential = new fftw_complex[nc * nc * nc]();
    
    // Set up a potential that linearly increases
    for (int i = 0; i < nc; ++i) {
        for (int j = 0; j < nc; ++j) {
            for (int k = 0; k < nc; ++k) {
                int index = (k + nc * (j + nc * i));
                potential[index][0] = static_cast<double>(i);
            }
        }
    }
    
    // Calculate the gradient and check if the edges wrap correctly
    auto gradient = sim.calculateGradient(potential);
    
    int edgeIndex = nc - 1;
    int wrappedIndex = sim.wrapIndex(edgeIndex + 1, nc); // Should be 0
    
    REQUIRE(wrappedIndex == 0);
    REQUIRE(gradient[edgeIndex][0] == Approx(-1.0).margin(1e-5));  // The gradient at the edge should be -1
    
    // Cleanup
    delete[] potential;
}

TEST_CASE("Particles are updated correctly", "[Simulation]") {
    // Set initial conditions
    Particle particle({0.5, 0.5, 0.5}, {0.0, 0.0, 0.0});
    std::array<double, 3> acceleration = {1.0, 0.0, 0.0}; // Assume a unit acceleration along the x-axis
    double delta_t = 1.0; // Time step of 1 second

    // Perform the update
    particle.update(delta_t, acceleration);

    // Check the results
    auto position = particle.getPosition();
    auto velocity = particle.getVelocity();

    REQUIRE(position[0] == Approx(0.5 + 1.0 * delta_t)); // Check if x position is correctly updated
    REQUIRE(velocity[0] == Approx(1.0 * delta_t)); // Check if x velocity is correctly updated
    // Since there was no acceleration along the y and z directions, their positions and velocities should not change
    REQUIRE(position[1] == Approx(0.5));
    REQUIRE(position[2] == Approx(0.5));
    REQUIRE(velocity[1] == Approx(0.0));
    REQUIRE(velocity[2] == Approx(0.0));
}
