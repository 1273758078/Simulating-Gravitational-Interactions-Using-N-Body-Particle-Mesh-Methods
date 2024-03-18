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