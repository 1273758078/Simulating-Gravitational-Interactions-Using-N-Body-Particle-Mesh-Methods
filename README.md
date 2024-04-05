# UniverseInABox
Assignment to write a simple Particle-Mesh gravitational simulation

# Building and Running the Simulation
This project simulates gravitational interactions between particles in an expanding universe using both serial and parallel computing techniques. It leverages the FFTW library for fast Fourier transforms and OpenMP/MPI for parallelism.

# Requirements:
- C++17 compiler (g++ recommended)
- FFTW3 library
- MPI (for parallel distributed memory execution)
- OpenMP (for shared memory parallelism)
- Catch2 (for testing, optional)

# Setup:

1. Install FFTW3:
sudo apt-get install libfftw3-dev

2. Clone the repository:
git clone git@github.com:UCL-COMP0210-23-24/universe-in-a-box-1273758078.git COMP0210Assignment2

3. Navigate to the project directory:
cd COMP0210Assignment2

4. Build the project using CMake:
cmake -B build
cmake --build build

# Running the Simulation:
- To visualize the universe's development:
./build/bin/NBody_Visualiser -nc <cells> -np <particles per cell> -t <total time> -dt <time step> -F <expansion factor> -o <output folder> -s <random seed>

Example:
./build/bin/NBody_Visualiser -nc 101 -np 10 -t 1.5 -dt 0.01 -F 1.02 -o Images -s 42

- To compare universes with different expansion factors:
mpirun -np <number of processes> ./build/bin/NBody_Comparison <output folder> <min expansion factor> <max expansion factor> <>

Example:
mpirun -np 4 ./build/bin/NBody_Comparison Correlations 1.0 1.04

# Visualizing Output:
The simulation outputs are saved in the specified output folder. You can view the evolving universe through the generated images.