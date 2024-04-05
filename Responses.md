1.2 Representing Particles
To address the tasks presented in the assignment, I have implemented a Particle class that encapsulates both position and velocity as three-dimensional std::array<double, 3> objects. This choice leverages the efficiency and determinism of std::array over other containers like std::vector, which is crucial for simulations involving a large number of particles. The position and velocity are private members, ensuring encapsulation and data hiding which are best practices in object-oriented programming.

The Particle constructor takes the initial position and velocity as arguments, directly initializing the class members. This follows the principle of constructors doing all necessary initializations, making objects of Particle class ready-to-use immediately after instantiation.

Member functions getPosition() and getVelocity() provide a way to access the particle's state without modifying it, adhering to the principle of const-correctness. This is important for scenarios where particle properties need to be accessed without risk of altering the simulation state inadvertently.

The setPosition() method allows individual manipulation of particle position, offering the flexibility to set custom distributions or for particles to be repositioned as needed by the simulation.

The update() method implements the basic physics of motion, using the supplied acceleration and time step (delta_t) to update the particle's velocity and position. It directly modifies the position_ and velocity_ members, which is a straightforward and efficient approach to implementing the particle's state change.

These implementations are found in Simulation.hpp and Simulation.cpp, which should be included in the project files where the Particle class is utilized. The performance expectation is met by ensuring that all operations are done in-place with minimal overhead, and ease of use is provided by clear, self-describing methods that reflect the particle's physical properties and behaviors.


1.10 Complexity
1. Density Calculation: The complexity of the density calculation is proportional to the total number of particles $n_p$​ because each particle needs to be placed into its respective grid cell once. Therefore, the time complexity for the density calculation is $O(n_p)$.

2. Potential Calculation: The calculation of the potential involves a Fast Fourier Transform (FFT) and its inverse, which typically have a complexity of $O(nlogn)$. In this case, $n$ is the number of points in the grid $N_c$. Thus, the time complexity for the potential calculation is $O(N_clogN_c)$.

3. Gravitational Field Calculation: The calculation of the gravitational field involves computing the spatial gradient of the potential at each point in the grid. Since this is done independently for each point in the grid, the time complexity is proportional to the number of cells, i.e., $O(N_c)$.

4. Particle Update Function: Updating the state of each particle involves calculating their new positions and velocities based on their accelerations. This is independent for each particle, hence the time complexity is proportional to the total number of particles, i.e., $O(n_p)$.

Given that $n_p$ must always be larger than $N_c$ to ensure a smooth density function, we can expect the particle-mesh method to scale better than the particle-particle method for large simulations. The PP method has a time complexity of $O(n_p^2)$ which grows much faster with the number of particles. In contrast, the primary bottleneck in the PM method is the FFT computation, which scales more modestly, thus it is more efficient for large-scale systems.



1.11 Application 1: Visualising A Developing Universe
./build/bin/NBody_Visualiser -nc 101 -np 10 -t 1.5 -dt 0.01 -F 1.0 -o Images/F1.0 -s 93170929
Images can be found in "/workspaces/COMP0210Assignment2/Images/F1.0"

./build/bin/NBody_Visualiser -nc 101 -np 10 -t 1.5 -dt 0.01 -F 1.02 -o Images/F1.02 -s 93170929
Images can be found in "/workspaces/COMP0210Assignment2/Images/F1.02"


2.1 Shared Memory Parallelism

1. Parallelisable Functions:
Calculate Density: The process of calculating density involves traversing all particles and updating the density buffer, which is a parallelizable process because each particle's processing is independent.
Calculate Potential: The calculation of potential energy involves performing a Fourier transform and applying a scaling factor to the transformed data, and this process can also be parallelized.
Calculate Gradient: Calculating the potential energy gradient (i.e. acceleration) is achieved by calculating the spatial derivative of potential energy at each grid point, which can also be parallelized.
UpdateParticles: Updating the position and velocity of all particles based on their acceleration is a parallelizable process.

the timing of each of the functions for different numbers of threads:

Benchmarking calculateDensity with 1 threads.
Time = 0.000173007 Info: 

Benchmarking calculatePotential with 1 threads.
Time = 0.167691 Info: 

Benchmarking calculateGradient with 1 threads.
Time = 0.00343542 Info: 

Benchmarking updateParticles with 1 threads.
Time = 1.4085e-05 Info: 

Benchmarking calculateDensity with 2 threads.
Time = 0.000207253 Info: 

Benchmarking calculatePotential with 2 threads.
Time = 0.00351921 Info: 

Benchmarking calculateGradient with 2 threads.
Time = 0.00171919 Info: 

Benchmarking updateParticles with 2 threads.
Time = 1.4056e-05 Info: 

Benchmarking calculateDensity with 4 threads.
Time = 0.000244718 Info: 

Benchmarking calculatePotential with 4 threads.
Time = 0.00387202 Info: 

Benchmarking calculateGradient with 4 threads.
Time = 0.000892769 Info: 

Benchmarking updateParticles with 4 threads.
Time = 1.1765e-05 Info: 

Benchmarking calculateDensity with 8 threads.
Time = 0.00030881 Info: 

Benchmarking calculatePotential with 8 threads.
Time = 0.0039106 Info: 

Benchmarking calculateGradient with 8 threads.
Time = 0.000850439 Info: 

Benchmarking updateParticles with 8 threads.
Time = 1.4669e-05 Info:

calculateDensity: Implemented with OpenMP, showing slight increases in execution time with more threads, which could indicate overhead from thread management or data access patterns.

updateParticles: Very short execution times with minimal changes across different thread counts. This suggests the operation is not a major bottleneck and thus doesn't benefit much from parallelisation.

For calculateDensity, the slight increase in execution time with more threads may necessitate investigating potential data contention, memory access patterns, or the overhead of managing a larger number of threads.

updateParticles being extremely fast already means that parallelisation does not significantly impact overall performance, indicating that optimization efforts might be better focused elsewhere.

2.2 Application 2: Comparing Universes with Distributed Memory
4. The correlationFunction was modified in several ways to address potential issues leading to incorrect results, specifically, infinite or negative infinite values. The modifications are as follows:

(1) Ensure CR[idx] is positive and non-zero before taking its logarithm: The code now checks if CR[idx] is greater than zero before applying the logarithm. If CR[idx] is zero or negative, which can happen due to calculation or rounding errors, it assigns a very small positive number (std::numeric_limits<double>::min()) to avoid taking the logarithm of zero or negative numbers, which is undefined or results in negative infinity.

(2) Handle edge conditions: The calculation of distance r and its corresponding index idx is adjusted to ensure that idx does not exceed the bounds of the CR array. Additionally, the condition checks if r is greater than zero to prevent division by zero errors when calculating CR[idx].

(3) Optimize the use of the logarithm function: When dealing with the logarithm's input value, if the value is zero, the code now opts to assign a minimal positive value before taking the logarithm. This approach avoids negative infinity results, which are not meaningful in the context of this calculation.

6. command line arguments:
mpirun -np 4 ./build/bin/NBody_Comparison Correlations 1.0 1.04