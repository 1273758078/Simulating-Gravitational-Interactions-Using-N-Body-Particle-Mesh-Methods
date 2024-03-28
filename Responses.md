1.2 Representing Particles



1.10 Complexity
1. Density Calculation: The complexity of the density calculation is proportional to the total number of particles $n_p$​ because each particle needs to be placed into its respective grid cell once. Therefore, the time complexity for the density calculation is $O(n_p)$.

2. Potential Calculation: The calculation of the potential involves a Fast Fourier Transform (FFT) and its inverse, which typically have a complexity of $O(nlogn)$. In this case, $n$ is the number of points in the grid $N_c$. Thus, the time complexity for the potential calculation is $O(N_clogN_c)$.

3. Gravitational Field Calculation: The calculation of the gravitational field involves computing the spatial gradient of the potential at each point in the grid. Since this is done independently for each point in the grid, the time complexity is proportional to the number of cells, i.e., $O(N_c)$.

4. Particle Update Function: Updating the state of each particle involves calculating their new positions and velocities based on their accelerations. This is independent for each particle, hence the time complexity is proportional to the total number of particles, i.e., $O(n_p)$.

Given that $n_p$ must always be larger than $N_c$ to ensure a smooth density function, we can expect the particle-mesh method to scale better than the particle-particle method for large simulations. The PP method has a time complexity of $O(n_p^2)$ which grows much faster with the number of particles. In contrast, the primary bottleneck in the PM method is the FFT computation, which scales more modestly, thus it is more efficient for large-scale systems.



2.2 Application 2: Comparing Universes with Distributed Memory
./build/bin/NBody_Visualiser -nc 101 -np 10 -t 1.5 -dt 0.01 -F 1.0 -o Images/F1.0 -s 93170929
Images can be found in "/workspaces/COMP0210Assignment2/Images/F1.0"

./build/bin/NBody_Visualiser -nc 101 -np 10 -t 1.5 -dt 0.01 -F 1.02 -o Images/F1.02 -s 93170929
Images can be found in "/workspaces/COMP0210Assignment2/Images/F1.02"