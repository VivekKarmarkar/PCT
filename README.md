# PCT
Scripts relating to the Point Cluster Technique algorithm

# First Iteration Highlights

## Algorithms
1D Unconstrained Optimization (Original PCT, 3 <= N <= 12)<br />
2D Unconstrained Optimization (New mass redistribution function, 3 <= N <= 12)
2D Constrained Optimization (New mass redistribution function which also accounts for angular error, Performed on each eigenvalue separately, 3 <= N <= 12)
Iterative Perturbation Theory based approach (N=6, tolerances introduced by hand)

## Setup
Markers placed on cylindrical surface with varying radius (truncated cone)
Marker configuration created using a random number generator
Noise model based on maximum noise values from the noise model in the original PCT paper

## Metrics
Center of Mass location offset norm
Eigenvalue norm percentage offset
Eigenvector angular offset
Cluster Point Global Position reconstruction error
