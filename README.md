# PCT
Scripts relating to the Point Cluster Technique algorithm

# First Iteration Highlights

## Algorithms
1D Unconstrained Optimization (Original PCT, 3 <= N <= 12) <br />
2D Unconstrained Optimization (New mass redistribution function, 3 <= N <= 12) <br />
2D Constrained Optimization (New mass redistribution function which also accounts for angular error, Performed on each eigenvalue separately, 3 <= N <= 12) <br />
Only some of the Optimizations are Levenberg-Marquardt <br />
Iterative Perturbation Theory based approach (N=6, tolerances introduced by hand) <br />

## Setup
Markers placed on cylindrical surface with varying radius (truncated cone) <br />
Marker configuration created using a random number generator <br />
Noise model based on maximum noise values from the noise model in the original PCT paper <br />

## Metrics
Center of Mass location offset norm <br />
Eigenvalue norm percentage offset <br />
Eigenvector angular offset <br />
Cluster Point Global Position reconstruction error <br />
