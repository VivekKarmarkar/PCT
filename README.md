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

# Second Iteration Highlights

## Algorithm Changes
Iterative Perturbation Theory based approach now supports N=4. Nonlinear mass solver using Levenberg-Marquardt for N=4 <br />
Iterative Perturbation Theory algorithm now constrains the CM location <br />
CM reflect function applied to the constrained Iterative Perturbation Theory CM as well as the Original PCT CM
Hybrid algorithm combining Iterative Perturbation Theory CM and optimized PCT eigenvectors introduced

## Setup Changes
Markers placed only on lateral and anterior aspects of the thigh. This has been modeled by restricted angular range <br />
Marker have been placed uniformly. This has been modeled by placing markers randomly in boxes of roughly equal area <br />
Markers are well separated and satisfy a minimum distance criterion for no inter-marker error propagation. This has been modeled by invalidating configurations not satisfying specific inter-marker distance constraints <br />
The following anatomical markers have been placed on the cylinder modeling the thigh: Lateral Epicondlye (LE), Medial Epicondlye (ME), Femoral Head (FH), Greater Trochanter (GT) <br />

## New Metrics
Markerless Anatomical Landmark reconsruction error <br />
Standard deviation in the global reconstruction errors of cluster point <br />
Planarity of marker configurations <br />
