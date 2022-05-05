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
Iterative Perturbation Theory based approach now supports N=4. Nonlinear mass solver uses Levenberg-Marquardt for N=4 <br />
Iterative Perturbation Theory algorithm now constrains the CM location <br />
CM reflect function applied to the constrained Iterative Perturbation Theory CM as well as the Original PCT CM <br />
Hybrid algorithm combining Iterative Perturbation Theory CM and optimized PCT eigenvectors introduced <br />

## Setup Changes
Markers placed only on lateral and anterior aspects of the thigh. This has been modeled by restricted angular range <br />
Marker have been placed uniformly. This has been modeled by placing markers randomly in boxes of roughly equal area <br />
Markers are well separated and satisfy a minimum distance criterion for no inter-marker error propagation. This has been modeled by invalidating configurations not satisfying specific inter-marker distance constraints <br />
The following anatomical markers have been placed on the cylinder modeling the thigh: Lateral Epicondlye (LE), Medial Epicondlye (ME), Femoral Head (FH), Greater Trochanter (GT) <br />

## New Metrics
Markerless Anatomical Landmark reconsruction error <br />
Standard deviation in the global reconstruction errors of cluster point <br />
Planarity of marker configurations <br />

# Third Iteration Highlights

## Algorithm Changes
SVD-LS algorithm added <br />
CM reflect function not applicaple for SVD-LS since SVD-LS predicts noisy PCT CM <br />

## Setup Changes
None <br />

## New Metrics
None <br />

# Fourth Iteration Highlights

## Algorithm Changes
2D Optimization permanently removed <br />
1D Unconstrained Optimization (Original PCT and PCT reflected) temporarily removed due to instabilities <br />
Iterative initial condition update applied for Perturbation Theory (PT), m0(t) = m*(t-dt) <br />
CM_reflect adapted for cylinder in non-standard pose <br />
Post-processing of discontinuity in PT reflection curves included <br />

## Setup Changes
Real Soft Tissue Artifact (STA) noise added to markers <br />
PCT noise model removed <br />
Entire time series simulated <br />
System driven by kinematic input generated from real Groud Truth Data (GTD) <br />

## New Metrics
### For each configuration
Center of Mass (CM) reconstruction error <br />
Eigenvector reconstruction error <br />
Cluster Point reconstruction error <br />
Anatomical Landmark (AL) reconstruction error <br />
Technical Frame (TF) pose <br />
Anatomical Frame (AF) pose <br/>

### For the average configuration
Center of Mass (CM) reconstruction error <br />
Eigenvector reconstruction error <br />
Cluster Point reconstruction error <br />
Anatomical Landmark (AL) reconstruction error <br />

# Fifth Iteration Highlights

## Algorithm Changes
Bug in nonlinear mass solver resolved <br />
Best fit mass redistribution computed instead of solving for the mass resdistribution at each time step. This is done to avoid instabilities <br />
CM_constraint fuction modified based on vector geometry <br />
CM_constraint function no longer embedded in the nonlinear mass solver <br />

# Sixth Iteration Highlights

## Algorithm Changes
Estimate Rotation Matrix using the known good Translation vector <br />
