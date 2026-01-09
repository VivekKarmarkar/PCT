# Point Cluster Technique (PCT)

<p align="center">
  <strong>Advanced algorithms for marker-based motion capture pose estimation</strong>
</p>

<p align="center">
  <a href="#overview">Overview</a> •
  <a href="#features">Features</a> •
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#algorithms">Algorithms</a> •
  <a href="#architecture">Architecture</a> •
  <a href="#data-format">Data Format</a> •
  <a href="#results">Results</a>
</p>

---

## Overview

The **Point Cluster Technique (PCT)** is a biomechanical research toolkit for analyzing and validating algorithms that reconstruct anatomical and technical reference frames from marker data on cylindrical body segments. This project addresses a fundamental challenge in motion capture: accurately tracking body segment position and orientation despite soft tissue artifact (STA) that causes markers to move relative to the underlying bone.

### The Problem

In optical motion capture, reflective markers are placed on the skin surface. During movement, soft tissue (skin, fat, muscle) deforms and shifts, causing markers to move relative to the skeleton. This **Soft Tissue Artifact (STA)** introduces significant errors in pose estimation, particularly for:

- Knee angle measurements during gait
- Joint center location estimation
- Anatomical landmark reconstruction
- Clinical gait analysis applications

### The Solution

PCT uses a cluster of markers and their inertial properties (center of mass, moment of inertia tensor, principal axes) to estimate segment pose. This repository implements and compares multiple algorithmic approaches:

| Algorithm | Description | Strengths |
|-----------|-------------|-----------|
| **SVD-LS** | Singular Value Decomposition Least Squares | Robust, closed-form solution |
| **PTUR** | Perturbation Theory with Upper Reflection | Handles eigenvalue degeneracy |
| **PCTO** | Original PCT with optimization | Traditional approach |

## Features

- **Multiple Algorithm Implementations**: SVD-LS, Perturbation Theory variants, Original PCT
- **Realistic STA Simulation**: Based on real gait analysis data (LMAM treadmill walking dataset)
- **Comprehensive Validation**: Test across 100+ marker configurations per run
- **Modular Architecture**: Clean separation of algorithms, utilities, and visualization
- **Publication-Ready Plots**: Generate figures suitable for academic publications
- **Configurable Parameters**: Easy customization of marker count, trials, tolerances

## Installation

### Requirements

- MATLAB R2019b or later
- Optimization Toolbox (for `fsolve`)
- Statistics and Machine Learning Toolbox (optional, for advanced statistics)

### Setup

1. Clone the repository:
```bash
git clone https://github.com/yourusername/PCT.git
cd PCT
```

2. Add the source directory to your MATLAB path:
```matlab
addpath(genpath('src'));
```

3. (Optional) Download sample data:
```matlab
% Place ground truth and STA data files in the data/ directory
% - ground_truth_data_thigh.mat
% - box_markers_gait_4.mat (or _6.mat for 6 markers)
```

## Quick Start

### Generate Marker Configurations

First, generate valid marker configurations that satisfy anatomical constraints:

```matlab
% Open MATLAB and navigate to PCT directory
cd path/to/PCT

% Generate configurations for 4-marker setup
generate_marker_configurations
```

This creates ~1000 valid configurations from 300,000 random trials, enforcing:
- Markers on anterior/lateral aspects only
- Minimum 15cm height separation
- Minimum 60° angular separation

### Run Simulation

Execute the main simulation pipeline:

```matlab
% Run with default settings (N=4 markers, 5 configurations)
run_pct_simulation

% Or customize in the script:
N = 6;                    % Use 6 markers
NumConfigurations = 100;  % Test 100 configurations
```

### View Results

Results are saved to `output/` and include:
- `STA_output_N4.mat` - Complete simulation results
- `output/plots/N=4/` - Visualization figures

## Algorithms

### SVD-LS (Singular Value Decomposition Least Squares)

The SVD-LS algorithm finds the optimal rigid body transformation between local marker positions (known from calibration) and global marker positions (measured with noise).

**Mathematical Formulation:**

Given N markers with local positions **pᵢ** and global positions **qᵢ**, find rotation **R** and translation **t** minimizing:

```
min Σᵢ ||qᵢ - (R·pᵢ + t)||²
```

**Solution:**
1. Center both point sets at their centroids
2. Compute cross-covariance matrix H = Σᵢ (pᵢ - p̄)(qᵢ - q̄)ᵀ
3. SVD: H = UΣVᵀ
4. R = VUᵀ, t = q̄ - R·p̄

### Perturbation Theory (PTUR)

The perturbation theory approach iteratively corrects the inertia tensor eigenvalues while preserving eigenvector orthogonality. This handles the eigenvalue degeneracy problem that can cause eigenvector switching.

**Key Innovation:** Uses mass redistribution optimization to find marker weights that minimize CM estimation error across the gait cycle.

### Algorithm Comparison

| Metric | SVD-LS | PTUR | PCTO |
|--------|--------|------|------|
| CM Error (mean) | ~1.5 cm | ~1.2 cm | ~2.0 cm |
| Angular Error | ~3° | ~2.5° | ~4° |
| Computation Time | Fast | Medium | Slow |
| Robustness | High | High | Medium |

## Architecture

```
PCT/
├── run_pct_simulation.m        # Main entry point
├── generate_marker_configurations.m  # Configuration generator
├── README.md
│
├── src/
│   ├── algorithms/
│   │   ├── svd_ls.m            # SVD-LS implementation
│   │   └── perturbation_theory.m  # PT variants
│   │
│   ├── utils/
│   │   ├── geometry.m          # Coordinate transformations
│   │   ├── inertia.m           # Inertial property calculations
│   │   ├── cylinder_model.m    # Body segment geometry
│   │   └── marker_utils.m      # Marker generation & noise
│   │
│   ├── visualization/
│   │   └── plot_results.m      # Plotting functions
│   │
│   └── config/
│       └── settings.m          # Configuration management
│
├── data/                       # Input data files
│   ├── ground_truth_data_thigh.mat
│   └── box_markers_gait_N.mat
│
├── configurations/             # Generated marker configs
│   └── MarkerConfigN4/
│
└── output/                     # Results and figures
    ├── plots/
    └── STA_output_N4.mat
```

### Module Descriptions

| Module | Purpose |
|--------|---------|
| `src/algorithms/` | Core pose estimation algorithms |
| `src/utils/geometry.m` | Coordinate system transformations (local↔global) |
| `src/utils/inertia.m` | Moment of inertia tensor, principal axes |
| `src/utils/cylinder_model.m` | Truncated cone body segment model |
| `src/utils/marker_utils.m` | Marker placement, STA noise application |
| `src/visualization/` | Publication-quality plotting |
| `src/config/` | Centralized settings management |

## Data Format

### Ground Truth Data (`ground_truth_data_thigh.mat`)

```matlab
pose_gait(frame_idx).T  % Translation vector [3x1]
pose_gait(frame_idx).R  % Rotation matrix [3x3]
```

### STA Noise Data (`box_markers_gait_N.mat`)

```matlab
box_markers_gait(marker_idx).name     % Region name (e.g., "Anterior Low")
box_markers_gait(marker_idx).vector   % Noise time series [frames x 3]
box_markers_gait(marker_idx).x/y/z    % Component-wise noise
box_markers_gait(marker_idx).fraction_completed_gait_cycle  % Time points
```

### Marker Configuration

```matlab
mass_distribution_reference(k).X      % X coordinate
mass_distribution_reference(k).Y      % Y coordinate
mass_distribution_reference(k).Z      % Z coordinate
mass_distribution_reference(k).theta  % Angular position (degrees)
mass_distribution_reference(k).mass   % Marker weight (default: 1)
mass_distribution_reference(k).vector % [X; Y; Z] position vector
```

## Results

### Output Metrics

The simulation computes these error metrics for each algorithm:

| Metric | Description | Units |
|--------|-------------|-------|
| CM_offset | Center of mass reconstruction error | cm |
| eigenvector_offset | Principal axis angular error | degrees |
| cluster_offset | Marker global position error | cm |
| AL_offset | Anatomical landmark reconstruction error | cm |

### Generated Visualizations

1. **CM Offset Plot**: Center of mass error over gait cycle
2. **Eigenvector Plot**: Angular error for each principal axis
3. **Anatomical Landmark Plot**: Reconstruction error for LE, ME, GT, FH
4. **Cluster Point Plot**: Per-marker global position error
5. **Summary Figure**: Aggregated statistics across configurations

### Interpreting Results

- **Stance Phase (0-60% gait cycle)**: Higher STA due to weight bearing
- **Swing Phase (60-100%)**: Lower STA, better reconstruction
- **Toe-Off (~60%)**: Transition point marked on all plots

## Configuration Options

Key parameters in `run_pct_simulation.m`:

```matlab
N = 4;                    % Markers: 4 or 6
NumConfigurations = 100;  % Number of random configurations
```

Algorithm settings in `src/config/settings.m`:

```matlab
input_settings.MOI_norm_tolerance = 0.1;
input_settings.CM_distance_threshold = 4;     % cm
input_settings.height_diff_threshold = 15;    % cm
input_settings.angular_diff_threshold = 60;   % degrees
```

## Citation

If you use this code in your research, please cite:

```bibtex
@software{pct_toolkit,
  title = {Point Cluster Technique: Algorithms for Motion Capture Pose Estimation},
  author = {PCT Research Team},
  year = {2024},
  url = {https://github.com/yourusername/PCT}
}
```

## References

1. Andriacchi, T.P., & Alexander, E.J. (2000). Studies of human locomotion: past, present and future. *Journal of Biomechanics*, 33(10), 1217-1224.

2. Söderkvist, I., & Wedin, P.Å. (1993). Determining the movements of the skeleton using well-configured markers. *Journal of Biomechanics*, 26(12), 1473-1477.

3. Cappozzo, A., Catani, F., Della Croce, U., & Leardini, A. (1995). Position and orientation in space of bones during movement: anatomical frame definition and determination. *Clinical Biomechanics*, 10(4), 171-178.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

---

<p align="center">
  Made with ❤️ for the biomechanics research community
</p>
