%% SETTINGS - Configuration and settings for PCT simulations
%
% This module provides configuration structures and default settings
% for running PCT algorithm simulations.
%
% Functions:
%   get_default_settings      - Get default algorithm settings
%   get_input_size_struct     - Generate input size structure
%   get_output_containers     - Initialize output data containers
%   get_paths                 - Get configurable paths
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = settings()
    funcs.get_default_settings = @get_default_settings;
    funcs.get_input_size_struct = @get_input_size_struct;
    funcs.get_output_containers = @get_output_containers;
    funcs.get_paths = @get_paths;
    funcs.get_algorithm_list = @get_algorithm_list;
end

function input_settings = get_default_settings()
    % GET_DEFAULT_SETTINGS Get default algorithm configuration
    %
    % Returns a structure with default settings for PCT algorithms
    % including axis reversal flags and optimization parameters.
    %
    % Output:
    %   input_settings - Default settings structure

    input_settings = struct;

    % Axis reversal settings
    % Controls handling of eigenvector sign ambiguity
    input_settings.initial_axes_reversal = false;
    input_settings.initial_axes_perturbed_reversal = true;
    input_settings.optimized_axes_reversal = true;

    % Mass redistribution settings
    input_settings.m0_switch_idx = 2;

    % Tolerance settings
    input_settings.MOI_norm_tolerance = 0.1;
    input_settings.CM_distance_threshold = 4;          % cm
    input_settings.plane_distance_threshold = 1;       % cm
    input_settings.height_diff_threshold = 15;         % cm
    input_settings.angular_diff_threshold = 60;        % degrees

    % Optimization settings
    input_settings.max_iterations = 10000;
    input_settings.function_tolerance = 1.0;
    input_settings.optimizer = 'levenberg-marquardt';

    % Gait cycle parameters
    input_settings.toe_off_fraction = 0.6;
end

function input_size = get_input_size_struct(N, num_configurations, box_markers_gait)
    % GET_INPUT_SIZE_STRUCT Generate input size structure
    %
    % Creates a structure containing all size-related parameters
    % for the simulation.
    %
    % Inputs:
    %   N                  - Number of markers
    %   num_configurations - Number of marker configurations to test
    %   box_markers_gait   - Gait cycle marker data (for frame info)
    %
    % Output:
    %   input_size - Size parameters structure

    input_size = struct;

    % Extract gait cycle information
    fraction_completed_gait_cycle = box_markers_gait(1).fraction_completed_gait_cycle;
    num_frames = length(fraction_completed_gait_cycle);

    % Populate structure
    input_size.N = N;
    input_size.fraction_completed_gait_cycle = fraction_completed_gait_cycle;
    input_size.NumFrames = num_frames;
    input_size.NumConfigurations = num_configurations;

    % Anatomical parameters
    input_size.Num_AL_Thigh = 4;      % Lateral/Medial Epicondyle, GT, FH
    input_size.Num_AL_Shank = 4;      % Similar for shank
    input_size.NumAngles = 3;         % 3 principal axes
    input_size.NumEigenvectors = 3;
end

function output_containers = get_output_containers(input_size)
    % GET_OUTPUT_CONTAINERS Initialize output data containers
    %
    % Creates pre-allocated structures for storing algorithm outputs
    % across all frames and configurations.
    %
    % Input:
    %   input_size - Size parameters structure
    %
    % Output:
    %   output_containers - Initialized output structure

    output_containers = struct;

    % Algorithm list
    algorithm_list = {'SVDLS', 'PCTUO', 'PTNR', 'PTUR', 'PTLR'};
    num_algorithms = length(algorithm_list);

    % Dimensions
    num_frames = input_size.NumFrames;
    num_AL = input_size.Num_AL_Thigh;
    N = input_size.N;
    num_eigenvectors = 3;

    % Initialize containers for each metric
    eigenvalueNorm_offset = struct;
    eigenvector_offset = struct;
    eps_star = struct;
    mass_redistribution = struct;
    CM_offset = struct;
    cluster_offset = struct;
    AL_offset = struct;
    deltaY_offset = struct;
    T_estimated = struct;
    R_estimated = struct;
    T_reference = struct;
    R_reference = struct;
    AL_estimated = struct;
    AL_reference = struct;
    AF_Translation_estimated = struct;
    AF_Rotation_estimated = struct;
    AF_Translation_reference = struct;
    AF_Rotation_reference = struct;
    inertiaTensor = struct;
    noisy_positions = struct;

    % Pre-allocate for each algorithm
    for n = 1:num_algorithms
        alg = algorithm_list{n};

        eigenvalueNorm_offset.(alg) = nan(num_frames, 1);
        eigenvector_offset.(alg) = nan(num_frames, num_eigenvectors);
        eps_star.(alg) = nan(num_frames, 1);
        mass_redistribution.(alg) = nan(num_frames, N);
        CM_offset.(alg) = nan(num_frames, 1);
        cluster_offset.(alg) = nan(num_frames, N);
        AL_offset.(alg) = nan(num_frames, num_AL);
        deltaY_offset.(alg) = nan(num_frames, 1);

        T_estimated.(alg) = nan(num_frames, 3);
        R_estimated.(alg) = nan(num_frames, 3, 3);
        T_reference.(alg) = nan(num_frames, 3);
        R_reference.(alg) = nan(num_frames, 3, 3);

        AL_estimated.(alg) = nan(num_frames, num_AL, 3);
        AL_reference.(alg) = nan(num_frames, num_AL, 3);

        AF_Translation_estimated.(alg) = nan(num_frames, 3);
        AF_Rotation_estimated.(alg) = nan(num_frames, 3, 3);
        AF_Translation_reference.(alg) = nan(num_frames, 3);
        AF_Rotation_reference.(alg) = nan(num_frames, 3, 3);

        inertiaTensor.(alg) = nan(num_frames, 3, 3);
        noisy_positions.(alg) = nan(num_frames, N, 3);
    end

    % Assign to output structure
    output_containers.eigenvalueNorm_offset = eigenvalueNorm_offset;
    output_containers.eigenvector_offset = eigenvector_offset;
    output_containers.eps_star = eps_star;
    output_containers.mass_redistribution = mass_redistribution;
    output_containers.CM_offset = CM_offset;
    output_containers.cluster_offset = cluster_offset;
    output_containers.AL_offset = AL_offset;
    output_containers.deltaY_offset = deltaY_offset;
    output_containers.T_estimated = T_estimated;
    output_containers.R_estimated = R_estimated;
    output_containers.T_reference = T_reference;
    output_containers.R_reference = R_reference;
    output_containers.AL_estimated = AL_estimated;
    output_containers.AL_reference = AL_reference;
    output_containers.AF_Translation_estimated = AF_Translation_estimated;
    output_containers.AF_Rotation_estimated = AF_Rotation_estimated;
    output_containers.AF_Translation_reference = AF_Translation_reference;
    output_containers.AF_Rotation_reference = AF_Rotation_reference;
    output_containers.inertiaTensor = inertiaTensor;
    output_containers.noisy_positions = noisy_positions;
end

function paths = get_paths(base_path, N)
    % GET_PATHS Get configurable paths for data and output
    %
    % Returns a structure with all file paths configured relative
    % to a base path. Allows easy configuration for different systems.
    %
    % Inputs:
    %   base_path - Base directory for PCT project
    %   N         - Number of markers (for path construction)
    %
    % Output:
    %   paths - Structure with all file paths

    if nargin < 1 || isempty(base_path)
        base_path = pwd;
    end
    if nargin < 2
        N = 4;
    end

    paths = struct;

    % Data paths
    paths.data_dir = fullfile(base_path, 'data');
    paths.ground_truth_thigh = fullfile(paths.data_dir, 'ground_truth_data_thigh.mat');
    paths.ground_truth_shank = fullfile(paths.data_dir, 'ground_truth_data_shank.mat');
    paths.sta_markers = fullfile(paths.data_dir, sprintf('box_markers_gait_%d.mat', N));
    paths.sta_markers_shank = fullfile(paths.data_dir, sprintf('box_markers_gait_shank_%d.mat', N));

    % Configuration paths
    paths.config_dir = fullfile(base_path, 'configurations');
    paths.marker_configs = fullfile(paths.config_dir, sprintf('MarkerConfigN%d', N));

    % Output paths
    paths.output_dir = fullfile(base_path, 'output');
    paths.plots_dir = fullfile(paths.output_dir, 'plots', sprintf('N=%d', N));
    paths.results_file = fullfile(paths.output_dir, sprintf('STA_output_N%d.mat', N));

    % Create directories if they don't exist
    dirs_to_create = {paths.data_dir, paths.config_dir, paths.output_dir, paths.plots_dir};
    for i = 1:length(dirs_to_create)
        if ~exist(dirs_to_create{i}, 'dir')
            mkdir(dirs_to_create{i});
        end
    end
end

function algorithms = get_algorithm_list(include_all)
    % GET_ALGORITHM_LIST Get list of available algorithms
    %
    % Returns cell array of algorithm identifiers.
    %
    % Input:
    %   include_all - Boolean to include all algorithms (default: false)
    %
    % Output:
    %   algorithms - Cell array of algorithm names

    if nargin < 1
        include_all = false;
    end

    if include_all
        algorithms = {'SVDLS', 'PCTUO', 'PTNR', 'PTUR', 'PTLR', 'PCTO'};
    else
        % Standard set for comparison
        algorithms = {'SVDLS', 'PTUR'};
    end
end
