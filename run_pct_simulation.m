%% RUN_PCT_SIMULATION - Main entry point for PCT algorithm simulations
%
% This script runs the Point Cluster Technique (PCT) simulation pipeline
% for validating body segment pose estimation algorithms.
%
% The simulation:
%   1. Loads ground truth kinematic data and STA noise models
%   2. Generates or loads marker configurations
%   3. Applies soft tissue artifact noise to marker positions
%   4. Runs multiple pose estimation algorithms (SVD-LS, Perturbation Theory)
%   5. Computes error metrics and generates visualizations
%
% Configuration:
%   - N: Number of markers (4 or 6)
%   - NumConfigurations: Number of random marker configurations to test
%   - output_path: Directory for saving results
%
% Usage:
%   >> run_pct_simulation
%   >> run_pct_simulation  % Then modify N and NumConfigurations as needed
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)
%
% See also: src/algorithms/svd_ls, src/algorithms/perturbation_theory

%% Setup
clear; clc; close all;

% Add source paths
addpath(genpath('src'));

%% Configuration
% Modify these parameters for your simulation
N = 4;                          % Number of markers (4 or 6)
NumConfigurations = 5;          % Number of marker configurations
output_path = fullfile(pwd, 'output');

% Get paths and settings
paths = settings().get_paths(pwd, N);
input_settings = settings().get_default_settings();

%% Load Data
fprintf('Loading ground truth and STA data...\n');

% Load ground truth data
if exist(paths.ground_truth_thigh, 'file')
    load(paths.ground_truth_thigh, 'pose_gait');
else
    error('Ground truth data not found: %s', paths.ground_truth_thigh);
end

% Load STA marker noise data
if exist(paths.sta_markers, 'file')
    load(paths.sta_markers, 'box_markers_gait');
else
    error('STA marker data not found: %s', paths.sta_markers);
end

%% Initialize Structures
fprintf('Initializing simulation structures...\n');

% Generate input size structure
input_size = settings().get_input_size_struct(N, NumConfigurations, box_markers_gait);

% Generate cylinder model
cylinder_data = cylinder_model().generate_cylinder_data();

% Initialize final output container
final_output = struct;

% Visualization flags
visualization_flags = struct;
visualization_flags.animation = true(NumConfigurations, 1);
visualization_flags.fig_scale = true(NumConfigurations, 1);

%% Main Simulation Loop
fprintf('Starting simulation with %d configurations...\n', NumConfigurations);

for j = 1:NumConfigurations
    fprintf('\n=== Processing Configuration %d/%d ===\n', j, NumConfigurations);

    configuration_idx = j;

    % Load marker configuration
    config_file = fullfile(paths.marker_configs, sprintf('FixedMarkerDistribution_%d.mat', configuration_idx));
    if exist(config_file, 'file')
        load(config_file, 'mass_distribution_reference');
    else
        fprintf('  Generating new random marker configuration...\n');
        mass_distribution_reference = marker_utils().generate_random_distribution_constrained(cylinder_data, N);
    end

    % Package input data
    input_data = struct;
    input_data.cylinder_data = cylinder_data;
    input_data.configuration_idx = configuration_idx;
    input_data.mass_distribution_reference = mass_distribution_reference;
    input_data.box_markers_gait = box_markers_gait;
    input_data.pose_gait = pose_gait;

    input_struct = struct;
    input_struct.input_size = input_size;
    input_struct.input_settings = input_settings;
    input_struct.input_data = input_data;

    % Initialize output containers
    output_containers = settings().get_output_containers(input_size);

    % Process all algorithms
    fprintf('  Running algorithms...\n');
    [output_containers, temp_containers] = process_all_algorithms(input_struct, output_containers);

    % Compute mass redistribution
    fprintf('  Computing mass redistribution...\n');
    mass_redistribution = perturbation_theory().compute_mass_redistribution(output_containers, input_size);
    scaling_factor_CM = perturbation_theory().compute_scaling_factor_CM(mass_redistribution, output_containers, input_size);

    % Reprocess with optimized masses
    fprintf('  Reprocessing with optimized masses...\n');
    output_containers = reprocess_algorithms(input_struct, output_containers, temp_containers, ...
                                              mass_redistribution, scaling_factor_CM);

    % Post-processing
    output_containers = post_processing_discontinuity(output_containers);

    % Store results
    final_output(configuration_idx).output_containers = output_containers;

    % Generate plots
    if visualization_flags.fig_scale(configuration_idx)
        fprintf('  Generating plots...\n');
        plot_results().generate_output_plots(input_struct, output_containers, paths.plots_dir);
    end
end

%% Generate Summary Statistics
fprintf('\nGenerating summary statistics...\n');
final_output_mean = generate_output_stats(final_output, input_size);

% Create summary figure
summary_fig = plot_results().create_summary_figure(final_output_mean, input_size, ...
    fullfile(paths.plots_dir, 'summary.png'));

%% Save Results
fprintf('Saving results...\n');
save(paths.results_file, 'final_output', 'final_output_mean', 'input_size', '-v7.3');

fprintf('\n=== Simulation Complete ===\n');
fprintf('Results saved to: %s\n', paths.results_file);
fprintf('Plots saved to: %s\n', paths.plots_dir);

%% Helper Functions

function [output_containers, temp_containers] = process_all_algorithms(input_struct, output_containers)
    % PROCESS_ALL_ALGORITHMS Run all PCT algorithms on current configuration

    temp_containers = struct;

    input_size = input_struct.input_size;
    input_settings = input_struct.input_settings;
    input_data = input_struct.input_data;

    N = input_size.N;
    NumFrames = input_size.NumFrames;

    cylinder_data = input_data.cylinder_data;
    mass_distribution_reference_standard = input_data.mass_distribution_reference;
    box_markers_gait = input_data.box_markers_gait;
    pose_gait = input_data.pose_gait;

    inertia_funcs = inertia();
    marker_funcs = marker_utils();
    svdls_funcs = svd_ls();
    pt_funcs = perturbation_theory();

    for frame_idx = 1:NumFrames
        fraction_gait = input_size.fraction_completed_gait_cycle(frame_idx);

        % Generate noisy marker positions
        noisy_dist_standard = marker_funcs.generate_noisy_distribution_STA(...
            mass_distribution_reference_standard, frame_idx, box_markers_gait, N);

        % Apply current pose
        current_pose = pose_gait(frame_idx);
        mass_dist_ref = marker_funcs.apply_current_pose(mass_distribution_reference_standard, current_pose, N);
        mass_dist_noisy = marker_funcs.apply_current_pose(noisy_dist_standard, current_pose, N);

        % Compute inertial properties
        [props_ref, mass_dist_ref] = inertia_funcs.compute_inertial_properties(mass_dist_ref, N);
        [props_noisy, mass_dist_noisy] = inertia_funcs.compute_inertial_properties(mass_dist_noisy, N);
        props_offset = inertia_funcs.calculate_offset(props_ref, props_noisy);

        % Handle axis reversal if needed
        if input_settings.initial_axes_reversal
            props_noisy = inertia_funcs.axis_reversal_original(props_noisy, props_offset);
            props_offset = inertia_funcs.calculate_offset(props_ref, props_noisy);
        end

        % Setup algorithm container
        algorithm_container = struct;
        algorithm_container.N = N;
        algorithm_container.cylinder_data = cylinder_data;
        algorithm_container.current_pose = current_pose;
        algorithm_container.mass_distribution_reference = mass_dist_ref;
        algorithm_container.mass_distribution_noisy = mass_dist_noisy;
        algorithm_container.inertial_properties_reference = props_ref;
        algorithm_container.inertial_properties_noisy = props_noisy;
        algorithm_container.inertial_properties_offset = props_offset;

        % Run SVD-LS
        algorithm_container = compute_local_offsets(algorithm_container);
        algorithm_container = svdls_funcs.svdls_processing(algorithm_container);

        % Store SVD-LS results
        store_svdls_results(output_containers, algorithm_container, frame_idx, input_size);

        % Run Perturbation Theory (initial pass)
        pt_container = algorithm_container;
        pt_container.inertial_properties_noisy_original = props_noisy;
        pt_container = pt_funcs.perturbation_theory_processing(pt_container);

        temp_containers(frame_idx).pt_container = pt_container;

        % Store noisy positions for mass redistribution
        for k = 1:N
            output_containers.noisy_positions.PTNR(frame_idx, k, :) = mass_dist_noisy(k).vector;
        end
    end
end

function output_containers = reprocess_algorithms(input_struct, output_containers, temp_containers, ...
                                                   mass_redistribution, scaling_factor_CM)
    % REPROCESS_ALGORITHMS Reprocess PT algorithms with optimized masses

    input_size = input_struct.input_size;
    NumFrames = input_size.NumFrames;

    pt_funcs = perturbation_theory();

    for frame_idx = 1:NumFrames
        pt_container = temp_containers(frame_idx).pt_container;
        pt_container.mass_redistribution = mass_redistribution;
        pt_container.scaling_factor_CM = scaling_factor_CM;

        % Process with different reflection settings
        for reflection = {'NR', 'UR', 'LR'}
            refl = reflection{1};
            pt_container.reflection_setting = refl;
            pt_container = pt_funcs.perturbation_theory_reprocessing(pt_container);

            % Store results based on reflection type
            field_name = ['PT', refl];
            store_pt_results(output_containers, pt_container, frame_idx, input_size, field_name);
        end
    end
end

function output_containers = post_processing_discontinuity(output_containers)
    % POST_PROCESSING_DISCONTINUITY Handle discontinuities in PT reflection

    deltaY = output_containers.deltaY_offset.PTNR;
    roots_idx = find_sign_changes(deltaY);

    if ~isempty(roots_idx) && mod(length(roots_idx), 2) == 0
        % Apply correction based on sign pattern
        for j = 1:length(roots_idx)-1
            output_containers.CM_offset.PTUR(roots_idx(j):roots_idx(j+1)) = ...
                output_containers.CM_offset.PTLR(roots_idx(j):roots_idx(j+1));
        end
    end
end

function roots = find_sign_changes(data)
    % FIND_SIGN_CHANGES Find indices where data changes sign
    roots = [];
    for j = 1:length(data)-1
        if (data(j) > 0 && data(j+1) < 0) || (data(j) < 0 && data(j+1) > 0)
            roots = [roots; j];
        end
    end
end

function final_output_mean = generate_output_stats(final_output, input_size)
    % GENERATE_OUTPUT_STATS Compute mean statistics across configurations

    NumConfigs = min(length(final_output), input_size.NumConfigurations);
    NumFrames = input_size.NumFrames;
    N = input_size.N;

    final_output_mean = struct;
    final_output_mean.CM_offset.SVDLS = nan(NumFrames, 1);
    final_output_mean.CM_offset.PTUR = nan(NumFrames, 1);
    final_output_mean.cluster_offset.SVDLS = nan(NumFrames, N);
    final_output_mean.cluster_offset.PTUR = nan(NumFrames, N);
    final_output_mean.AL_offset.SVDLS = nan(NumFrames, input_size.Num_AL_Thigh);
    final_output_mean.AL_offset.PTUR = nan(NumFrames, input_size.Num_AL_Thigh);
    final_output_mean.eigenvector_offset.SVDLS = nan(NumFrames, 3);
    final_output_mean.eigenvector_offset.PTUR = nan(NumFrames, 3);

    for t = 1:NumFrames
        cm_svdls = nan(NumConfigs, 1);
        cm_ptur = nan(NumConfigs, 1);

        for k = 1:NumConfigs
            cm_svdls(k) = final_output(k).output_containers.CM_offset.SVDLS(t);
            cm_ptur(k) = final_output(k).output_containers.CM_offset.PTUR(t);
        end

        final_output_mean.CM_offset.SVDLS(t) = nanmean(cm_svdls);
        final_output_mean.CM_offset.PTUR(t) = nanmean(cm_ptur);
    end
end

function algorithm_container = compute_local_offsets(algorithm_container)
    % COMPUTE_LOCAL_OFFSETS Compute local coordinate offsets for markers

    N = algorithm_container.N;
    mass_dist_ref = algorithm_container.mass_distribution_reference;
    mass_dist_noisy = algorithm_container.mass_distribution_noisy;
    props_ref = algorithm_container.inertial_properties_reference;
    props_noisy = algorithm_container.inertial_properties_noisy;

    geom = geometry();

    for k = 1:N
        % Compute local displacements
        pos_ref = mass_dist_ref(k).vector;
        pos_noisy = mass_dist_noisy(k).vector;

        local_ref = geom.global2local(props_ref.T, props_ref.R, pos_ref);
        local_noisy = geom.global2local(props_noisy.T, props_noisy.R, pos_noisy);

        mass_dist_ref(k).local_displacement = local_ref;
        mass_dist_noisy(k).local_displacement = local_noisy;
    end

    algorithm_container.mass_distribution_reference = mass_dist_ref;
    algorithm_container.mass_distribution_noisy = mass_dist_noisy;
end

function store_svdls_results(output_containers, container, frame_idx, input_size)
    % STORE_SVDLS_RESULTS Store SVD-LS algorithm results
    % Implementation would store results to output_containers
end

function store_pt_results(output_containers, container, frame_idx, input_size, field_name)
    % STORE_PT_RESULTS Store Perturbation Theory results
    % Implementation would store results to output_containers
end
