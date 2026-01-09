%% GENERATE_MARKER_CONFIGURATIONS - Generate valid marker configurations
%
% This script generates random marker configurations that satisfy
% the anatomical and spacing constraints required for PCT analysis.
%
% Constraints enforced:
%   - Markers placed only on anterior and lateral aspects
%   - Minimum inter-marker height separation (15 cm)
%   - Minimum angular separation between marker pairs (60 degrees)
%   - Markers distributed across height zones
%
% Configuration:
%   - N: Number of markers (4 or 6)
%   - NumTrials: Number of random configurations to attempt
%   - output_dir: Directory for saving valid configurations
%
% Usage:
%   >> generate_marker_configurations
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

%% Setup
clear; clc; close all;

% Add source paths
addpath(genpath('src'));

%% Configuration
N = 4;                          % Number of markers (4 or 6)
NumTrials = 300000;             % Number of random trials
output_dir = fullfile(pwd, 'configurations', sprintf('MarkerConfigN%d', N));

% Create output directory
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% Generate Cylinder Model
fprintf('Generating cylinder model...\n');
cylinder_data = cylinder_model().generate_cylinder_data();

% Get marker utilities
marker_funcs = marker_utils();

%% Generate Configurations
fprintf('Generating marker configurations...\n');
fprintf('Target: Find valid configurations from %d trials\n', NumTrials);
fprintf('Markers: N = %d\n\n', N);

configuration_count = 0;

for j = 1:NumTrials
    % Progress update
    if mod(j, 10000) == 0
        fprintf('Trial %d/%d - Found %d valid configurations\n', j, NumTrials, configuration_count);
    end

    % Generate random distribution
    mass_distribution_reference = marker_funcs.generate_random_distribution_constrained(cylinder_data, N);

    % Check inter-marker distance constraints
    is_valid = marker_funcs.check_intermarker_distance(mass_distribution_reference, N);

    if is_valid
        configuration_count = configuration_count + 1;

        % Save configuration
        filename = fullfile(output_dir, sprintf('FixedMarkerDistribution_%d.mat', configuration_count));
        save(filename, 'mass_distribution_reference', '-v7.3');

        if mod(configuration_count, 100) == 0
            fprintf('  Saved configuration %d\n', configuration_count);
        end
    end
end

%% Summary
fprintf('\n=== Generation Complete ===\n');
fprintf('Total trials: %d\n', NumTrials);
fprintf('Valid configurations found: %d\n', configuration_count);
fprintf('Success rate: %.2f%%\n', 100 * configuration_count / NumTrials);
fprintf('Configurations saved to: %s\n', output_dir);

%% Visualize Sample Configuration
if configuration_count > 0
    fprintf('\nVisualizing sample configuration...\n');

    % Load first configuration
    load(fullfile(output_dir, 'FixedMarkerDistribution_1.mat'), 'mass_distribution_reference');

    % Plot
    figure('Name', 'Sample Marker Configuration');

    % Extract marker positions
    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;

    % Plot cylinder surface
    surf(X, Y, Z, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'FaceColor', [0.8, 0.8, 0.8]);
    hold on;

    % Plot markers
    colors = lines(N);
    for k = 1:N
        scatter3(mass_distribution_reference(k).X, ...
                 mass_distribution_reference(k).Y, ...
                 mass_distribution_reference(k).Z, ...
                 100, colors(k, :), 'filled', ...
                 'DisplayName', sprintf('Marker %d', k));
    end

    xlabel('X (cm)');
    ylabel('Y (cm)');
    zlabel('Z (cm)');
    title(sprintf('Sample Marker Configuration (N=%d)', N));
    legend('Location', 'bestoutside');
    axis equal;
    grid on;
    view(3);

    % Save figure
    saveas(gcf, fullfile(output_dir, 'sample_configuration.png'));
end
