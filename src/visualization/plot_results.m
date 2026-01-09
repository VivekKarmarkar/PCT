%% PLOT_RESULTS - Visualization functions for PCT results
%
% This module provides comprehensive plotting functions for
% visualizing PCT algorithm results and performance metrics.
%
% Functions:
%   plot_cm_offset           - Plot CM reconstruction error over gait cycle
%   plot_eigenvector_offset  - Plot eigenvector angular errors
%   plot_anatomical_landmarks - Plot AL reconstruction errors
%   plot_cluster_points      - Plot cluster point errors
%   plot_knee_angle          - Plot estimated vs reference knee angle
%   assign_color_info        - Get algorithm color scheme
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = plot_results()
    funcs.plot_cm_offset = @plot_cm_offset;
    funcs.plot_eigenvector_offset = @plot_eigenvector_offset;
    funcs.plot_anatomical_landmarks = @plot_anatomical_landmarks;
    funcs.plot_cluster_points = @plot_cluster_points;
    funcs.plot_knee_angle = @plot_knee_angle;
    funcs.assign_color_info = @assign_color_info;
    funcs.generate_output_plots = @generate_output_plots;
    funcs.create_summary_figure = @create_summary_figure;
end

function color_info = assign_color_info(use_default)
    % ASSIGN_COLOR_INFO Get algorithm color scheme
    %
    % Returns a structure with colors for each algorithm
    % for consistent visualization.
    %
    % Input:
    %   use_default - Boolean to use default MATLAB colors
    %
    % Output:
    %   color_info - Struct with color values for each algorithm

    color_info = struct;
    color_info.reference = 'k';

    if use_default
        color_info.SVDLS = [1, 0, 1];       % Magenta
        color_info.PTUR = [0, 0, 1];        % Blue
        color_info.PCTO = [1, 0, 0];        % Red
    else
        % Custom publication colors
        color_info.SVDLS = [6/255, 199/255, 151/255];    % Teal
        color_info.PTUR = [0, 114/255, 189/255];         % Blue
        color_info.PCTO = [162/255, 20/255, 47/255];     % Dark red
    end
end

function fig = plot_cm_offset(output_containers, input_size, config_idx, save_path)
    % PLOT_CM_OFFSET Plot Center of Mass reconstruction error
    %
    % Creates a plot showing CM offset over the gait cycle for
    % multiple algorithms.
    %
    % Inputs:
    %   output_containers - Algorithm output data
    %   input_size        - Size parameters including gait cycle info
    %   config_idx        - Configuration index for title
    %   save_path         - Optional path to save figure
    %
    % Output:
    %   fig - Figure handle

    fraction_gait = input_size.fraction_completed_gait_cycle;
    TOE_OFF = 0.6;

    color_info = assign_color_info(false);
    algorithms = {'SVDLS', 'PTUR'};

    fig = figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);

    % Toe-off line
    xline(TOE_OFF, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Toe Off');
    hold on;

    % Plot each algorithm
    for i = 1:length(algorithms)
        alg = algorithms{i};
        data = output_containers.CM_offset.(alg);
        plot(fraction_gait, data, 'Color', color_info.(alg), ...
             'LineWidth', 1.5, 'DisplayName', alg);
    end

    % Formatting
    grid on;
    legend('Location', 'best');
    xlabel('Fraction of Gait Cycle', 'FontSize', 12);
    ylabel('CM Offset (cm)', 'FontSize', 12);
    title(sprintf('Center of Mass Reconstruction Error - Config %d', config_idx), ...
          'FontSize', 14);

    % Save if path provided
    if nargin >= 4 && ~isempty(save_path)
        saveas(fig, save_path);
    end
end

function fig = plot_eigenvector_offset(output_containers, input_size, config_idx, save_path)
    % PLOT_EIGENVECTOR_OFFSET Plot eigenvector angular errors
    %
    % Creates subplots showing angular error for each principal axis
    % over the gait cycle.

    fraction_gait = input_size.fraction_completed_gait_cycle;
    TOE_OFF = 0.6;

    color_info = assign_color_info(false);
    algorithms = {'SVDLS', 'PTUR'};

    fig = figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);

    for j = 1:3
        subplot(3, 1, j);

        xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
        hold on;

        for i = 1:length(algorithms)
            alg = algorithms{i};
            data = output_containers.eigenvector_offset.(alg)(:, j);
            plot(fraction_gait, data, 'Color', color_info.(alg), ...
                 'LineWidth', 1.5, 'DisplayName', alg);
        end

        grid on;
        ylabel('Error (degrees)', 'FontSize', 11);
        title(sprintf('Principal Axis %d', j), 'FontSize', 12);

        if j == 1
            legend('Location', 'best');
        end
        if j == 3
            xlabel('Fraction of Gait Cycle', 'FontSize', 12);
        end
    end

    sgtitle(sprintf('Eigenvector Angular Offset - Config %d', config_idx), ...
            'FontSize', 14);

    if nargin >= 4 && ~isempty(save_path)
        saveas(fig, save_path);
    end
end

function fig = plot_anatomical_landmarks(output_containers, input_size, config_idx, save_path)
    % PLOT_ANATOMICAL_LANDMARKS Plot anatomical landmark errors
    %
    % Creates subplots for each anatomical landmark showing
    % reconstruction error over the gait cycle.

    fraction_gait = input_size.fraction_completed_gait_cycle;
    TOE_OFF = 0.6;

    color_info = assign_color_info(false);
    algorithms = {'SVDLS', 'PTUR'};

    AL_NAMES = {'Lateral Epicondyle', 'Medial Epicondyle', ...
                'Greater Trochanter', 'Femoral Head'};

    fig = figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);

    for j = 1:4
        subplot(2, 2, j);

        xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
        hold on;

        for i = 1:length(algorithms)
            alg = algorithms{i};
            data = output_containers.AL_offset.(alg)(:, j);
            plot(fraction_gait, data, 'Color', color_info.(alg), ...
                 'LineWidth', 1.5, 'DisplayName', alg);
        end

        grid on;
        title(AL_NAMES{j}, 'FontSize', 12);

        if j == 2
            legend('Location', 'best');
        end
        if j > 2
            xlabel('Fraction of Gait Cycle', 'FontSize', 11);
        end
        if mod(j, 2) == 1
            ylabel('Error (cm)', 'FontSize', 11);
        end
    end

    sgtitle(sprintf('Anatomical Landmark Reconstruction - Config %d', config_idx), ...
            'FontSize', 14);

    if nargin >= 4 && ~isempty(save_path)
        saveas(fig, save_path);
    end
end

function fig = plot_cluster_points(output_containers, input_size, config_idx, save_path)
    % PLOT_CLUSTER_POINTS Plot cluster point reconstruction errors
    %
    % Creates subplots for each marker showing global position
    % reconstruction error.

    fraction_gait = input_size.fraction_completed_gait_cycle;
    TOE_OFF = 0.6;
    N = input_size.N;

    color_info = assign_color_info(false);
    algorithms = {'SVDLS', 'PTUR'};

    fig = figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);

    for j = 1:N
        subplot(2, ceil(N/2), j);

        xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
        hold on;

        for i = 1:length(algorithms)
            alg = algorithms{i};
            data = output_containers.cluster_offset.(alg)(:, j);
            plot(fraction_gait, data, 'Color', color_info.(alg), ...
                 'LineWidth', 1.5, 'DisplayName', alg);
        end

        grid on;
        title(sprintf('Marker %d', j), 'FontSize', 12);

        if j == ceil(N/2)
            legend('Location', 'best');
        end
        if j > ceil(N/2)
            xlabel('Fraction of Gait Cycle', 'FontSize', 11);
        end
        if mod(j, ceil(N/2)) == 1 || j == 1
            ylabel('Error (cm)', 'FontSize', 11);
        end
    end

    sgtitle(sprintf('Cluster Point Reconstruction - Config %d', config_idx), ...
            'FontSize', 14);

    if nargin >= 4 && ~isempty(save_path)
        saveas(fig, save_path);
    end
end

function fig = plot_knee_angle(knee_angle_data, input_size, config_idx, save_path)
    % PLOT_KNEE_ANGLE Plot estimated vs reference knee angle
    %
    % Compares knee angle estimates from different algorithms
    % against the reference values.

    fraction_gait = input_size.fraction_completed_gait_cycle;
    TOE_OFF = 0.6;

    color_info = assign_color_info(false);

    ANGLE_NAMES = {'Flexion', 'Adduction', 'External Rotation'};

    fig = figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);

    for j = 1:3
        subplot(3, 1, j);

        xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
        hold on;

        % Reference
        plot(fraction_gait, knee_angle_data.reference(:, j), ...
             'k', 'LineWidth', 2, 'DisplayName', 'Reference');

        % SVDLS estimate
        if isfield(knee_angle_data.estimated, 'SVDLS')
            plot(fraction_gait, knee_angle_data.estimated.SVDLS(:, j), ...
                 'Color', color_info.SVDLS, 'LineWidth', 1.5, 'DisplayName', 'SVDLS');
        end

        % PTUR estimate
        if isfield(knee_angle_data.estimated, 'PTUR')
            plot(fraction_gait, knee_angle_data.estimated.PTUR(:, j), ...
                 'Color', color_info.PTUR, 'LineWidth', 1.5, 'DisplayName', 'PTUR');
        end

        grid on;
        ylabel([ANGLE_NAMES{j}, ' (deg)'], 'FontSize', 11);

        if j == 1
            legend('Location', 'best');
        end
        if j == 3
            xlabel('Fraction of Gait Cycle', 'FontSize', 12);
        end
    end

    sgtitle(sprintf('Knee Angle Estimation - Config %d', config_idx), ...
            'FontSize', 14);

    if nargin >= 4 && ~isempty(save_path)
        saveas(fig, save_path);
    end
end

function generate_output_plots(input_struct, output_containers, output_path)
    % GENERATE_OUTPUT_PLOTS Generate all standard output plots
    %
    % Creates and saves all standard visualization plots for a
    % single configuration.
    %
    % Inputs:
    %   input_struct      - Input parameters
    %   output_containers - Algorithm outputs
    %   output_path       - Base path for saving figures

    input_size = input_struct.input_size;
    config_idx = input_struct.input_data.configuration_idx;

    if nargin < 3
        output_path = pwd;
    end

    % Create output directories if needed
    if ~exist(fullfile(output_path, 'CM_plots'), 'dir')
        mkdir(fullfile(output_path, 'CM_plots'));
    end
    if ~exist(fullfile(output_path, 'Eigenvector_plots'), 'dir')
        mkdir(fullfile(output_path, 'Eigenvector_plots'));
    end
    if ~exist(fullfile(output_path, 'AL_plots'), 'dir')
        mkdir(fullfile(output_path, 'AL_plots'));
    end
    if ~exist(fullfile(output_path, 'ClusterPt_plots'), 'dir')
        mkdir(fullfile(output_path, 'ClusterPt_plots'));
    end

    % Generate plots
    cm_path = fullfile(output_path, 'CM_plots', sprintf('CM_plot_%d.png', config_idx));
    fig1 = plot_cm_offset(output_containers, input_size, config_idx, cm_path);
    close(fig1);

    ev_path = fullfile(output_path, 'Eigenvector_plots', sprintf('Eigenvector_plot_%d.png', config_idx));
    fig2 = plot_eigenvector_offset(output_containers, input_size, config_idx, ev_path);
    close(fig2);

    al_path = fullfile(output_path, 'AL_plots', sprintf('AL_plot_%d.png', config_idx));
    fig3 = plot_anatomical_landmarks(output_containers, input_size, config_idx, al_path);
    close(fig3);

    cp_path = fullfile(output_path, 'ClusterPt_plots', sprintf('ClusterPt_plot_%d.png', config_idx));
    fig4 = plot_cluster_points(output_containers, input_size, config_idx, cp_path);
    close(fig4);
end

function fig = create_summary_figure(final_output_mean, input_size, save_path)
    % CREATE_SUMMARY_FIGURE Create summary figure with mean statistics
    %
    % Creates a comprehensive figure showing mean performance metrics
    % across all configurations.

    fraction_gait = input_size.fraction_completed_gait_cycle;
    TOE_OFF = 0.6;

    color_info = assign_color_info(false);

    fig = figure('Units', 'normalized', 'OuterPosition', [0, 0, 1, 1]);

    % CM offset subplot
    subplot(2, 2, 1);
    xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    plot(fraction_gait, final_output_mean.CM_offset.SVDLS, ...
         'Color', color_info.SVDLS, 'LineWidth', 1.5, 'DisplayName', 'SVDLS');
    plot(fraction_gait, final_output_mean.CM_offset.PTUR, ...
         'Color', color_info.PTUR, 'LineWidth', 1.5, 'DisplayName', 'PTUR');
    grid on;
    legend('Location', 'best');
    xlabel('Fraction of Gait Cycle');
    ylabel('Mean CM Offset (cm)');
    title('Center of Mass Error');

    % Mean eigenvector error
    subplot(2, 2, 2);
    xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    mean_ev_svdls = mean(final_output_mean.eigenvector_offset.SVDLS, 2);
    mean_ev_ptur = mean(final_output_mean.eigenvector_offset.PTUR, 2);
    plot(fraction_gait, mean_ev_svdls, ...
         'Color', color_info.SVDLS, 'LineWidth', 1.5, 'DisplayName', 'SVDLS');
    plot(fraction_gait, mean_ev_ptur, ...
         'Color', color_info.PTUR, 'LineWidth', 1.5, 'DisplayName', 'PTUR');
    grid on;
    legend('Location', 'best');
    xlabel('Fraction of Gait Cycle');
    ylabel('Mean Angular Error (deg)');
    title('Mean Eigenvector Error');

    % Mean cluster point error
    subplot(2, 2, 3);
    xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    mean_cp_svdls = mean(final_output_mean.cluster_offset.SVDLS, 2);
    mean_cp_ptur = mean(final_output_mean.cluster_offset.PTUR, 2);
    plot(fraction_gait, mean_cp_svdls, ...
         'Color', color_info.SVDLS, 'LineWidth', 1.5, 'DisplayName', 'SVDLS');
    plot(fraction_gait, mean_cp_ptur, ...
         'Color', color_info.PTUR, 'LineWidth', 1.5, 'DisplayName', 'PTUR');
    grid on;
    legend('Location', 'best');
    xlabel('Fraction of Gait Cycle');
    ylabel('Mean Error (cm)');
    title('Mean Cluster Point Error');

    % Mean AL error
    subplot(2, 2, 4);
    xline(TOE_OFF, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    mean_al_svdls = mean(final_output_mean.AL_offset.SVDLS, 2);
    mean_al_ptur = mean(final_output_mean.AL_offset.PTUR, 2);
    plot(fraction_gait, mean_al_svdls, ...
         'Color', color_info.SVDLS, 'LineWidth', 1.5, 'DisplayName', 'SVDLS');
    plot(fraction_gait, mean_al_ptur, ...
         'Color', color_info.PTUR, 'LineWidth', 1.5, 'DisplayName', 'PTUR');
    grid on;
    legend('Location', 'best');
    xlabel('Fraction of Gait Cycle');
    ylabel('Mean Error (cm)');
    title('Mean Anatomical Landmark Error');

    sgtitle('PCT Algorithm Performance Summary', 'FontSize', 16);

    if nargin >= 3 && ~isempty(save_path)
        saveas(fig, save_path);
    end
end
