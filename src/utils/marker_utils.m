%% MARKER_UTILS - Marker distribution and noise utilities
%
% This module provides functions for generating, validating, and
% manipulating marker distributions on body segment surfaces.
%
% Functions:
%   generate_random_distribution_constrained - Generate valid marker config
%   check_intermarker_distance              - Validate marker spacing
%   generate_noisy_distribution_STA         - Apply STA noise to markers
%   eval_marker                             - Assign marker to spatial box
%   eval_marker_four                        - Assign marker (4-marker config)
%   apply_current_pose                      - Transform markers to pose
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = marker_utils()
    funcs.generate_random_distribution_constrained = @generate_random_distribution_constrained;
    funcs.check_intermarker_distance = @check_intermarker_distance;
    funcs.generate_noisy_distribution_STA = @generate_noisy_distribution_STA;
    funcs.eval_marker = @eval_marker;
    funcs.eval_marker_four = @eval_marker_four;
    funcs.apply_current_pose = @apply_current_pose;
end

function mass_distribution = generate_random_distribution_constrained(cylinder_data, N)
    % GENERATE_RANDOM_DISTRIBUTION_CONSTRAINED Generate constrained marker distribution
    %
    % Generates a random marker configuration on the cylinder surface
    % with constraints:
    %   - Markers placed only on anterior and lateral aspects
    %   - Markers distributed across height zones
    %   - Anatomically realistic placement
    %
    % Inputs:
    %   cylinder_data - Cylinder surface data
    %   N             - Number of markers (must be even)
    %
    % Output:
    %   mass_distribution - Struct array with marker positions

    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;
    np = cylinder_data.np;
    n = cylinder_data.n;

    % Compute angular coordinates
    theta = rad2deg(atan(Y ./ X));
    first_jump_idx = find(Y(1, :) == max(Y(1, :)));
    second_jump_idx = find(Y(1, :) == min(Y(1, :)));
    theta(:, first_jump_idx:second_jump_idx) = theta(:, first_jump_idx:second_jump_idx) + 180;
    theta(:, second_jump_idx:end) = theta(:, second_jump_idx:end) + 360;

    % Constrained region (anterior + lateral aspects only)
    np_constrained = floor(np / 3) + 1;
    height_division_size = floor(n / (N / 2));

    % Height boundary parameters
    dZ = 0.1;
    deltaZ_max = 6;
    deltaZ_random = 0.5;
    Nmax = deltaZ_max / dZ;
    Nrandom = deltaZ_random / dZ;
    Nsubtract = Nmax + Nrandom + 1;

    % Generate random plane indices for anterior and lateral markers
    plane_idx_anterior = randi([1, floor(np_constrained / 2)], 1, N / 2);
    plane_idx_lateral = randi([floor(np_constrained / 2) + 1, np_constrained], 1, N / 2);

    % Height zone boundaries
    height_lower_idx = Nsubtract;
    height_upper_idx = n - Nsubtract;
    height_idx = nan(N / 2 + 1, 1);
    height_idx(1) = height_lower_idx;
    for k = 1:(N / 2 - 1)
        height_idx(k + 1) = k * height_division_size;
    end
    height_idx(end) = height_upper_idx;

    % Generate marker distribution
    mass_distribution = struct;
    for k = 1:N
        if k <= N / 2
            % Anterior markers
            current_plane_idx = plane_idx_anterior(k);
            current_height_idx = randi([height_idx(k), height_idx(k + 1)], 1, 1);
        else
            % Lateral markers
            k_shifted = k - N / 2;
            current_plane_idx = plane_idx_lateral(k_shifted);
            current_height_idx = randi([height_idx(k_shifted), height_idx(k_shifted + 1)], 1, 1);
        end

        mass_distribution(k).mass = 1;
        mass_distribution(k).plane_idx = current_plane_idx;
        mass_distribution(k).height_idx = current_height_idx;
        mass_distribution(k).X = X(current_height_idx, current_plane_idx);
        mass_distribution(k).Y = Y(current_height_idx, current_plane_idx);
        mass_distribution(k).Z = Z(current_height_idx, current_plane_idx);
        mass_distribution(k).r = sqrt(mass_distribution(k).X^2 + mass_distribution(k).Y^2);
        mass_distribution(k).theta = theta(current_height_idx, current_plane_idx);
        mass_distribution(k).vector = [mass_distribution(k).X; ...
                                        mass_distribution(k).Y; ...
                                        mass_distribution(k).Z];
    end
end

function is_valid = check_intermarker_distance(mass_distribution, N)
    % CHECK_INTERMARKER_DISTANCE Validate inter-marker spacing
    %
    % Checks that markers are sufficiently separated to avoid
    % inter-marker error propagation. Validates both height
    % separation and angular separation.
    %
    % Inputs:
    %   mass_distribution - Marker positions
    %   N                 - Number of markers
    %
    % Output:
    %   is_valid - Boolean indicating valid configuration

    % Threshold parameters
    HEIGHT_DIFF_THRESHOLD = 15;   % cm
    THETA_DIFF_THRESHOLD = 60;    % degrees

    is_valid = false;

    % Extract heights and angles
    height_list_anterior = nan(N / 2, 1);
    height_list_lateral = nan(N / 2, 1);
    theta_diff_list = nan(N / 2, 1);

    for k = 1:N / 2
        k_shifted = k + N / 2;
        height_list_anterior(k) = mass_distribution(k).Z;
        height_list_lateral(k) = mass_distribution(k_shifted).Z;
        theta_diff_list(k) = mass_distribution(k_shifted).theta - mass_distribution(k).theta;
    end

    % Check minimum separations
    height_anterior_diff_min = min(diff(height_list_anterior));
    height_lateral_diff_min = min(diff(height_list_lateral));
    theta_diff_min = min(theta_diff_list);

    % Validate constraints
    height_ok = (height_anterior_diff_min > HEIGHT_DIFF_THRESHOLD) && ...
                (height_lateral_diff_min > HEIGHT_DIFF_THRESHOLD);
    angle_ok = (theta_diff_min > THETA_DIFF_THRESHOLD);

    if height_ok && angle_ok
        is_valid = true;
    end
end

function noisy_distribution = generate_noisy_distribution_STA(mass_distribution_ref, frame_idx, box_markers_gait, N)
    % GENERATE_NOISY_DISTRIBUTION_STA Apply STA noise to markers
    %
    % Applies Soft Tissue Artifact noise from real gait data to
    % the reference marker positions.
    %
    % Inputs:
    %   mass_distribution_ref - Reference marker positions
    %   frame_idx             - Current frame index
    %   box_markers_gait      - STA noise data per marker region
    %   N                     - Number of markers
    %
    % Output:
    %   noisy_distribution - Markers with STA noise applied

    noisy_distribution = struct;

    for j = 1:N
        % Get STA noise for this marker region at this frame
        STA_data = box_markers_gait(j).vector(frame_idx, :)';

        % Convert from mm to cm and apply coordinate rotation
        noise_mm = rotx(90) * STA_data;
        noise_cm = noise_mm / 10;

        % Apply noise to reference position
        noisy_vector = mass_distribution_ref(j).vector + noise_cm;

        noisy_distribution(j).mass = 1;
        noisy_distribution(j).X = noisy_vector(1);
        noisy_distribution(j).Y = noisy_vector(2);
        noisy_distribution(j).Z = noisy_vector(3);
        noisy_distribution(j).vector = noisy_vector;
    end
end

function box_val = eval_marker(theta, h)
    % EVAL_MARKER Assign marker to spatial box (6-box configuration)
    %
    % Classifies a marker into one of 6 spatial regions based on
    % its angular position and height.
    %
    % Inputs:
    %   theta - Angular position (degrees)
    %   h     - Height (cm)
    %
    % Output:
    %   box_val - Box index (1-6)
    %       1: Anterior Low,  2: Anterior Mid,  3: Anterior High
    %       4: Lateral Low,   5: Lateral Mid,   6: Lateral High

    if theta > 0 && theta < 60
        % Anterior region
        if h < 15
            box_val = 1;
        elseif h >= 15 && h < 30
            box_val = 2;
        else
            box_val = 3;
        end
    else
        % Lateral region
        if h < 15
            box_val = 4;
        elseif h >= 15 && h < 30
            box_val = 5;
        else
            box_val = 6;
        end
    end
end

function box_val = eval_marker_four(theta, h)
    % EVAL_MARKER_FOUR Assign marker to spatial box (4-box configuration)
    %
    % Classifies a marker into one of 4 spatial regions.
    %
    % Inputs:
    %   theta - Angular position (degrees)
    %   h     - Height (cm)
    %
    % Output:
    %   box_val - Box index (1-4)
    %       1: Anterior Low,  2: Anterior High
    %       3: Lateral Low,   4: Lateral High

    if theta > 0 && theta < 60
        % Anterior region
        if h < 22.5
            box_val = 1;
        else
            box_val = 2;
        end
    else
        % Lateral region
        if h < 22.5
            box_val = 3;
        else
            box_val = 4;
        end
    end
end

function mass_distribution_pose = apply_current_pose(mass_distribution_standard, current_pose, N)
    % APPLY_CURRENT_POSE Transform markers to specified pose
    %
    % Transforms marker positions from standard pose to a specified
    % global pose using rotation and translation.
    %
    % Inputs:
    %   mass_distribution_standard - Markers in standard pose
    %   current_pose               - Struct with .T and .R
    %   N                          - Number of markers
    %
    % Output:
    %   mass_distribution_pose - Transformed marker positions

    mass_distribution_pose = struct;

    T = current_pose.T;
    R = current_pose.R;

    for j = 1:N
        vec_standard = mass_distribution_standard(j).vector;
        vec_STA = rotx(-90) * vec_standard;
        vec_pose = T + R * vec_STA;

        mass_distribution_pose(j).mass = 1;
        mass_distribution_pose(j).X = vec_pose(1);
        mass_distribution_pose(j).Y = vec_pose(2);
        mass_distribution_pose(j).Z = vec_pose(3);
        mass_distribution_pose(j).vector = vec_pose;
    end
end
