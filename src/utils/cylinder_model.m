%% CYLINDER_MODEL - Cylinder geometry generation utilities
%
% This module provides functions for generating cylinder surface models
% representing body segments (thigh, shank) as truncated cones.
%
% Functions:
%   generate_cylinder_data       - Generate standard cylinder model
%   generate_cylinder_data_thigh - Generate thigh-specific cylinder
%   generate_cylinder_data_shank - Generate shank-specific cylinder
%   generate_cylinder_data_pose  - Transform cylinder to pose
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = cylinder_model()
    funcs.generate_cylinder_data = @generate_cylinder_data;
    funcs.generate_cylinder_data_thigh = @generate_cylinder_data_thigh;
    funcs.generate_cylinder_data_shank = @generate_cylinder_data_shank;
    funcs.generate_cylinder_data_pose = @generate_cylinder_data_pose;
end

function cylinder_data = generate_cylinder_data()
    % GENERATE_CYLINDER_DATA Generate standard truncated cone cylinder model
    %
    % Creates a parametric cylinder model with linearly varying radius,
    % representing a generic body segment surface.
    %
    % Output:
    %   cylinder_data - Struct containing:
    %       .X, .Y, .Z - Surface coordinate matrices [n x np]
    %       .np        - Number of circumferential points
    %       .n         - Number of height points

    % Height discretization
    t = 0:0.1:50;
    n = length(t);

    % Circumferential discretization
    np = 100;

    % Linearly varying radius (truncated cone)
    r = 1 + 0.05 * t;

    % Generate cylinder surface
    [X, Y, Z] = cylinder(r, np - 1);
    Z_scaled = 50 * Z;

    % Package output
    cylinder_data = struct;
    cylinder_data.X = X;
    cylinder_data.Y = Y;
    cylinder_data.Z = Z_scaled;
    cylinder_data.np = np;
    cylinder_data.n = n;
end

function cylinder_data = generate_cylinder_data_thigh()
    % GENERATE_CYLINDER_DATA_THIGH Generate thigh-specific cylinder model
    %
    % Creates a cylinder model with dimensions appropriate for the
    % thigh segment based on anthropometric data.
    %
    % Output:
    %   cylinder_data - Struct with thigh cylinder surface data

    % Thigh-specific parameters
    height_max = 45;  % cm
    t = 0:0.1:height_max;
    n = length(t);
    np = 100;

    % Thigh radius profile (tapers from hip to knee)
    r_proximal = 8;   % cm (at hip)
    r_distal = 5;     % cm (at knee)
    r = r_proximal + (r_distal - r_proximal) * t / height_max;

    [X, Y, Z] = cylinder(r, np - 1);
    Z_scaled = height_max * Z;

    cylinder_data = struct;
    cylinder_data.X = X;
    cylinder_data.Y = Y;
    cylinder_data.Z = Z_scaled;
    cylinder_data.np = np;
    cylinder_data.n = n;
    cylinder_data.segment = 'thigh';
end

function cylinder_data = generate_cylinder_data_shank()
    % GENERATE_CYLINDER_DATA_SHANK Generate shank-specific cylinder model
    %
    % Creates a cylinder model with dimensions appropriate for the
    % shank (lower leg) segment based on anthropometric data.
    %
    % Output:
    %   cylinder_data - Struct with shank cylinder surface data

    % Shank-specific parameters
    height_max = 40;  % cm
    t = 0:0.1:height_max;
    n = length(t);
    np = 100;

    % Shank radius profile (tapers from knee to ankle)
    r_proximal = 5;   % cm (at knee)
    r_distal = 3;     % cm (at ankle)
    r = r_proximal + (r_distal - r_proximal) * t / height_max;

    [X, Y, Z] = cylinder(r, np - 1);
    Z_scaled = height_max * Z;

    cylinder_data = struct;
    cylinder_data.X = X;
    cylinder_data.Y = Y;
    cylinder_data.Z = Z_scaled;
    cylinder_data.np = np;
    cylinder_data.n = n;
    cylinder_data.segment = 'shank';
end

function cylinder_data_pose = generate_cylinder_data_pose(cylinder_data, current_pose)
    % GENERATE_CYLINDER_DATA_POSE Transform cylinder to specified pose
    %
    % Transforms the cylinder surface from the standard pose to a
    % specified global pose using rotation and translation.
    %
    % Inputs:
    %   cylinder_data - Standard cylinder data structure
    %   current_pose  - Struct with .T (translation) and .R (rotation)
    %
    % Output:
    %   cylinder_data_pose - Transformed cylinder surface data

    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;
    n = cylinder_data.n;
    np = cylinder_data.np;

    T = current_pose.T;
    R = current_pose.R;

    % Initialize transformed coordinates
    X_pose = nan(n, np);
    Y_pose = nan(n, np);
    Z_pose = nan(n, np);

    % Transform each surface point
    for k = 1:n
        for j = 1:np
            % Standard pose vector
            vec_standard = [X(k, j); Y(k, j); Z(k, j)];

            % Apply STA coordinate system rotation
            vec_STA = rotx(-90) * vec_standard;

            % Apply pose transformation
            vec_pose = T + R * vec_STA;

            X_pose(k, j) = vec_pose(1);
            Y_pose(k, j) = vec_pose(2);
            Z_pose(k, j) = vec_pose(3);
        end
    end

    cylinder_data_pose = struct;
    cylinder_data_pose.X = X_pose;
    cylinder_data_pose.Y = Y_pose;
    cylinder_data_pose.Z = Z_pose;
    cylinder_data_pose.np = np;
    cylinder_data_pose.n = n;
end
