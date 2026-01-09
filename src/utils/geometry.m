%% GEOMETRY - Geometric utility functions for PCT
%
% This module provides geometric transformation and coordinate system
% utilities used throughout the PCT algorithm implementation.
%
% Functions:
%   local2global_perturbed  - Transform local coordinates to global frame
%   global2local            - Transform global coordinates to local frame
%   leftquat                - Left quaternion multiplication matrix
%   rightquat               - Right quaternion multiplication matrix
%   eval_theta              - Calculate azimuthal angle from x,z coordinates
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = geometry()
    funcs.local2global_perturbed = @local2global_perturbed;
    funcs.global2local = @global2local;
    funcs.leftquat = @leftquat;
    funcs.rightquat = @rightquat;
    funcs.eval_theta = @eval_theta;
end

function position_global = local2global_perturbed(T, R, position_local)
    % LOCAL2GLOBAL_PERTURBED Transform local coordinates to global frame
    %
    % Applies rotation R and translation T to convert position from
    % local coordinate frame to global coordinate frame.
    %
    % Inputs:
    %   T              - Translation vector [3x1]
    %   R              - Rotation matrix [3x3]
    %   position_local - Position in local frame [3x1]
    %
    % Output:
    %   position_global - Position in global frame [3x1]

    position_global = T + R * position_local;
end

function position_local = global2local(T, R, position_global)
    % GLOBAL2LOCAL Transform global coordinates to local frame
    %
    % Applies inverse rotation and translation to convert position from
    % global coordinate frame to local coordinate frame.
    %
    % Inputs:
    %   T               - Translation vector [3x1]
    %   R               - Rotation matrix [3x3]
    %   position_global - Position in global frame [3x1]
    %
    % Output:
    %   position_local - Position in local frame [3x1]

    position_local = R' * (position_global - T);
end

function L = leftquat(v)
    % LEFTQUAT Create left quaternion multiplication matrix
    %
    % Constructs the 4x4 matrix for left quaternion multiplication
    % used in quaternion-based rotation estimation.
    %
    % Input:
    %   v - 3D vector [3x1]
    %
    % Output:
    %   L - Left quaternion matrix [4x4]

    L = [0     -v(1)  -v(2)  -v(3);
         v(1)   0     -v(3)   v(2);
         v(2)   v(3)   0     -v(1);
         v(3)  -v(2)   v(1)   0   ];
end

function R = rightquat(v)
    % RIGHTQUAT Create right quaternion multiplication matrix
    %
    % Constructs the 4x4 matrix for right quaternion multiplication
    % used in quaternion-based rotation estimation.
    %
    % Input:
    %   v - 3D vector [3x1]
    %
    % Output:
    %   R - Right quaternion matrix [4x4]

    R = [0     -v(1)  -v(2)  -v(3);
         v(1)   0      v(3)  -v(2);
         v(2)  -v(3)   0      v(1);
         v(3)   v(2)  -v(1)   0   ];
end

function theta_val = eval_theta(x, z)
    % EVAL_THETA Calculate azimuthal angle from x,z coordinates
    %
    % Computes the azimuthal angle in degrees for a point given its
    % x and z coordinates. Handles all four quadrants correctly.
    %
    % Inputs:
    %   x - X coordinate
    %   z - Z coordinate
    %
    % Output:
    %   theta_val - Azimuthal angle in degrees [0, 360)

    z_reflected = -z;

    if x > 0 && z_reflected > 0
        theta_val = rad2deg(atan(z_reflected / x));
    elseif x < 0 && z_reflected > 0
        theta_val = 180 + rad2deg(atan(z_reflected / x));
    elseif x < 0 && z_reflected < 0
        theta_val = 180 + rad2deg(atan(z_reflected / x));
    else
        theta_val = 360 + rad2deg(atan(z_reflected / x));
    end
end
