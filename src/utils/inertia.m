%% INERTIA - Inertial properties computation utilities
%
% This module provides functions for computing and manipulating
% inertial properties of marker clusters.
%
% Functions:
%   compute_inertial_properties   - Compute full inertial properties
%   calculate_offset              - Calculate offset between reference and noisy
%   axis_reversal_original        - Handle eigenvector sign ambiguity
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = inertia()
    funcs.compute_inertial_properties = @compute_inertial_properties;
    funcs.calculate_offset = @calculate_offset;
    funcs.axis_reversal_original = @axis_reversal_original;
end

function [inertial_properties, mass_distribution] = compute_inertial_properties(mass_distribution, N)
    % COMPUTE_INERTIAL_PROPERTIES Compute inertial properties of marker cluster
    %
    % Calculates center of mass, moment of inertia tensor, and principal
    % axes for a distribution of markers with associated masses.
    %
    % Inputs:
    %   mass_distribution - Struct array with marker positions and masses
    %   N                 - Number of markers
    %
    % Outputs:
    %   inertial_properties - Struct containing:
    %       .CM              - Center of mass coordinates
    %       .MOI             - Moment of inertia tensor [3x3]
    %       .PrincipalAxis   - Principal axis directions
    %       .PrincipalMOI    - Principal moments of inertia
    %       .R               - Rotation matrix (eigenvectors)
    %       .T               - Translation vector (CM position)
    %   mass_distribution   - Updated with relative positions

    inertial_properties = struct;

    % Compute center of mass
    M_sum = 0;
    CM_X_sum = 0;
    CM_Y_sum = 0;
    CM_Z_sum = 0;

    for k = 1:N
        m_k = mass_distribution(k).mass;
        M_sum = M_sum + m_k;
        CM_X_sum = CM_X_sum + m_k * mass_distribution(k).X;
        CM_Y_sum = CM_Y_sum + m_k * mass_distribution(k).Y;
        CM_Z_sum = CM_Z_sum + m_k * mass_distribution(k).Z;
    end

    CM_X = CM_X_sum / M_sum;
    CM_Y = CM_Y_sum / M_sum;
    CM_Z = CM_Z_sum / M_sum;

    % Compute moment of inertia tensor
    I_XX = 0; I_YY = 0; I_ZZ = 0;
    I_XY = 0; I_XZ = 0; I_YZ = 0;

    for k = 1:N
        % Compute relative positions
        mass_distribution(k).X_rel = mass_distribution(k).X - CM_X;
        mass_distribution(k).Y_rel = mass_distribution(k).Y - CM_Y;
        mass_distribution(k).Z_rel = mass_distribution(k).Z - CM_Z;

        m_k = mass_distribution(k).mass;
        x_rel = mass_distribution(k).X_rel;
        y_rel = mass_distribution(k).Y_rel;
        z_rel = mass_distribution(k).Z_rel;

        % Diagonal elements
        I_XX = I_XX + m_k * (y_rel^2 + z_rel^2);
        I_YY = I_YY + m_k * (x_rel^2 + z_rel^2);
        I_ZZ = I_ZZ + m_k * (x_rel^2 + y_rel^2);

        % Off-diagonal elements (products of inertia)
        I_XY = I_XY - m_k * x_rel * y_rel;
        I_XZ = I_XZ - m_k * x_rel * z_rel;
        I_YZ = I_YZ - m_k * y_rel * z_rel;
    end

    % Construct inertia tensor
    I = [I_XX, I_XY, I_XZ;
         I_XY, I_YY, I_YZ;
         I_XZ, I_YZ, I_ZZ];

    % Compute eigendecomposition
    tol_orthogonal = 1e-10;
    try
        [V, D] = eig(I);
    catch
        V = nan(3, 3);
        D = diag(nan(3, 1));
    end

    diag_MOI = diag(D);

    % Ensure right-handed coordinate system
    if dot(cross(V(:,1), V(:,2)), V(:,3)) + 1 < tol_orthogonal
        V(:,3) = -V(:,3);
    end

    % Populate output structure
    inertial_properties.CM.x = CM_X;
    inertial_properties.CM.y = CM_Y;
    inertial_properties.CM.z = CM_Z;
    inertial_properties.CM.vector = [CM_X; CM_Y; CM_Z];
    inertial_properties.MOI = I;

    inertial_properties.PrincipalAxis = struct;
    inertial_properties.PrincipalMOI = struct;
    for j = 1:3
        inertial_properties.PrincipalAxis(j).vector = V(:, j);
        inertial_properties.PrincipalMOI(j).value = diag_MOI(j);
    end

    inertial_properties.PrincipalMOINorm = norm(diag_MOI);
    inertial_properties.R = V;
    inertial_properties.T = inertial_properties.CM.vector;
end

function offset = calculate_offset(ref_props, noisy_props)
    % CALCULATE_OFFSET Calculate offset between reference and noisy properties
    %
    % Computes angular differences between principal axes and
    % differences in principal moments of inertia.
    %
    % Inputs:
    %   ref_props   - Reference inertial properties
    %   noisy_props - Noisy inertial properties
    %
    % Output:
    %   offset - Struct containing angular and MOI differences

    offset = struct;

    for j = 1:3
        ref_axis = ref_props.PrincipalAxis(j).vector;
        noisy_axis = noisy_props.PrincipalAxis(j).vector;

        % Angular offset in degrees
        offset.angle(j) = rad2deg(acos(dot(ref_axis, noisy_axis)));

        % MOI difference
        offset.PrincipalMOI(j) = ref_props.PrincipalMOI(j).value - ...
                                  noisy_props.PrincipalMOI(j).value;
    end

    offset.PrincipalMOINorm = ref_props.PrincipalMOINorm - noisy_props.PrincipalMOINorm;
    offset.CMVector = ref_props.CM.vector - noisy_props.CM.vector;
end

function props_reversed = axis_reversal_original(props, offset)
    % AXIS_REVERSAL_ORIGINAL Handle eigenvector sign ambiguity
    %
    % Eigenvectors have sign ambiguity (v and -v are both valid).
    % This function reverses axes where the angular offset exceeds 90
    % degrees to ensure consistent orientation.
    %
    % Inputs:
    %   props  - Inertial properties to potentially reverse
    %   offset - Calculated offset with angle information
    %
    % Output:
    %   props_reversed - Properties with corrected axis orientations

    orthogonal_tol = 1e-10;
    props_reversed = props;
    reversal_flag = false;

    for m = 1:3
        if offset.angle(m) > 90
            reversal_flag = true;
            props_reversed.PrincipalAxis(m).vector = -props.PrincipalAxis(m).vector;
            props_reversed.R(:, m) = -props.R(:, m);
        end
    end

    if reversal_flag
        if abs(det(props_reversed.R) - 1) < orthogonal_tol
            disp('Axis reversal successful');
        end
    end
end
