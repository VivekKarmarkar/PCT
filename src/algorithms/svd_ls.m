%% SVD_LS - Singular Value Decomposition Least Squares Algorithm
%
% This module implements the SVD-LS algorithm for pose estimation
% from marker cluster data. The algorithm finds the optimal rigid
% body transformation between local and global marker positions.
%
% Functions:
%   svdls_processing     - Main SVD-LS processing pipeline
%   pose_estimate_SVD_LS - Compute pose using SVD
%   add_inertial_properties_SVDLS - Package results
%   calculate_SVDLS_offset        - Compute estimation errors
%
% Reference:
%   Soderkvist & Wedin (1993) - Determining the movements of the skeleton
%   using well-configured markers
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = svd_ls()
    funcs.svdls_processing = @svdls_processing;
    funcs.pose_estimate_SVD_LS = @pose_estimate_SVD_LS;
    funcs.add_inertial_properties_SVDLS = @add_inertial_properties_SVDLS;
    funcs.calculate_SVDLS_offset = @calculate_SVDLS_offset;
    funcs.evaluate_PCT_SVDLS = @evaluate_PCT_SVDLS;
end

function algorithm_container = svdls_processing(algorithm_container)
    % SVDLS_PROCESSING Main SVD-LS algorithm processing pipeline
    %
    % Processes marker data through the SVD-LS algorithm to estimate
    % the rigid body pose and compute reconstruction errors.
    %
    % Input/Output:
    %   algorithm_container - Struct containing marker data and results

    cylinder_data = algorithm_container.cylinder_data;
    inertial_properties_reference = algorithm_container.inertial_properties_reference;

    % Estimate pose using SVD-LS
    [T_SVDLS, R_SVDLS] = pose_estimate_SVD_LS(algorithm_container);

    % Package inertial properties
    inertial_properties_SVDLS = add_inertial_properties_SVDLS(T_SVDLS, R_SVDLS);

    % Calculate offset from reference
    inertial_properties_SVDLS_offset = calculate_SVDLS_offset(...
        inertial_properties_reference, inertial_properties_SVDLS);

    algorithm_container.inertial_properties_SVDLS = inertial_properties_SVDLS;
    algorithm_container.inertial_properties_SVDLS_offset = inertial_properties_SVDLS_offset;

    % Compute local offsets for cluster points
    algorithm_container = compute_Final_localOffset_SVDLS(algorithm_container);

    % Compute anatomical landmarks
    reference_pose = struct;
    reference_pose.T = inertial_properties_reference.T;
    reference_pose.R = inertial_properties_reference.R;

    optimized_pose = struct;
    optimized_pose.T = inertial_properties_SVDLS.T;
    optimized_pose.R = inertial_properties_SVDLS.R;

    current_pose = algorithm_container.current_pose;

    [anatomical_landmarks_data, anatomical_frame] = compute_anatomical_landmarks_thigh_data(...
        cylinder_data, reference_pose, optimized_pose, current_pose);

    algorithm_container.anatomical_landmarks_thigh_data = anatomical_landmarks_data;
    algorithm_container.anatomical_frame = anatomical_frame;
end

function [T, R] = pose_estimate_SVD_LS(algorithm_container)
    % POSE_ESTIMATE_SVD_LS Compute pose using Singular Value Decomposition
    %
    % Implements the SVD-based least squares solution for finding
    % the optimal rotation and translation between two point sets.
    %
    % Input:
    %   algorithm_container - Contains local and global marker positions
    %
    % Outputs:
    %   T - Translation vector [3x1]
    %   R - Rotation matrix [3x3]

    N = algorithm_container.N;

    % Extract global and local vectors
    global_vecs = struct;
    local_vecs = struct;

    for c = 1:N
        global_vecs(c).vec = algorithm_container.mass_distribution_noisy(c).vector;
        local_vecs(c).vec = algorithm_container.mass_distribution_reference(c).local_displacement;
    end

    % Compute centroids
    global_centroid = zeros(3, 1);
    local_centroid = zeros(3, 1);

    for c = 1:N
        global_centroid = global_centroid + global_vecs(c).vec;
        local_centroid = local_centroid + local_vecs(c).vec;
    end

    global_centroid = global_centroid / N;
    local_centroid = local_centroid / N;

    % Center the point sets
    global_centered = nan(3, N);
    local_centered = nan(3, N);

    for c = 1:N
        global_centered(:, c) = global_vecs(c).vec - global_centroid;
        local_centered(:, c) = local_vecs(c).vec - local_centroid;
    end

    % Compute cross-covariance matrix
    H = local_centered * global_centered';

    % SVD decomposition
    [U, ~, V] = svd(H);

    % Compute rotation
    R = V * U';

    % Ensure proper rotation (det = +1)
    if det(R) < 0
        V(:, 3) = -V(:, 3);
        R = V * U';
    end

    % Compute translation
    T = global_centroid - R * local_centroid;
end

function inertial_properties = add_inertial_properties_SVDLS(T, R)
    % ADD_INERTIAL_PROPERTIES_SVDLS Package SVD-LS results
    %
    % Packages the estimated translation and rotation into the
    % standard inertial properties format.
    %
    % Inputs:
    %   T - Translation vector [3x1]
    %   R - Rotation matrix [3x3]
    %
    % Output:
    %   inertial_properties - Standardized results structure

    inertial_properties = struct;
    inertial_properties.T = T;
    inertial_properties.R = R;

    % Extract CM from translation
    inertial_properties.CM.vector = T;
    inertial_properties.CM.x = T(1);
    inertial_properties.CM.y = T(2);
    inertial_properties.CM.z = T(3);

    % Extract principal axes from rotation matrix columns
    inertial_properties.PrincipalAxis = struct;
    for j = 1:3
        inertial_properties.PrincipalAxis(j).vector = R(:, j);
    end
end

function offset = calculate_SVDLS_offset(ref_props, svdls_props)
    % CALCULATE_SVDLS_OFFSET Compute SVD-LS estimation errors
    %
    % Calculates the angular and positional errors between
    % reference and SVD-LS estimated poses.
    %
    % Inputs:
    %   ref_props   - Reference inertial properties
    %   svdls_props - SVD-LS estimated properties
    %
    % Output:
    %   offset - Structure with error metrics

    offset = struct;

    % Angular offsets for each principal axis
    for j = 1:3
        ref_axis = ref_props.PrincipalAxis(j).vector;
        est_axis = svdls_props.PrincipalAxis(j).vector;

        cos_angle = dot(ref_axis, est_axis);
        cos_angle = max(min(cos_angle, 1), -1);  % Clamp for numerical stability
        offset.angle(j) = rad2deg(acos(cos_angle));
    end

    % CM offset
    offset.CMVector = ref_props.CM.vector - svdls_props.CM.vector;
end

function performance = evaluate_PCT_SVDLS(svdls_offset, original_offset)
    % EVALUATE_PCT_SVDLS Evaluate SVD-LS algorithm performance
    %
    % Compares SVD-LS estimation errors against original (noisy) errors
    % to assess algorithm effectiveness.
    %
    % Inputs:
    %   svdls_offset    - Errors after SVD-LS processing
    %   original_offset - Original (before processing) errors
    %
    % Output:
    %   performance - Performance metrics structure

    performance = struct;

    % Angular performance
    performance.angle.newList = svdls_offset.angle;
    performance.angle.originalList = original_offset.angle;

    % CM distance performance
    performance.CM_distance.newNorm = norm(svdls_offset.CMVector);
    performance.CM_distance.originalNorm = norm(original_offset.CMVector);
end

function algorithm_container = compute_Final_localOffset_SVDLS(algorithm_container)
    % COMPUTE_FINAL_LOCALOFFSET_SVDLS Compute cluster point reconstruction errors
    %
    % Calculates the global position reconstruction error for each
    % marker using the SVD-LS estimated pose.

    N = algorithm_container.N;
    mass_distribution_noisy = algorithm_container.mass_distribution_noisy;
    mass_distribution_ref = algorithm_container.mass_distribution_reference;
    inertial_properties_SVDLS = algorithm_container.inertial_properties_SVDLS;

    T_est = inertial_properties_SVDLS.T;
    R_est = inertial_properties_SVDLS.R;

    mass_distribution_SVDLS = mass_distribution_noisy;

    for k = 1:N
        % Local position from reference
        local_pos = mass_distribution_ref(k).local_displacement;

        % Reconstruct global position
        global_reconstructed = T_est + R_est * local_pos;

        % Original noisy global position
        global_noisy = mass_distribution_noisy(k).vector;

        % Compute offset
        offset = norm(global_reconstructed - global_noisy);

        mass_distribution_SVDLS(k).vector_global_reconstructed = global_reconstructed;
        mass_distribution_SVDLS(k).global_offset = offset;
    end

    algorithm_container.mass_distribution_SVDLS = mass_distribution_SVDLS;
end

function [anatomical_landmarks, anatomical_frame] = compute_anatomical_landmarks_thigh_data(...
        cylinder_data, reference_pose, optimized_pose, current_pose)
    % COMPUTE_ANATOMICAL_LANDMARKS_THIGH_DATA Compute anatomical landmarks
    %
    % Reconstructs anatomical landmark positions using the estimated pose
    % and computes the anatomical reference frame.

    % Thigh anatomical landmarks in local coordinates
    % LE: Lateral Epicondyle, ME: Medial Epicondyle
    % GT: Greater Trochanter, FH: Femoral Head
    AL_local = struct;
    AL_local(1).name = 'Lateral Epicondyle';
    AL_local(1).position = [5; 0; 0];
    AL_local(2).name = 'Medial Epicondyle';
    AL_local(2).position = [-5; 0; 0];
    AL_local(3).name = 'Greater Trochanter';
    AL_local(3).position = [3; 0; 45];
    AL_local(4).name = 'Femoral Head';
    AL_local(4).position = [-3; 0; 48];

    Num_AL = length(AL_local);
    anatomical_landmarks = struct;

    T_ref = reference_pose.T;
    R_ref = reference_pose.R;
    T_opt = optimized_pose.T;
    R_opt = optimized_pose.R;

    for k = 1:Num_AL
        local_pos = AL_local(k).position;

        % Reference global position
        global_ref = T_ref + R_ref * local_pos;

        % Reconstructed global position
        global_recon = T_opt + R_opt * local_pos;

        % Offset
        offset = norm(global_recon - global_ref);

        anatomical_landmarks(k).name = AL_local(k).name;
        anatomical_landmarks(k).local_position = local_pos;
        anatomical_landmarks(k).vector_global_reference = global_ref;
        anatomical_landmarks(k).vector_global_reconstructed = global_recon;
        anatomical_landmarks(k).global_offset = offset;
    end

    % Compute anatomical frame
    anatomical_frame = struct;
    anatomical_frame.T = T_opt;
    anatomical_frame.R = R_opt;
    anatomical_frame.T_ref = T_ref;
    anatomical_frame.R_ref = R_ref;
end
