%% PERTURBATION_THEORY - Perturbation Theory Based Pose Estimation
%
% This module implements the iterative perturbation theory approach
% for PCT pose estimation. The algorithm iteratively corrects the
% eigenvalues while preserving the eigenvector structure.
%
% Functions:
%   perturbation_theory_processing    - Initial PT processing
%   perturbation_theory_reprocessing  - Reprocessing with mass redistribution
%   iterative_perturbation            - Core iterative algorithm
%   compute_mass_redistribution       - Optimal mass solver
%   compute_scaling_factor_CM         - CM constraint scaling
%
% Author: PCT Research Team
% Version: 2.0 (Refactored)

function funcs = perturbation_theory()
    funcs.perturbation_theory_processing = @perturbation_theory_processing;
    funcs.perturbation_theory_reprocessing = @perturbation_theory_reprocessing;
    funcs.iterative_perturbation = @iterative_perturbation;
    funcs.compute_mass_redistribution = @compute_mass_redistribution;
    funcs.compute_scaling_factor_CM = @compute_scaling_factor_CM;
    funcs.evaluate_PCT_perturbation = @evaluate_PCT_perturbation;
end

function algorithm_container = perturbation_theory_processing(algorithm_container)
    % PERTURBATION_THEORY_PROCESSING Initial perturbation theory processing
    %
    % Performs the first pass of perturbation theory to compute the
    % perturbed inertia tensor eigenvalues.
    %
    % Input/Output:
    %   algorithm_container - Struct with marker data and results

    inertial_props_noisy = algorithm_container.inertial_properties_noisy;
    inertial_props_ref = algorithm_container.inertial_properties_reference;

    % Extract starting eigenvalues and eigenvectors
    lambda_start = nan(3, 1);
    lambda_target = nan(3, 1);
    vectors_start = nan(3, 3);

    for j = 1:3
        lambda_start(j) = inertial_props_noisy.PrincipalMOI(j).value;
        lambda_target(j) = inertial_props_ref.PrincipalMOI(j).value;
        vectors_start(:, j) = inertial_props_noisy.PrincipalAxis(j).vector;
    end

    % Run iterative perturbation
    [vectors_final, lambda_final, inertiaTensor_final, ~] = iterative_perturbation(...
        vectors_start, lambda_start, lambda_target);

    % Store results
    algorithm_container.vectors_final = vectors_final;
    algorithm_container.inertiaTensor_final = inertiaTensor_final;
    algorithm_container.eigenvalues_final = lambda_final;
end

function [vectors_final, lambda_final, I_final, convergence_info] = iterative_perturbation(...
        vectors_start, lambda_start, lambda_target)
    % ITERATIVE_PERTURBATION Core iterative perturbation algorithm
    %
    % Iteratively perturbs the inertia tensor to match target eigenvalues
    % while maintaining eigenvector orthogonality.
    %
    % Inputs:
    %   vectors_start - Initial eigenvector matrix [3x3]
    %   lambda_start  - Initial eigenvalues [3x1]
    %   lambda_target - Target eigenvalues [3x1]
    %
    % Outputs:
    %   vectors_final     - Final eigenvector matrix
    %   lambda_final      - Final eigenvalues
    %   I_final           - Final inertia tensor
    %   convergence_info  - Convergence diagnostics

    % Algorithm parameters
    MAX_ITERATIONS = 1000;
    TOLERANCE = 1e-10;
    RELAXATION = 0.1;

    % Initialize
    vectors = vectors_start;
    lambda = lambda_start;
    delta_lambda = lambda_target - lambda;

    convergence_info = struct;
    convergence_info.iterations = 0;
    convergence_info.converged = false;

    for iter = 1:MAX_ITERATIONS
        % Check convergence
        if max(abs(delta_lambda)) < TOLERANCE
            convergence_info.converged = true;
            convergence_info.iterations = iter;
            break;
        end

        % Compute perturbation step
        delta_lambda_step = RELAXATION * delta_lambda;

        % Update eigenvalues
        lambda = lambda + delta_lambda_step;

        % Reconstruct inertia tensor
        D = diag(lambda);
        I_current = vectors * D * vectors';

        % Ensure symmetry
        I_current = (I_current + I_current') / 2;

        % Re-diagonalize
        [vectors, D_new] = eig(I_current);
        lambda = diag(D_new);

        % Ensure right-handed coordinate system
        if det(vectors) < 0
            vectors(:, 3) = -vectors(:, 3);
        end

        % Update delta
        delta_lambda = lambda_target - lambda;
    end

    convergence_info.iterations = iter;
    vectors_final = vectors;
    lambda_final = lambda;
    I_final = vectors * diag(lambda) * vectors';
end

function algorithm_container = perturbation_theory_reprocessing(algorithm_container)
    % PERTURBATION_THEORY_REPROCESSING Reprocess with mass redistribution
    %
    % Second pass processing that incorporates the computed mass
    % redistribution and CM reflection settings.

    N = algorithm_container.N;
    inertial_props_ref = algorithm_container.inertial_properties_reference;
    inertial_props_noisy = algorithm_container.inertial_properties_noisy_original;
    mass_distribution_noisy = algorithm_container.mass_distribution_noisy;
    mass_distribution_ref = algorithm_container.mass_distribution_reference;
    reflection_setting = algorithm_container.reflection_setting;
    m_star = algorithm_container.mass_redistribution;
    scaling_factor_CM = algorithm_container.scaling_factor_CM;

    vectors_final = algorithm_container.vectors_final;
    inertiaTensor_final = algorithm_container.inertiaTensor_final;
    eigenvalues_final = algorithm_container.eigenvalues_final;

    % Compute CM positions
    cm_vector_noisy = inertial_props_noisy.CM.vector;
    [cm_vector_constrained, cm_vector_reflected] = cm_constrain_and_reflect(...
        m_star, mass_distribution_noisy, scaling_factor_CM, inertial_props_noisy, N);

    % Determine reflection based on setting
    cm_vector_relative = cm_vector_noisy - cm_vector_reflected;
    deltaY_offset = dot(cm_vector_relative, vectors_final(:, 1));
    algorithm_container.deltaY_offset = deltaY_offset;

    cm_reflect_bool = false;
    if strcmp(reflection_setting, 'UR')
        if deltaY_offset > 0
            cm_reflect_bool = true;
        end
    elseif strcmp(reflection_setting, 'LR')
        if deltaY_offset < 0
            cm_reflect_bool = true;
        end
    end

    % Update mass distribution
    mass_distribution_perturbed = mass_distribution_noisy;
    for i = 1:N
        mass_distribution_perturbed(i).mass = m_star(i);
    end

    % Compute perturbed inertial properties
    inertial_props_perturbed = struct;
    inertial_props_perturbed.T = cm_vector_constrained;
    if cm_reflect_bool
        inertial_props_perturbed.T = cm_vector_reflected;
    end

    % Compute rotation from translation
    new_R = compute_rotation_given_translation(...
        inertial_props_perturbed.T, N, mass_distribution_noisy, mass_distribution_ref);

    inertial_props_perturbed.R = new_R;
    inertial_props_perturbed.CM.vector = inertial_props_perturbed.T;

    % Extract principal axes from rotation matrix
    inertial_props_perturbed.PrincipalAxis = struct;
    for j = 1:3
        inertial_props_perturbed.PrincipalAxis(j).vector = new_R(:, j);
    end

    % Calculate offset
    inertial_props_perturbed_offset = calculate_offset_perturbed(...
        inertial_props_ref, inertial_props_perturbed);

    algorithm_container.mass_distribution_perturbed = mass_distribution_perturbed;
    algorithm_container.inertial_properties_perturbed = inertial_props_perturbed;
    algorithm_container.inertial_properties_perturbed_offset = inertial_props_perturbed_offset;

    % Compute local offsets
    algorithm_container = compute_Final_localOffset_perturbed(algorithm_container);

    % Compute anatomical landmarks
    cylinder_data = algorithm_container.cylinder_data;
    current_pose = algorithm_container.current_pose;

    reference_pose = struct;
    reference_pose.T = inertial_props_ref.T;
    reference_pose.R = inertial_props_ref.R;

    optimized_pose = struct;
    optimized_pose.T = inertial_props_perturbed.T;
    optimized_pose.R = inertial_props_perturbed.R;

    [anatomical_landmarks, anatomical_frame] = compute_anatomical_landmarks_data(...
        cylinder_data, reference_pose, optimized_pose, current_pose);

    algorithm_container.anatomical_landmarks_thigh_data = anatomical_landmarks;
    algorithm_container.anatomical_frame = anatomical_frame;
end

function mass_redistribution = compute_mass_redistribution(output_containers, input_size)
    % COMPUTE_MASS_REDISTRIBUTION Solve for optimal marker masses
    %
    % Uses Levenberg-Marquardt optimization to find marker masses
    % that minimize CM estimation error across all frames.
    %
    % Inputs:
    %   output_containers - Output data from initial processing
    %   input_size        - Size parameters
    %
    % Output:
    %   mass_redistribution - Optimal marker masses [Nx1]

    N = input_size.N;

    % Define objective function
    objective = @(m) mass_objective_function(m, output_containers, input_size);

    % Initial guess (uniform masses)
    m0 = ones(N, 1);

    % Optimization options
    options = optimoptions('fsolve', ...
        'Algorithm', 'levenberg-marquardt', ...
        'MaxFunctionEvaluations', 10000, ...
        'FunctionTolerance', 1.0, ...
        'Display', 'off');

    % Solve
    mass_redistribution = fsolve(objective, m0, options);
end

function residual = mass_objective_function(m, output_containers, input_size)
    % MASS_OBJECTIVE_FUNCTION Objective for mass redistribution optimization
    %
    % Computes residual for matching weighted CM to target CM.

    N = input_size.N;
    NumFrames = input_size.NumFrames;

    residual = zeros(N, 1);

    for j = 1:NumFrames
        for i = 1:N
            noisy_pos = reshape(output_containers.noisy_positions.PTNR(j, i, :), [3, 1]);
            residual(i) = residual(i) + m(i) * norm(noisy_pos);
        end
    end

    residual = residual / NumFrames;
end

function scaling_factor = compute_scaling_factor_CM(m, output_containers, input_size)
    % COMPUTE_SCALING_FACTOR_CM Compute CM constraint scaling factor
    %
    % Finds the optimal scaling factor to constrain the weighted CM
    % to lie within an acceptable distance from the centroid.

    N = input_size.N;
    NumFrames = input_size.NumFrames;

    % Compute weighted and unweighted CM for each frame
    vector_cm_PT_list = nan(NumFrames, 3);
    vector_cm_N_list = nan(NumFrames, 3);
    distance_threshold_list = nan(NumFrames, 1);

    for j = 1:NumFrames
        cm_weighted = zeros(3, 1);
        cm_unweighted = zeros(3, 1);
        m_sum = 0;

        for i = 1:N
            noisy_pos = reshape(output_containers.noisy_positions.PTNR(j, i, :), [3, 1]);
            cm_weighted = cm_weighted + m(i) * noisy_pos;
            cm_unweighted = cm_unweighted + noisy_pos;
            m_sum = m_sum + m(i);
        end

        vector_cm_PT_list(j, :) = cm_weighted' / m_sum;
        vector_cm_N_list(j, :) = cm_unweighted' / N;

        % Compute distance threshold based on gait cycle
        frac = input_size.fraction_completed_gait_cycle(j);
        distance_threshold_list(j) = compute_cm_dist_threshold_STA(frac);
    end

    % Find maximum valid scaling factor
    scaling_factor = 1.0;
    while scaling_factor >= 0.1
        cm_constrained = (1 - scaling_factor) * vector_cm_N_list + ...
                         scaling_factor * vector_cm_PT_list;
        cm_relative = cm_constrained - vector_cm_N_list;
        distances = vecnorm(cm_relative, 2, 2);

        if all(distances < distance_threshold_list)
            break;
        end
        scaling_factor = scaling_factor - 0.1;
    end
end

function threshold = compute_cm_dist_threshold_STA(fraction_completed)
    % COMPUTE_CM_DIST_THRESHOLD_STA Compute CM distance threshold
    %
    % Returns the allowable CM offset based on gait cycle phase.

    if fraction_completed <= 0.5
        threshold = 0.7 + 2.2 * fraction_completed;
    else
        threshold = 2.9 - 2.2 * fraction_completed;
    end
end

function [cm_constrained, cm_reflected] = cm_constrain_and_reflect(m, mass_dist, scale, props, N)
    % CM_CONSTRAIN_AND_REFLECT Compute constrained and reflected CM
    %
    % Applies mass redistribution constraint and reflection.

    % Compute weighted CM
    cm_weighted = zeros(3, 1);
    cm_unweighted = zeros(3, 1);
    m_sum = 0;

    for i = 1:N
        pos = mass_dist(i).vector;
        cm_weighted = cm_weighted + m(i) * pos;
        cm_unweighted = cm_unweighted + pos;
        m_sum = m_sum + m(i);
    end

    cm_weighted = cm_weighted / m_sum;
    cm_unweighted = cm_unweighted / N;

    % Apply scaling constraint
    cm_constrained = (1 - scale) * cm_unweighted + scale * cm_weighted;

    % Reflection across centroid plane
    cm_reflected = 2 * cm_unweighted - cm_constrained;
end

function R = compute_rotation_given_translation(T, N, mass_dist_noisy, mass_dist_ref)
    % COMPUTE_ROTATION_GIVEN_TRANSLATION Estimate rotation from known translation
    %
    % Uses quaternion-based optimization to find the rotation matrix
    % that best aligns the marker positions given the translation.

    D = zeros(4, 4);

    for k = 1:N
        pos_global = mass_dist_noisy(k).vector;
        pos_global_translated = pos_global - T;
        pos_local = mass_dist_ref(k).local_displacement;

        A = leftquat(pos_global_translated) - rightquat(pos_local);
        D = D + A' * A;
    end

    % Find eigenvector corresponding to smallest eigenvalue
    [evec, eval] = eig(D);
    [~, idx] = sort(diag(eval));
    q_star = evec(:, idx(1));

    % Convert quaternion to rotation matrix
    R = quat2rotm(q_star');
end

function L = leftquat(v)
    % LEFTQUAT Left quaternion multiplication matrix
    L = [0     -v(1)  -v(2)  -v(3);
         v(1)   0     -v(3)   v(2);
         v(2)   v(3)   0     -v(1);
         v(3)  -v(2)   v(1)   0   ];
end

function R = rightquat(v)
    % RIGHTQUAT Right quaternion multiplication matrix
    R = [0     -v(1)  -v(2)  -v(3);
         v(1)   0      v(3)  -v(2);
         v(2)  -v(3)   0      v(1);
         v(3)   v(2)  -v(1)   0   ];
end

function offset = calculate_offset_perturbed(ref_props, perturbed_props)
    % CALCULATE_OFFSET_PERTURBED Calculate perturbation estimation errors

    offset = struct;

    for j = 1:3
        ref_axis = ref_props.PrincipalAxis(j).vector;
        pert_axis = perturbed_props.PrincipalAxis(j).vector;

        cos_angle = dot(ref_axis, pert_axis);
        cos_angle = max(min(cos_angle, 1), -1);
        offset.angle(j) = rad2deg(acos(cos_angle));
    end

    offset.CMVector = ref_props.CM.vector - perturbed_props.CM.vector;
end

function performance = evaluate_PCT_perturbation(perturbed_offset, original_offset)
    % EVALUATE_PCT_PERTURBATION Evaluate perturbation algorithm performance

    performance = struct;
    performance.angle.newList = perturbed_offset.angle;
    performance.angle.originalList = original_offset.angle;
    performance.CM_distance.newNorm = norm(perturbed_offset.CMVector);
    performance.CM_distance.originalNorm = norm(original_offset.CMVector);
end

function algorithm_container = compute_Final_localOffset_perturbed(algorithm_container)
    % COMPUTE_FINAL_LOCALOFFSET_PERTURBED Compute cluster point errors

    N = algorithm_container.N;
    mass_dist_perturbed = algorithm_container.mass_distribution_perturbed;
    mass_dist_ref = algorithm_container.mass_distribution_reference;
    inertial_props = algorithm_container.inertial_properties_perturbed;

    T_est = inertial_props.T;
    R_est = inertial_props.R;

    for k = 1:N
        local_pos = mass_dist_ref(k).local_displacement;
        global_recon = T_est + R_est * local_pos;
        global_noisy = mass_dist_perturbed(k).vector;

        offset = norm(global_recon - global_noisy);
        mass_dist_perturbed(k).vector_global_reconstructed = global_recon;
        mass_dist_perturbed(k).global_offset = offset;
    end

    algorithm_container.mass_distribution_perturbed = mass_dist_perturbed;
end

function [anatomical_landmarks, anatomical_frame] = compute_anatomical_landmarks_data(...
        cylinder_data, reference_pose, optimized_pose, current_pose)
    % COMPUTE_ANATOMICAL_LANDMARKS_DATA Compute anatomical landmarks

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
        global_ref = T_ref + R_ref * local_pos;
        global_recon = T_opt + R_opt * local_pos;
        offset = norm(global_recon - global_ref);

        anatomical_landmarks(k).name = AL_local(k).name;
        anatomical_landmarks(k).local_position = local_pos;
        anatomical_landmarks(k).vector_global_reference = global_ref;
        anatomical_landmarks(k).vector_global_reconstructed = global_recon;
        anatomical_landmarks(k).global_offset = offset;
    end

    anatomical_frame = struct;
    anatomical_frame.T = T_opt;
    anatomical_frame.R = R_opt;
    anatomical_frame.T_ref = T_ref;
    anatomical_frame.R_ref = R_ref;
end
