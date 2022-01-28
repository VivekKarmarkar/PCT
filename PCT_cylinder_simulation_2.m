N = 6;
unconstrained_optimization_required = false;
visualize_flag = false;
visualize_good_results = false;
mass_distribution_type = 'otherwise';
initial_axes_reversal = false;
initial_axes_perturbed_reversal = true;
optimized_axes_reversal = true;

folderName_settings = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\MarkerConfigurationConstrainedSettings_" + string(N);
NumTrials = length(dir(folderName_settings)) - 2;
Num_AL_Thigh = 4;
bad_list = [];

xyz_points = nan(N,3);
plane_bool = false(NumTrials,1);

MOInorm_tol = 0.1;
cm_dist_threshold = 4;
plane_distThreshold = 1;

cylinder_data = generate_cylinder_data;

trial_idx_list = 1:NumTrials;
fixed_noise_array = nan(NumTrials, N);

angleOne_old_offset = nan(NumTrials, 1);
angleTwo_old_offset = nan(NumTrials, 1);
angleThree_old_offset = nan(NumTrials, 1);

angleOne_new_offset = nan(NumTrials, 1);
angleTwo_new_offset = nan(NumTrials, 1);
angleThree_new_offset = nan(NumTrials, 1);

angleOne_new_offsetOne = nan(NumTrials, 1);
angleTwo_new_offsetOne = nan(NumTrials, 1);
angleThree_new_offsetOne = nan(NumTrials, 1);

angleOne_new_offsetDelta = nan(NumTrials, 1);
angleTwo_new_offsetDelta = nan(NumTrials, 1);
angleThree_new_offsetDelta = nan(NumTrials, 1);

angleOne_new_offset2D = nan(NumTrials, 1);
angleTwo_new_offset2D = nan(NumTrials, 1);
angleThree_new_offset2D = nan(NumTrials, 1);

angleOne_new_offset2D_constrained = nan(NumTrials, 1);
angleTwo_new_offset2D_constrained = nan(NumTrials, 1);
angleThree_new_offset2D_constrained = nan(NumTrials, 1);

angleOne_new_offset2D_constrained_new_1 = nan(NumTrials, 1);
angleTwo_new_offset2D_constrained_new_1 = nan(NumTrials, 1);
angleThree_new_offset2D_constrained_new_1 = nan(NumTrials, 1);

angleOne_new_offset2D_constrained_new_2 = nan(NumTrials, 1);
angleTwo_new_offset2D_constrained_new_2 = nan(NumTrials, 1);
angleThree_new_offset2D_constrained_new_2 = nan(NumTrials, 1);

angleOne_new_offset2D_constrained_new_3 = nan(NumTrials, 1);
angleTwo_new_offset2D_constrained_new_3 = nan(NumTrials, 1);
angleThree_new_offset2D_constrained_new_3 = nan(NumTrials, 1);

angleOne_new_offset3D = nan(NumTrials, 1);
angleTwo_new_offset3D = nan(NumTrials, 1);
angleThree_new_offset3D = nan(NumTrials, 1);

angleOne_new_offsetPerturbed = nan(NumTrials, 1);
angleTwo_new_offsetPerturbed = nan(NumTrials, 1);
angleThree_new_offsetPerturbed = nan(NumTrials, 1);

angleOne_new_offsetReflected = nan(NumTrials, 1);
angleTwo_new_offsetReflected = nan(NumTrials, 1);
angleThree_new_offsetReflected = nan(NumTrials, 1);

angleOne_new_offsetHybrid = nan(NumTrials, 1);
angleTwo_new_offsetHybrid = nan(NumTrials, 1);
angleThree_new_offsetHybrid = nan(NumTrials, 1);

CMnorm_old_offset = nan(NumTrials, 1);
CMnorm_new_offset = nan(NumTrials, 1);
CMnorm_new_offsetOne = nan(NumTrials, 1);
CMnorm_new_offsetDelta = nan(NumTrials, 1);
CMnorm_new_offset2D = nan(NumTrials, 1);
CMnorm_new_offset2D_constrained = nan(NumTrials, 1);
CMnorm_new_offset2D_constrained_new_1 = nan(NumTrials, 1);
CMnorm_new_offset2D_constrained_new_2 = nan(NumTrials, 1);
CMnorm_new_offset2D_constrained_new_3 = nan(NumTrials, 1);
CMnorm_new_offset3D = nan(NumTrials, 1);
CMnorm_new_offsetPerturbed = nan(NumTrials, 1);
CMnorm_new_offsetReflected = nan(NumTrials, 1);
CMnorm_new_offsetHybrid = nan(NumTrials, 1);

CMnorm_new_offset_reflected = nan(NumTrials, 1);

eigenvalueChange_old_offset = nan(NumTrials, 1);
eigenvalueChange_new_offset = nan(NumTrials, 1);
eigenvalueChange_new_offsetOne = nan(NumTrials, 1);
eigenvalueChange_new_offsetDelta = nan(NumTrials, 1);
eigenvalueChange_new_offset2D = nan(NumTrials, 1);
eigenvalueChange_new_offset2D_constrained = nan(NumTrials, 1);
eigenvalueChange_new_offset2D_constrained_new = nan(NumTrials, 1);
eigenvalueChange_new_offset2D_constrained_new_1 = nan(NumTrials, 1);
eigenvalueChange_new_offset2D_constrained_new_2 = nan(NumTrials, 1);
eigenvalueChange_new_offset2D_constrained_new_3 = nan(NumTrials, 1);
eigenvalueChange_new_offset3D = nan(NumTrials, 1);
eigenvalueChange_new_offsetPerturbed = nan(NumTrials, 1);
eigenvalueChange_new_offsetReflected = nan(NumTrials, 1);
eigenvalueChange_new_offsetHybrid = nan(NumTrials, 1);

global_offset_data = nan(NumTrials, N, 1);
global_offset_data_noisy = nan(NumTrials, N, 1);
global_offset_dataOne = nan(NumTrials, N, 1);
global_offset_dataDelta = nan(NumTrials, N, 1);
global_offset_data2D = nan(NumTrials, N, 1);
global_offset_data2D_constrained = nan(NumTrials, N, 1);
global_offset_data2D_constrained_new_1= nan(NumTrials, N, 1);
global_offset_data2D_constrained_new_2= nan(NumTrials, N, 1);
global_offset_data2D_constrained_new_3= nan(NumTrials, N, 1);
global_offset_data3D = nan(NumTrials, N, 1);
global_offset_data_perturbed = nan(NumTrials, N, 1);
global_offset_data_reflected = nan(NumTrials, N, 1);
global_offset_data_hybrid = nan(NumTrials, N, 1);

global_offset_data_AL = nan(NumTrials, Num_AL_Thigh, 1);
global_offset_data_perturbed_AL = nan(NumTrials, Num_AL_Thigh, 1);
global_offset_data_reflected_AL = nan(NumTrials, Num_AL_Thigh, 1);
global_offset_data_hybrid_AL = nan(NumTrials, Num_AL_Thigh, 1);

eps_iter_list = 1:0.01:5;
eps_MOI_offset = nan(NumTrials, length(eps_iter_list));
eps_First_MOI_offset = nan(NumTrials, length(eps_iter_list));

BestAlgorithm = strings(NumTrials,1);
BestConfiguration = struct;
BestConfiguration.minMOI = nan(NumTrials,1);
BestConfiguration.normMOI = nan(NumTrials,1);
BestConfiguration.firstRatioMOI = nan(NumTrials,1);
BestConfiguration.secondRatioMOI = nan(NumTrials,1);
BestConfiguration.thirdRatioMOI = nan(NumTrials,1);
BestConfiguration.normMOIdiff = nan(NumTrials,1);
BestConfiguration.firstMOIdiff = nan(NumTrials,1);
BestConfiguration.secondMOIdiff = nan(NumTrials,1);
BestConfiguration.thirdMOIdiff = nan(NumTrials,1);

for m=1:NumTrials
    fileName_settings = folderName_settings + "\FixedMarkerDistribution_" + string(m);
    load(fileName_settings);

    mass_distribution_noisy = generate_noisy_ditribution(mass_distribution_reference, N, cylinder_data);
    [inertial_properties_reference, mass_distribution_reference] = compute_inertial_properties(mass_distribution_reference, N);
    [inertial_properties_noisy, mass_distribution_noisy] = compute_inertial_properties(mass_distribution_noisy, N);
    sanity_check_bool_reference = sanity_check_axes(mass_distribution_reference, inertial_properties_reference, N);
    sanity_check_bool_noisy = sanity_check_axes(mass_distribution_noisy, inertial_properties_noisy, N);

    inertial_properties_offset = calculate_offset(inertial_properties_reference, inertial_properties_noisy);
    if initial_axes_reversal
        inertial_properties_noisy_reversed = axis_reversal_original(inertial_properties_noisy, inertial_properties_offset);
        inertial_properties_noisy = inertial_properties_noisy_reversed;
        inertial_properties_offset = calculate_offset(inertial_properties_reference, inertial_properties_noisy);
    end
    
    inertial_properties_noisy_copy = inertial_properties_noisy;
    inertial_properties_offset_copy = calculate_offset(inertial_properties_reference, inertial_properties_noisy_copy);
    disp("Original offset");
    disp(inertial_properties_offset_copy.angle);
    if initial_axes_perturbed_reversal
        inertial_properties_noisy_reversed = axis_reversal_original(inertial_properties_noisy_copy, inertial_properties_offset_copy);
        inertial_properties_noisy_copy = inertial_properties_noisy_reversed;
        inertial_properties_offset_copy = calculate_offset(inertial_properties_reference, inertial_properties_noisy_copy);
    end
    disp("New offset");
    disp(inertial_properties_offset_copy.angle);
    
    PCT_tol = 0.1;
    if abs(inertial_properties_offset.PrincipalMOINorm) >= PCT_tol
        unconstrained_optimization_required = true;
    end

    visualize_input = struct;
    visualize_input.TrialNumber = m;
    visualize_input.visualize_flag = visualize_flag;
    visualize_input.mass_distribution_reference = mass_distribution_reference;
    visualize_input.mass_distribution_noisy = mass_distribution_noisy;
    visualize_input.inertial_properties_reference = inertial_properties_reference;
    visualize_input.inertial_properties_noisy = inertial_properties_noisy;
    visualize_input.inertial_properties_offset = inertial_properties_offset;
    visualize_input.N = N;
    visualize_input.cylinder_data = cylinder_data;
    visualize_input.optimized_axis_reversal_bool = optimized_axes_reversal;
    
    [visualize_input, local_offset, local_angular_offset, local_offset_vector] = compute_localOffset(visualize_input);
    visualize_input.local_offset = local_offset;
    visualize_input.local_angular_offset = local_angular_offset;
    visualize_input.local_offset_vector = local_offset_vector;
    
    visualize_inputOne = visualize_input;
    visualize_inputDelta = visualize_input;
    visualize_input2D = visualize_input;
    visualize_input2D_constrained = visualize_input;
    visualize_input2D_constrained_new_1 = visualize_input;
    visualize_input2D_constrained_new_2 = visualize_input;
    visualize_input2D_constrained_new_3 = visualize_input;
    visualize_input3D = visualize_input;
    visualize_inputCopy = visualize_input;
    visualize_inputPerturbed = visualize_input;
    visualize_inputReflected = visualize_input;
    visualize_inputHybrid = visualize_input;
    
    visualize_inputPerturbed.inertial_properties_noisy = inertial_properties_noisy_copy;
    visualize_inputPerturbed.inertial_properties_offset = inertial_properties_offset_copy;
    [visualize_inputPerturbed, local_offset_copy, local_angular_offset_copy, local_offset_vector_copy] = compute_localOffset(visualize_inputPerturbed);
    visualize_inputPerturbed.local_offset = local_offset_copy;
    visualize_inputPerturbed.local_angular_offset = local_angular_offset_copy;
    visualize_inputPerturbed.local_offset_vector = local_offset_vector_copy;

    eps_start = 1.0;
    eps_start_vector2D = [1.0; 1.0];
    eps_start_vector3D = [1.0; 1.0; 1.0];
    options = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', 'Display', 'iter-detailed');
    func = @(eps)compute_EpsilonOffset(eps, visualize_input);
    funcOne = @(eps)compute_First_EpsilonOffset(eps, visualize_inputOne);
    func2D = @(eps_vector2D)compute_EpsilonOffset_2D(eps_vector2D, visualize_input2D);
    func3D = @(eps_vector3D)compute_EpsilonOffset_3D(eps_vector3D, visualize_input3D);
    
    func2DConstrained_1 = @(eps_vector2D)compute_EpsilonOffset_2D_modified(eps_vector2D, visualize_input2D_constrained_new_1);
    func2DConstrained_2 = @(eps_vector2D)compute_EpsilonOffset_2D_modified(eps_vector2D, visualize_input2D_constrained_new_2);
    func2DConstrained_3 = @(eps_vector2D)compute_EpsilonOffset_2D_modified(eps_vector2D, visualize_input2D_constrained_new_3);
    
    eps_star = lsqnonlin(func, eps_start, [], [], options);
    eps_starOne = lsqnonlin(funcOne, eps_start, [], [], options);
    eps_star2D = lsqnonlin(func2D, eps_start_vector2D, [], [], options);
    eps_star3D = lsqnonlin(func3D, eps_start_vector3D, [], [], options);
    
    output_eps_star = compute_total_EpsilonOffset(eps_star, visualize_input);
    visualize_input.output_eps_star = output_eps_star;
    visualize_input = compute_Final_localOffset(visualize_input);
    
    output_eps_starOne = compute_total_EpsilonOffset(eps_starOne, visualize_inputOne);
    visualize_inputOne.output_eps_star = output_eps_starOne;
    visualize_inputOne = compute_Final_localOffset(visualize_inputOne);
    
    output_eps_star2D = compute_total_EpsilonOffset_2D(eps_star2D, visualize_input2D);
    visualize_input2D.output_eps_star = output_eps_star2D;
    visualize_input2D = compute_Final_localOffset(visualize_input2D);
    
    output_eps_star3D = compute_total_EpsilonOffset_3D(eps_star3D, visualize_input3D);
    visualize_input3D.output_eps_star = output_eps_star3D;
    visualize_input3D = compute_Final_localOffset(visualize_input3D);
    
    delta = 0.5;
    eps_delta = delta*eps_star + (1-delta)*eps_starOne;
    if eps_starOne > eps_star
        lb = eps_delta;
        ub = [];
    else
        lb = [];
        ub = eps_delta;
    end
    eps_starDelta = lsqnonlin(func, eps_start, lb, ub, options);
    output_eps_starDelta = compute_total_EpsilonOffset(eps_starDelta, visualize_inputDelta);
    visualize_inputDelta.output_eps_star = output_eps_starDelta;
    visualize_inputDelta = compute_Final_localOffset(visualize_inputDelta);
    
    options_nlcon = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    
    nonlcon_1 = @(eps_vector2D)compute_First_EpsilonOffset_2D_modified(eps_vector2D, visualize_input2D_constrained_new_1, output_eps_star.offset.PrincipalMOI(1));
    [eps_star2D_constrained_new_1, fval_nonlcon_1, exitflag_nonlcon_1, output_nonlcon_1] = fmincon(func2DConstrained_1, eps_start_vector2D, [], [], [], [], [], [], nonlcon_1, options_nlcon);
    disp(fval_nonlcon_1);
    output_eps_star2D_constrained_new_1 = compute_total_EpsilonOffset_2D(eps_star2D_constrained_new_1, visualize_input2D_constrained_new_1);
    visualize_input2D_constrained_new_1.output_eps_star = output_eps_star2D_constrained_new_1;
    visualize_input2D_constrained_new_1 = compute_Final_localOffset(visualize_input2D_constrained_new_1);
    
    nonlcon_2 = @(eps_vector2D)compute_Second_EpsilonOffset_2D_modified(eps_vector2D, visualize_input2D_constrained_new_2, output_eps_star.offset.PrincipalMOI(2));
    [eps_star2D_constrained_new_2, fval_nonlcon_2, exitflag_nonlcon_2, output_nonlcon_2] = fmincon(func2DConstrained_2, eps_start_vector2D, [], [], [], [], [], [], nonlcon_2, options_nlcon);
    disp(fval_nonlcon_2);
    output_eps_star2D_constrained_new_2 = compute_total_EpsilonOffset_2D(eps_star2D_constrained_new_2, visualize_input2D_constrained_new_2);
    visualize_input2D_constrained_new_2.output_eps_star = output_eps_star2D_constrained_new_2;
    visualize_input2D_constrained_new_2 = compute_Final_localOffset(visualize_input2D_constrained_new_2);
    
    nonlcon_3 = @(eps_vector2D)compute_Third_EpsilonOffset_2D_modified(eps_vector2D, visualize_input2D_constrained_new_3, output_eps_star.offset.PrincipalMOI(3));
    [eps_star2D_constrained_new_3, fval_nonlcon_3, exitflag_nonlcon_3, output_nonlcon_3] = fmincon(func2DConstrained_3, eps_start_vector2D, [], [], [], [], [], [], nonlcon_3, options_nlcon);
    disp(fval_nonlcon_3);
    output_eps_star2D_constrained_new_3 = compute_total_EpsilonOffset_2D(eps_star2D_constrained_new_3, visualize_input2D_constrained_new_3);
    visualize_input2D_constrained_new_3.output_eps_star = output_eps_star2D_constrained_new_3;
    visualize_input2D_constrained_new_3 = compute_Final_localOffset(visualize_input2D_constrained_new_3);
    
    lambda_start = nan(3,1);
    lambda_target = nan(3,1);
    vectors_start = nan(3,3);
    for j=1:3
        lambda_start(j) = inertial_properties_noisy_copy.PrincipalMOI(j).value;
        lambda_target(j) = inertial_properties_reference.PrincipalMOI(j).value;
        vectors_start(:,j) = inertial_properties_noisy_copy.PrincipalAxis(j).vector;
    end
    [vectors_final, lambda_final, inertiaTensor_final, iterate_data] = iterative_perturbation(vectors_start, lambda_start, lambda_target);
    m0 = ones(N,1);
    func = @(m)nLsysMD(m, mass_distribution_noisy, inertiaTensor_final, N);
    funcRef = @(m)nLsysMD(m, mass_distribution_reference, inertial_properties_reference.MOI, N);
    funcSum = @(m)nLsysMDsum(m, mass_distribution_noisy, inertiaTensor_final, N);
    cm_nonlcon = @(m)cm_constraint(m, mass_distribution_noisy, inertial_properties_noisy, N);
    
    mass_tol = 10^-10;
    cm_reflect_bool = false;
    cm_reflect_bool_eps = false;
    if all(abs(fsolve(funcRef, m0) - m0) < mass_tol*ones(N,1))
        try
            m_star = fsolve(func, m0);
            cm_constraint_val = cm_constraint(m_star, mass_distribution_noisy, inertial_properties_noisy, N);
            if cm_constraint_val > 0
                m_star_cm = fmincon(funcSum, m_star, [], [], [], [], zeros(N,1), [], cm_nonlcon, options_nlcon);
                cm_constraint_val_new = cm_constraint(m_star_cm, mass_distribution_noisy, inertial_properties_noisy, N);
                if cm_constraint_val_new <= 0.1
                    m_star = m_star_cm;
                else
                    m_star = m0;
                end
            end
        catch
            m_star = m0;
        end
    else
        m_star = m0;
    end
    
    cm_vector_reflected = cm_reflect(m_star, mass_distribution_noisy, inertial_properties_noisy, N);
    cm_vector_noisy = inertial_properties_noisy.CM.vector;
    if cm_vector_reflected(3) < cm_vector_noisy(3)
        cm_reflect_bool = true;
    end
    
    m_eps = nan(N,1);
    for m_idx=1:N
        m_eps(m_idx) = visualize_input.output_eps_star.mass_distribution(m_idx).mass;
    end 
    cm_vector_reference = inertial_properties_reference.CM.vector;
    cm_vector_reflected_eps = cm_reflect(m_eps, mass_distribution_noisy, inertial_properties_noisy, N);
    if cm_vector_reflected_eps(3) < cm_vector_noisy(3)
        CMnorm_new_offset_reflected(m) = norm(cm_vector_reflected_eps - cm_vector_reference);
        cm_reflect_bool_eps = true;
    else
        CMnorm_new_offset_reflected(m) = norm(visualize_input.output_eps_star.offset.CMVector);
    end
    
    disp(m_star);
    mass_distribution_perturbed = mass_distribution_noisy;
    for i=1:N
        mass_distribution_perturbed(i).mass = m_star(i);
        xyz_points(i,:) = transpose(mass_distribution_reference(i).vector);
    end
    
    ptCloud = pointCloud(xyz_points);
    [plane_model, plane_inlierIndices, plane_outlierIndices] = pcfitplane(ptCloud, plane_distThreshold);
    if isempty(plane_outlierIndices)
        plane_bool(m) = true;
    end
    
    [inertial_properties_perturbed, mass_distribution_perturbed] = compute_inertial_properties(mass_distribution_perturbed, N);
    inertial_properties_perturbed = compute_inertial_properties_reflected(cm_reflect_bool, cm_vector_reflected, inertial_properties_perturbed);
    inertial_properties_perturbed = add_inertial_properties_perturbation(vectors_final, inertiaTensor_final, inertial_properties_perturbed);
    inertial_properties_perturbed_offset = calculate_offset(inertial_properties_reference, inertial_properties_perturbed);
    visualize_inputPerturbed.mass_distribution_perturbed = mass_distribution_perturbed;
    visualize_inputPerturbed.inertial_properties_perturbed = inertial_properties_perturbed;
    visualize_inputPerturbed = compute_Final_localOffset_perturbed(visualize_inputPerturbed);
    
    output_eps_star_reflected = output_eps_star;
    inertial_properties_reflected = compute_inertial_properties_reflected(cm_reflect_bool_eps, cm_vector_reflected_eps, output_eps_star_reflected.inertial_properties);
    output_eps_star_reflected.inertial_properties = inertial_properties_reflected;
    visualize_inputReflected.output_eps_star = output_eps_star_reflected;
    visualize_inputReflected = compute_Final_localOffset(visualize_inputReflected);
    
    inertial_properties_optimized = visualize_input.output_eps_star.inertial_properties;
    inertial_properties_hybrid = add_inertial_properties_hybrid(inertial_properties_optimized, inertial_properties_perturbed);
    inertial_properties_hybrid_offset = calculate_offset(inertial_properties_reference, inertial_properties_hybrid);
    visualize_inputHybrid.mass_distribution_hybrid = mass_distribution_perturbed;
    visualize_inputHybrid.inertial_properties_hybrid = inertial_properties_hybrid;
    visualize_inputHybrid = compute_Final_localOffset_hybrid(visualize_inputHybrid);
    
    reference_pose = struct;
    reference_pose.T = inertial_properties_reference.T;
    reference_pose.R = inertial_properties_reference.R;
    
    optimized_pose = struct;
    optimized_pose.T = visualize_input.output_eps_star.inertial_properties.T;
    optimized_pose.R = visualize_input.output_eps_star.inertial_properties.R;
    
    optimized_pose_reflected = struct;
    optimized_pose_reflected.T = inertial_properties_reflected.T;
    optimized_pose_reflected.R = inertial_properties_reflected.R;
    
    optimized_pose_perturbed = struct;
    optimized_pose_perturbed.T = inertial_properties_perturbed.T;
    optimized_pose_perturbed.R = inertial_properties_perturbed.R;
    
    optimized_pose_hybrid = struct;
    optimized_pose_hybrid.T = inertial_properties_hybrid.T;
    optimized_pose_hybrid.R = inertial_properties_hybrid.R;
    
    anatomical_landmarks_thigh_data = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose);
    visualize_input.output_eps_star.anatomical_landmarks_thigh_data = anatomical_landmarks_thigh_data;
    
    anatomical_landmarks_thigh_data_perturbed = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose_perturbed);
    visualize_inputPerturbed.anatomical_landmarks_thigh_data = anatomical_landmarks_thigh_data_perturbed;
    
    anatomical_landmarks_thigh_data_reflected = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose_reflected);
    visualize_inputReflected.anatomical_landmarks_thigh_data = anatomical_landmarks_thigh_data_reflected;
    
    anatomical_landmarks_thigh_data_hybrid = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose_hybrid);
    visualize_inputHybrid.anatomical_landmarks_thigh_data = anatomical_landmarks_thigh_data_hybrid;
    
    for j=1:N
        global_offset_data(m, j, 1) = visualize_input.output_eps_star.mass_distribution(j).global_offset;
        global_offset_dataOne(m, j, 1) = visualize_inputOne.output_eps_star.mass_distribution(j).global_offset;
        global_offset_dataDelta(m, j, 1) = visualize_inputDelta.output_eps_star.mass_distribution(j).global_offset;
        global_offset_data2D(m, j, 1) = visualize_input2D.output_eps_star.mass_distribution(j).global_offset;
        global_offset_data2D_constrained_new_1(m, j, 1) = visualize_input2D_constrained_new_1.output_eps_star.mass_distribution(j).global_offset;
        global_offset_data2D_constrained_new_2(m, j, 1) = visualize_input2D_constrained_new_2.output_eps_star.mass_distribution(j).global_offset;
        global_offset_data2D_constrained_new_3(m, j, 1) = visualize_input2D_constrained_new_3.output_eps_star.mass_distribution(j).global_offset;
        global_offset_data3D(m, j, 1) = visualize_input3D.output_eps_star.mass_distribution(j).global_offset;
        global_offset_data_perturbed(m, j, 1) = visualize_inputPerturbed.mass_distribution_perturbed(j).global_offset;
        global_offset_data_reflected(m, j, 1) = visualize_inputReflected.output_eps_star.mass_distribution(j).global_offset;
        global_offset_data_hybrid(m, j, 1) = visualize_inputHybrid.mass_distribution_hybrid(j).global_offset;
        global_offset_data_noisy(m, j, 1) = visualize_input.mass_distribution_noisy(j).global_offset;
        
        fixed_noise_array(m, j) = mass_distribution_noisy(j).noise_idx;
    end
    
    AL_thigh = struct;
    for j=1:Num_AL_Thigh
        global_offset_data_AL(m, j, 1) = visualize_input.output_eps_star.anatomical_landmarks_thigh_data(j).global_offset;
        global_offset_data_perturbed_AL(m, j, 1) = visualize_inputPerturbed.anatomical_landmarks_thigh_data(j).global_offset;
        global_offset_data_reflected_AL(m, j, 1) = visualize_inputReflected.anatomical_landmarks_thigh_data(j).global_offset;
        global_offset_data_hybrid_AL(m, j, 1) = visualize_inputHybrid.anatomical_landmarks_thigh_data(j).global_offset;
    end
    
    all_global_types = ["Unconstrained 1D";
                        "Unconstrained 2D";
                        "Constrained 2D First";
                        "Constrained 2D Second";
                        "Constrained 2D Third";
                        "Unoptimized";
                        "Perturbation Theory";
                        "Unconstrained 1D reflected";
                        "Hybrid"];
    
    H = horzcat(global_offset_data(m, :)', global_offset_data2D(m, :)', global_offset_data2D_constrained_new_1(m, :)', global_offset_data2D_constrained_new_2(m, :)', global_offset_data2D_constrained_new_3(m, :)', global_offset_data_noisy(m, :)', global_offset_data_perturbed(m, :)', global_offset_data_reflected(m, :)', global_offset_data_hybrid(m, :)');
    [H_min, I_min] = min(H,[],2);
    I_min_val = unique(I_min);
    I_min_size = size(I_min_val);
    if I_min_size(1) == 1
        if max(H(:,I_min_val)) < 1
            BestAlgorithm(m) = all_global_types(I_min_val);
            BestConfiguration.minMOI(m) = min([inertial_properties_reference.PrincipalMOI(1).value; inertial_properties_reference.PrincipalMOI(2).value; inertial_properties_reference.PrincipalMOI(3).value]);
            BestConfiguration.normMOI(m) = inertial_properties_reference.PrincipalMOINorm;
            BestConfiguration.firstRatioMOI(m) = inertial_properties_reference.PrincipalMOI(1).value/inertial_properties_reference.PrincipalMOINorm;
            BestConfiguration.secondRatioMOI(m) = inertial_properties_reference.PrincipalMOI(2).value/inertial_properties_reference.PrincipalMOINorm;
            BestConfiguration.thirdRatioMOI(m) = inertial_properties_reference.PrincipalMOI(3).value/inertial_properties_reference.PrincipalMOINorm;
            
            BestConfiguration.normMOIdiff(m) = inertial_properties_offset.PrincipalMOINorm;
            BestConfiguration.firstMOIdiff(m) = inertial_properties_offset.PrincipalMOI(1);
            BestConfiguration.secondMOIdiff(m) = inertial_properties_offset.PrincipalMOI(2);
            BestConfiguration.thirdMOIdiff(m) = inertial_properties_offset.PrincipalMOI(3);
        else
            BestAlgorithm(m) = "None";
            BestConfiguration.minMOI(m) = min([inertial_properties_reference.PrincipalMOI(1).value; inertial_properties_reference.PrincipalMOI(2).value; inertial_properties_reference.PrincipalMOI(3).value]);
            BestConfiguration.normMOI(m) = inertial_properties_reference.PrincipalMOINorm;
            BestConfiguration.firstRatioMOI(m) = inertial_properties_reference.PrincipalMOI(1).value/inertial_properties_reference.PrincipalMOINorm;
            BestConfiguration.secondRatioMOI(m) = inertial_properties_reference.PrincipalMOI(2).value/inertial_properties_reference.PrincipalMOINorm;
            BestConfiguration.thirdRatioMOI(m) = inertial_properties_reference.PrincipalMOI(3).value/inertial_properties_reference.PrincipalMOINorm;
           
            BestConfiguration.normMOIdiff(m) = inertial_properties_offset.PrincipalMOINorm;
            BestConfiguration.firstMOIdiff(m) = inertial_properties_offset.PrincipalMOI(1);
            BestConfiguration.secondMOIdiff(m) = inertial_properties_offset.PrincipalMOI(2);
            BestConfiguration.thirdMOIdiff(m) = inertial_properties_offset.PrincipalMOI(3);
        end
    else
        BestAlgorithm(m) = "None";
    end
    
    if max(global_offset_data_perturbed(m,:)) > 3.5
        bad_list = [bad_list;m];
        visualize_inputPerturbed.T_check = visualize_input.output_eps_star.inertial_properties.T;
        visualize_inputPerturbed.R_check = visualize_input.output_eps_star.inertial_properties.R;
        visualize_inputPerturbed = compute_Final_localOffset_perturbed(visualize_inputPerturbed);
    end

    current_PCT_performance = evaluate_PCT(output_eps_star, inertial_properties_offset);
    
    angleOne_old_offset(m) = current_PCT_performance.angle.oldList(1);
    angleTwo_old_offset(m) = current_PCT_performance.angle.oldList(2);
    angleThree_old_offset(m) = current_PCT_performance.angle.oldList(3);
    
    angleOne_new_offset(m) = current_PCT_performance.angle.newList(1);
    angleTwo_new_offset(m) = current_PCT_performance.angle.newList(2);
    angleThree_new_offset(m) = current_PCT_performance.angle.newList(3);
    
    CMnorm_old_offset(m) = current_PCT_performance.CM_distance.oldNorm;
    CMnorm_new_offset(m) = current_PCT_performance.CM_distance.newNorm;
    
    current_PCT_performanceOne = evaluate_PCT(output_eps_starOne, inertial_properties_offset);
    angleOne_new_offsetOne(m) = current_PCT_performanceOne.angle.newList(1);
    angleTwo_new_offsetOne(m) = current_PCT_performanceOne.angle.newList(2);
    angleThree_new_offsetOne(m) = current_PCT_performanceOne.angle.newList(3);
    CMnorm_new_offsetOne(m) = current_PCT_performanceOne.CM_distance.newNorm;
    
    current_PCT_performanceDelta = evaluate_PCT(output_eps_starDelta, inertial_properties_offset);
    angleOne_new_offsetDelta(m) = current_PCT_performanceDelta.angle.newList(1);
    angleTwo_new_offsetDelta(m) = current_PCT_performanceDelta.angle.newList(2);
    angleThree_new_offsetDelta(m) = current_PCT_performanceDelta.angle.newList(3);
    CMnorm_new_offsetDelta(m) = current_PCT_performanceDelta.CM_distance.newNorm;
    
    current_PCT_performance2D = evaluate_PCT(output_eps_star2D, inertial_properties_offset);
    angleOne_new_offset2D(m) = current_PCT_performance2D.angle.newList(1);
    angleTwo_new_offset2D(m) = current_PCT_performance2D.angle.newList(2);
    angleThree_new_offset2D(m) = current_PCT_performance2D.angle.newList(3);
    CMnorm_new_offset2D(m) = current_PCT_performance2D.CM_distance.newNorm;
    
    current_PCT_performance2D_constrained_new_1 = evaluate_PCT(output_eps_star2D_constrained_new_1, inertial_properties_offset);
    angleOne_new_offset2D_constrained_new_1(m) = current_PCT_performance2D_constrained_new_1.angle.newList(1);
    angleTwo_new_offset2D_constrained_new_1(m) = current_PCT_performance2D_constrained_new_1.angle.newList(2);
    angleThree_new_offset2D_constrained_new_1(m) = current_PCT_performance2D_constrained_new_1.angle.newList(3);
    CMnorm_new_offset2D_constrained_new_1(m) = current_PCT_performance2D_constrained_new_1.CM_distance.newNorm;
    
    current_PCT_performance2D_constrained_new_2 = evaluate_PCT(output_eps_star2D_constrained_new_2, inertial_properties_offset);
    angleOne_new_offset2D_constrained_new_2(m) = current_PCT_performance2D_constrained_new_2.angle.newList(1);
    angleTwo_new_offset2D_constrained_new_2(m) = current_PCT_performance2D_constrained_new_2.angle.newList(2);
    angleThree_new_offset2D_constrained_new_2(m) = current_PCT_performance2D_constrained_new_2.angle.newList(3);
    CMnorm_new_offset2D_constrained_new_2(m) = current_PCT_performance2D_constrained_new_2.CM_distance.newNorm;
    
    current_PCT_performance2D_constrained_new_3 = evaluate_PCT(output_eps_star2D_constrained_new_3, inertial_properties_offset);
    angleOne_new_offset2D_constrained_new_3(m) = current_PCT_performance2D_constrained_new_3.angle.newList(1);
    angleTwo_new_offset2D_constrained_new_3(m) = current_PCT_performance2D_constrained_new_3.angle.newList(2);
    angleThree_new_offset2D_constrained_new_3(m) = current_PCT_performance2D_constrained_new_3.angle.newList(3);
    CMnorm_new_offset2D_constrained_new_3(m) = current_PCT_performance2D_constrained_new_3.CM_distance.newNorm;
    
    current_PCT_performance3D = evaluate_PCT(output_eps_star3D, inertial_properties_offset);
    angleOne_new_offset3D(m) = current_PCT_performance3D.angle.newList(1);
    angleTwo_new_offset3D(m) = current_PCT_performance3D.angle.newList(2);
    angleThree_new_offset3D(m) = current_PCT_performance3D.angle.newList(3);
    CMnorm_new_offset3D(m) = current_PCT_performance3D.CM_distance.newNorm;
    
    current_PCT_performancePerturbed = evaluate_PCT_perturbation(inertial_properties_perturbed_offset, inertial_properties_offset_copy);
    angleOne_new_offsetPerturbed(m) = current_PCT_performancePerturbed.angle.newList(1);
    angleTwo_new_offsetPerturbed(m) = current_PCT_performancePerturbed.angle.newList(2);
    angleThree_new_offsetPerturbed(m) = current_PCT_performancePerturbed.angle.newList(3);
    CMnorm_new_offsetPerturbed(m) = current_PCT_performancePerturbed.CM_distance.newNorm;
    
    current_PCT_performanceReflected = evaluate_PCT(output_eps_star_reflected, inertial_properties_offset);
    angleOne_new_offsetReflected(m) = current_PCT_performanceReflected.angle.newList(1);
    angleTwo_new_offsetReflected(m) = current_PCT_performanceReflected.angle.newList(2);
    angleThree_new_offsetReflected(m) = current_PCT_performanceReflected.angle.newList(3);
    CMnorm_new_offsetReflected(m) = current_PCT_performanceReflected.CM_distance.newNorm;
    
    current_PCT_performanceHybrid = evaluate_PCT_hybrid(inertial_properties_hybrid_offset, inertial_properties_offset);
    angleOne_new_offsetHybrid(m) = current_PCT_performanceHybrid.angle.newList(1);
    angleTwo_new_offsetHybrid(m) = current_PCT_performanceHybrid.angle.newList(2);
    angleThree_new_offsetHybrid(m) = current_PCT_performanceHybrid.angle.newList(3);
    CMnorm_new_offsetHybrid(m) = current_PCT_performanceHybrid.CM_distance.newNorm;
    
    if current_PCT_performance.CM_distance.newNorm < 1.0 && visualize_good_results
        disp("hello");
        visualize_input.visualize_flag = true;
        visualize_distribution(visualize_input);
    end
    
    visualize_distribution(visualize_input);
    
    eigenvalueChange_old_offset(m) = (abs(inertial_properties_offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offset(m) = (abs(output_eps_star.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offsetOne(m) = (abs(output_eps_starOne.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offsetDelta(m) = (abs(output_eps_starDelta.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offset2D(m) = (abs(output_eps_star2D.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offset2D_constrained_new_1(m) = (abs(output_eps_star2D_constrained_new_1.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offset2D_constrained_new_2(m) = (abs(output_eps_star2D_constrained_new_2.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offset2D_constrained_new_3(m) = (abs(output_eps_star2D_constrained_new_3.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offset3D(m) = (abs(output_eps_star3D.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offsetPerturbed(m) = (abs(inertial_properties_perturbed.PrincipalMOINorm - inertial_properties_reference.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offsetReflected(m) = (abs(output_eps_star_reflected.offset.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
    eigenvalueChange_new_offsetHybrid(m) = (abs(inertial_properties_hybrid.PrincipalMOINorm - inertial_properties_reference.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
end

folderName_str = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Images\Individual Presentation 27 Images\PerturbationConstrained_copy_" + string(N) + "\";

fileName_angle_str = folderName_str + "Angle_offset.png";
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(3,1,1)
yline(2, 'black--', 'DisplayName', 'Threshold')
hold on;
plot(angleOne_old_offset, 'r', 'DisplayName', 'Non Optimized');
plot(angleOne_new_offset, 'b', 'DisplayName', 'Optimized');
plot(angleOne_new_offsetReflected, 'g', 'DisplayName', 'Optimized Reflected');
plot(angleOne_new_offsetPerturbed, 'cyan', 'DisplayName', 'Perturbation Theory');
plot(angleOne_new_offsetHybrid, 'k', 'DisplayName', 'Hybrid');
grid;
legend;
title("First Axis Angle Offset");

subplot(3,1,2)
yline(2, 'black--', 'DisplayName', 'Threshold')
hold on;
plot(angleTwo_old_offset, 'r', 'DisplayName', 'Non Optimized');
plot(angleTwo_new_offset, 'b', 'DisplayName', 'Optimized');
plot(angleTwo_new_offsetReflected, 'g', 'DisplayName', 'Optimized Reflected');
plot(angleTwo_new_offsetPerturbed, 'cyan', 'DisplayName', 'Perturbation Theory');
plot(angleTwo_new_offsetHybrid, 'k', 'DisplayName', 'Hybrid');
grid;
title("Second Axis Angle Offset");
legend;

subplot(3,1,3)
yline(2, 'black--', 'DisplayName', 'Threshold')
hold on;
plot(angleThree_old_offset, 'r', 'DisplayName', 'Non Optimized');
plot(angleThree_new_offset, 'b', 'DisplayName', 'Optimized');
plot(angleThree_new_offsetReflected, 'g', 'DisplayName', 'Optimized Reflected');
plot(angleThree_new_offsetPerturbed, 'cyan', 'DisplayName', 'Perturbation Theory');
plot(angleThree_new_offsetHybrid, 'k', 'DisplayName', 'Hybrid');
grid;
title("Third Axis Angle Offset");
legend;

suptitle("Number of Markers = " + string(N));
saveas(gcf, fileName_angle_str);
close(gcf);

fileName_cm_graph_str = folderName_str + "CM_Graph.png";
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(2,1,1)
yline(3, 'black--', 'DisplayName', 'PCT Paper Approx Max')
hold on;
plot(CMnorm_old_offset, 'r', 'DisplayName', 'Non Optimized');
plot(CMnorm_new_offset, 'b', 'DisplayName', 'Optimized');
plot(CMnorm_new_offset_reflected, 'g', 'DisplayName', 'Optimized Reflected');
plot(CMnorm_new_offsetPerturbed, 'cyan', 'DisplayName', 'Perturbation Theory');
plot(CMnorm_new_offsetHybrid, 'k', 'DisplayName', 'Hybrid');
grid;
title("CM Location Offset");
legend;

subplot(2,1,2)
yline(20, 'black--', 'DisplayName', 'PCT Paper Max')
hold on;
scatter(trial_idx_list, eigenvalueChange_old_offset, 1, 'r', 'DisplayName', 'Non Optimized');
scatter(trial_idx_list, eigenvalueChange_new_offset, 1, 'b', 'DisplayName', 'Optimized');
scatter(trial_idx_list, eigenvalueChange_new_offsetReflected, 1, 'g', 'DisplayName', 'Optimized Reflected');
scatter(trial_idx_list, eigenvalueChange_new_offsetPerturbed, 1, 'cyan', 'DisplayName', 'Perturbation Theory');
scatter(trial_idx_list, eigenvalueChange_new_offsetHybrid, 1, 'k', 'DisplayName', 'Hybrid');
grid;
title("Eigen Value Norm % Change");
legend;

saveas(gcf, fileName_cm_graph_str);
close(gcf);

visualize_EpsilonTransform(visualize_input);
visualize_distribution(visualize_input);

visualize_EpsilonTransform(visualize_inputOne);
visualize_distribution(visualize_inputOne);

visualize_EpsilonTransform(visualize_inputDelta);
visualize_distribution(visualize_inputDelta);

visualize_EpsilonTransform(visualize_input2D);
visualize_distribution(visualize_input2D);

visualize_EpsilonTransform(visualize_input3D);
visualize_distribution(visualize_input3D);

T = table(BestAlgorithm);
G = groupsummary(T,"BestAlgorithm");
disp(G);

BestConfiguration.name = BestAlgorithm;
T_updated = struct2table(BestConfiguration);

global_offset_data_max = nan(NumTrials,1);
global_offset_data2D_max = nan(NumTrials,1);
global_offset_data2D_constrained_new_1_max = nan(NumTrials,1);
global_offset_data2D_constrained_new_2_max = nan(NumTrials,1);
global_offset_data2D_constrained_new_3_max = nan(NumTrials,1);
global_offset_data_perturbed_max = nan(NumTrials,1);
global_offset_data_reflected_max = nan(NumTrials,1);
global_offset_data_hybrid_max = nan(NumTrials,1);

std_global_offset_data = nan(NumTrials, 1);
std_global_offset_data_perturbed = nan(NumTrials, 1);
std_global_offset_data_reflected = nan(NumTrials, 1);
std_global_offset_data_hybrid = nan(NumTrials, 1);

for k=1:NumTrials
    global_offset_data_max(k,1) = max(global_offset_data(k,:));
    global_offset_data2D_max(k,1) = max(global_offset_data2D(k,:));
    global_offset_data2D_constrained_new_1_max(k,1) = max(global_offset_data2D_constrained_new_1(k,:));
    global_offset_data2D_constrained_new_2_max(k,1) = max(global_offset_data2D_constrained_new_2(k,:));
    global_offset_data2D_constrained_new_3_max(k,1) = max(global_offset_data2D_constrained_new_3(k,:));
    global_offset_data_perturbed_max(k,1) = max(global_offset_data_perturbed(k,:));
    global_offset_data_reflected_max(k,1) = max(global_offset_data_reflected(k,:));
    global_offset_data_hybrid_max(k,1) = max(global_offset_data_hybrid(k,:));
    
    std_global_offset_data(k,1) = std(global_offset_data(k,:));
    std_global_offset_data_perturbed(k,1) = std(global_offset_data_perturbed(k,:));
    std_global_offset_data_reflected(k,1) = std(global_offset_data_reflected(k,:));
    std_global_offset_data_hybrid(k,1) = std(global_offset_data_hybrid(k,:));
end

nancount_pct = length(find(isnan(global_offset_data_max)));
nancount_perturbed = length(find(isnan(global_offset_data_perturbed_max)));
nancount_hybrid = length(find(isnan(global_offset_data_hybrid_max)));
disp(["Nancount PCT: ", nancount_pct]);
disp(["Nancount Perturbed: ", nancount_perturbed]);
disp(["Nancount Hybrid: ", nancount_hybrid]);

fileName_str = folderName_str + "Overall_Histogram.png";
figure('units','normalized','outerposition',[0 0 1 1]);
nbins = 20;

subplot(3,3,1);
histogram(global_offset_data_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained");
ylabel("Frequency");

subplot(3,3,2);
histogram(global_offset_data_reflected_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_reflected_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_reflected_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained: Reflected");
ylabel("Frequency");

subplot(3,3,3);
histogram(global_offset_data2D_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data2D_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data2D_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Unconstrained");

subplot(3,3,4);
histogram(global_offset_data2D_constrained_new_1_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data2D_constrained_new_1_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data2D_constrained_new_1_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Constrained First Eigenvalue");
ylabel("Frequency");

subplot(3,3,5);
histogram(global_offset_data2D_constrained_new_2_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data2D_constrained_new_2_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data2D_constrained_new_2_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Constrained Second Eigenvalue");

subplot(3,3,6);
histogram(global_offset_data2D_constrained_new_3_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data2D_constrained_new_3_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data2D_constrained_new_3_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Constrained Third Eigenvalue");

%subplot(3,3,[7,8,9]);
subplot(3,3,7);
histogram(global_offset_data_perturbed_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_perturbed_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_perturbed_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Perturbation Theory");
xlabel("Error (cm)");
ylabel("Frequency");

subplot(3,3,9);
histogram(global_offset_data_hybrid_max, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_hybrid_max), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_hybrid_max), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Hybrid");
xlabel("Error (cm)");

suptitle("Histograms of maximum global position reconstruction errors for cluster points");
saveas(gcf, fileName_str);
close(gcf);

fileName_cm_str = folderName_str + "CM_Histogram.png";
figure('units','normalized','outerposition',[0 0 1 1]);
nbins = 20;

subplot(3,3,1);
histogram(CMnorm_new_offset, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained");
ylabel("Frequency");

subplot(3,3,2);
histogram(CMnorm_new_offset_reflected, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset_reflected), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset_reflected), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained: Reflected");
ylabel("Frequency");

subplot(3,3,3);
histogram(CMnorm_new_offset2D, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset2D), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset2D), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Unconstrained");

subplot(3,3,4);
histogram(CMnorm_new_offset2D_constrained_new_1, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset2D_constrained_new_1), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset2D_constrained_new_1), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Constrained First Eigenvalue");
ylabel("Frequency");

subplot(3,3,5);
histogram(CMnorm_new_offset2D_constrained_new_2, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset2D_constrained_new_2), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset2D_constrained_new_2), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Constrained Second Eigenvalue");

subplot(3,3,6);
histogram(CMnorm_new_offset2D_constrained_new_3, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset2D_constrained_new_3), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset2D_constrained_new_3), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("2D Constrained Third Eigenvalue");

%subplot(3,3,[7,8,9]);
subplot(3,3,7);
histogram(CMnorm_new_offsetPerturbed, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetPerturbed), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetPerturbed), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Perturbation Theory");
xlabel("Error (cm)");
ylabel("Frequency");

subplot(3,3,9);
histogram(CMnorm_new_offsetHybrid, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetHybrid), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetHybrid), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Hybrid");
xlabel("Error (cm)");

suptitle("Histograms of Center of Mass offset norm from reference");
saveas(gcf, fileName_cm_str);
close(gcf);

fileName_plane_str = folderName_str + "Overall_Histogram_Plane.png";
figure('units','normalized','outerposition',[0 0 1 1]);
nbins = 20;

subplot(4,2,1);
histogram(global_offset_data_max(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_max(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_max(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained: Planar Distributions");
ylabel("Frequency");

subplot(4,2,2);
histogram(global_offset_data_max(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_max(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_max(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained: Non-Planar Distributions");

subplot(4,2,3);
histogram(global_offset_data_reflected_max(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_reflected_max(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_reflected_max(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained Reflected: Planar Distributions");
ylabel("Frequency");

subplot(4,2,4);
histogram(global_offset_data_reflected_max(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_reflected_max(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_reflected_max(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained Reflected: Non-Planar Distributions");

subplot(4,2,5);
histogram(global_offset_data_perturbed_max(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_perturbed_max(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_perturbed_max(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Perturbation Theory: Planar Distributions");
ylabel("Frequency");

subplot(4,2,6);
histogram(global_offset_data_perturbed_max(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_perturbed_max(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_perturbed_max(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Perturbation Theory: Non-Planar Distributions");

subplot(4,2,7);
histogram(global_offset_data_hybrid_max(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_hybrid_max(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_hybrid_max(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Hybrid: Planar Distributions");
xlabel("Error (cm)");
ylabel("Frequency");

subplot(4,2,8);
histogram(global_offset_data_hybrid_max(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(global_offset_data_hybrid_max(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(global_offset_data_hybrid_max(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Hybrid: Non-Planar Distributions");
xlabel("Error (cm)");

suptitle("Histograms of maximum global position reconstruction errors for cluster points");
saveas(gcf, fileName_plane_str);
close(gcf);

fileName_plane_cm_str = folderName_str + "CM_Histogram_Plane.png";
figure('units','normalized','outerposition',[0 0 1 1]);
nbins = 20;

subplot(4,2,1);
histogram(CMnorm_new_offset(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained: Planar Distributions");
ylabel("Frequency");

subplot(4,2,2);
histogram(CMnorm_new_offset(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained: Non-Planar Distributions");

subplot(4,2,3);
histogram(CMnorm_new_offset_reflected(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset_reflected(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset_reflected(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained Reflected: Planar Distributions");
ylabel("Frequency");

subplot(4,2,4);
histogram(CMnorm_new_offset_reflected(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset_reflected(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset_reflected(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("1D Unconstrained Reflected: Non-Planar Distributions");

subplot(4,2,5);
histogram(CMnorm_new_offsetPerturbed(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetPerturbed(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetPerturbed(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Perturbation Theory: Planar Distributions");
ylabel("Frequency");

subplot(4,2,6);
histogram(CMnorm_new_offsetPerturbed(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetPerturbed(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetPerturbed(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Perturbation Theory: Non-Planar Distributions");

subplot(4,2,7);
histogram(CMnorm_new_offsetHybrid(plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetHybrid(plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetHybrid(plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Hybrid: Planar Distributions");
xlabel("Error (cm)");
ylabel("Frequency");

subplot(4,2,8);
histogram(CMnorm_new_offsetHybrid(~plane_bool), nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetHybrid(~plane_bool)), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetHybrid(~plane_bool)), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
title("Hybrid: Non-Planar Distributions");
xlabel("Error (cm)");

suptitle("Histograms of Center of Mass offset norm from reference");
saveas(gcf, fileName_plane_cm_str);
close(gcf);

fileName_cm_new = folderName_str + "CM_Histogram_New.png";
figure('units','normalized','outerposition',[0 0 1 1]);
nbins = 20;

xlim_min = min([0; min(CMnorm_new_offset); min(CMnorm_new_offset_reflected); min(CMnorm_new_offsetPerturbed)]);
xlim_max = max([max(CMnorm_new_offset); max(CMnorm_new_offset_reflected); max(CMnorm_new_offsetPerturbed)]) + 1;

subplot(4,1,1);
histogram(CMnorm_new_offset, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_min, xlim_max]);
title("1D Unconstrained: Original");
ylabel("Frequency");

subplot(4,1,2);
histogram(CMnorm_new_offset_reflected, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offset_reflected), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offset_reflected), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_min, xlim_max]);
title("1D Unconstrained: Reflected");
ylabel("Frequency");

subplot(4,1,3);
histogram(CMnorm_new_offsetPerturbed, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetPerturbed), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetPerturbed), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_min, xlim_max]);
title("Perturbation Theory");
ylabel("Frequency");

subplot(4,1,4);
histogram(CMnorm_new_offsetHybrid, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(CMnorm_new_offsetHybrid), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(CMnorm_new_offsetHybrid), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_min, xlim_max]);
title("Hybrid");
xlabel("Error (cm)");
ylabel("Frequency");

suptitle("Histograms of Center of Mass offset norm from reference");
saveas(gcf, fileName_cm_new);
close(gcf);

fileName_al_str = folderName_str + "Anatomical_Landmarks_Graph.png";
figure('units','normalized','outerposition',[0 0 1 1]);

for j=1:Num_AL_Thigh
    subplot(2,2,j);
    plot(global_offset_data_AL(:, j), 'r', 'DisplayName', 'PCT Optimized');
    hold on;
    plot(global_offset_data_reflected_AL(:, j), 'g', 'DisplayName', 'PCT Optimized Reflected');
    plot(global_offset_data_perturbed_AL(:, j), 'b', 'DisplayName', 'Perturbation Theory');
    plot(global_offset_data_hybrid_AL(:, j), 'k', 'DisplayName', 'Hybrid');
    current_title_str = visualize_inputPerturbed.anatomical_landmarks_thigh_data(j).name;
    title(current_title_str);
    
    if j==2
        legend;
    end
    
    if ~(mod(j,2)==0)
        ylabel("Error (mm)");
    end
    
    if j>2
        xlabel("Configuration / Trial Number");
    end
    
end

suptitle("Reconstruction Error for the Global Positions of Markerless Femoral Anatomical Landmarks");
saveas(gcf, fileName_al_str);
close(gcf);

fileName_al_hist_str = folderName_str + "Anatomical_Landmarks_Histogram.png";
figure('units','normalized','outerposition',[0 0 1 1]);
nbins = 20;

for j=1:Num_AL_Thigh
    subplot(2,2,j);
    histogram(global_offset_data_AL(:, j), nbins, 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'None', 'DisplayName', 'PCT Optimized');
    hold on;
    histogram(global_offset_data_reflected_AL(:, j), nbins, 'FaceColor', 'g', 'FaceAlpha', 0.1, 'EdgeColor', 'None', 'DisplayName', 'PCT Optimized Reflected');
    histogram(global_offset_data_perturbed_AL(:, j), nbins, 'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'None', 'DisplayName', 'Perturbation Theory');
    histogram(global_offset_data_hybrid_AL(:, j), nbins, 'FaceColor', 'k', 'FaceAlpha', 0.1, 'EdgeColor', 'None', 'DisplayName', 'Hybrid');
    current_title_str = visualize_inputPerturbed.anatomical_landmarks_thigh_data(j).name;
    title(current_title_str);
    if j==2
        legend;
    end
end

suptitle("Reconstruction Error Histogram for the Global Positions of Markerless Femoral Anatomical Landmarks");
saveas(gcf, fileName_al_hist_str);
close(gcf);

fileName_std_str = folderName_str + "Overall_Histogram_STD.png";
figure('units','normalized','outerposition',[0 0 1 1]);
nbins = 20;

xlim_std_min = min([0; min(std_global_offset_data); min(std_global_offset_data_reflected); min(std_global_offset_data_perturbed)]);
xlim_std_max = max([max(std_global_offset_data); max(std_global_offset_data_reflected); max(std_global_offset_data_perturbed)]);

subplot(4,1,1);
histogram(std_global_offset_data, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(std_global_offset_data), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(std_global_offset_data), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_std_min, xlim_std_max]);
title("1D Unconstrained");
ylabel("Frequency");

subplot(4,1,2);
histogram(std_global_offset_data_reflected, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(std_global_offset_data_reflected), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(std_global_offset_data_reflected), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_std_min, xlim_std_max]);
title("1D Unconstrained: Reflected");
ylabel("Frequency");

subplot(4,1,3);
histogram(std_global_offset_data_perturbed, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(std_global_offset_data_perturbed), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(std_global_offset_data_perturbed), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_std_min, xlim_std_max]);
title("Perturbation Theory");

subplot(4,1,4);
histogram(std_global_offset_data_hybrid, nbins, 'DisplayName', 'Histogram');
hold on;
xline(nanmean(std_global_offset_data_hybrid), 'k', 'LineWidth', 1.5, 'DisplayName', 'Mean');
xline(nanmedian(std_global_offset_data_hybrid), 'k--', 'LineWidth', 1.5, 'DisplayName', 'Median');
legend;
xlim([xlim_std_min, xlim_std_max]);
title("Hybrid");
xlabel("Error (cm)");

suptitle("Histograms of standard deviation in global position reconstruction errors for cluster points");
saveas(gcf, fileName_std_str);
close(gcf);


function visualize_input = compute_Final_localOffset(visualize_input)
    T = visualize_input.output_eps_star.inertial_properties.T;
    R = visualize_input.output_eps_star.inertial_properties.R;
    N = visualize_input.N;
    for k=1:N
        current_global_vector = visualize_input.output_eps_star.mass_distribution(k).vector;
        current_ref_vector = visualize_input.mass_distribution_reference(k).local_displacement;
        original_global_vector = visualize_input.mass_distribution_reference(k).vector;
        current_local_vector = global2local(T, R, current_global_vector);
        global_vector_reconstructed = local2global(T, R, current_ref_vector);
        global_offset = norm(global_vector_reconstructed - original_global_vector);
        visualize_input.output_eps_star.mass_distribution(k).local_displacement = current_local_vector;
        visualize_input.output_eps_star.mass_distribution(k).local_offset = norm(current_local_vector - current_ref_vector);
        visualize_input.output_eps_star.mass_distribution(k).global_reconstructed = global_vector_reconstructed;
        visualize_input.output_eps_star.mass_distribution(k).global_offset = global_offset;
    end
    inertial_properties_noisy_reversed = axis_reversal_original(visualize_input.inertial_properties_noisy, visualize_input.inertial_properties_offset);
    visualize_input.inertial_properties_noisy = inertial_properties_noisy_reversed;
    visualize_input.inertial_properties_offset = calculate_offset(visualize_input.inertial_properties_reference, visualize_input.inertial_properties_noisy);
    T_noisy = visualize_input.inertial_properties_noisy.T;
    R_noisy = visualize_input.inertial_properties_noisy.R;
    for k=1:N
        current_ref_vector = visualize_input.mass_distribution_reference(k).local_displacement;
        original_global_vector = visualize_input.mass_distribution_reference(k).vector;
        global_vector_reconstructed_noisy = local2global(T_noisy, R_noisy, current_ref_vector);
        global_offset_noisy = norm(global_vector_reconstructed_noisy - original_global_vector);
        visualize_input.mass_distribution_noisy(k).global_reconstructed = global_vector_reconstructed_noisy;
        visualize_input.mass_distribution_noisy(k).global_offset = global_offset_noisy;
    end
end

function output_original_reversed = axis_reversal_original(output_original, inertial_properties_offset)
    orthogonal_tol = 10^(-10);
    output_original_reversed = output_original;
    reversal_flag = false;
    for m=1:3
        angleOriginal = inertial_properties_offset.angle(m);
        if angleOriginal > 90
            reversal_flag = true;
            output_original_reversed.PrincipalAxis(m).vector = -output_original.PrincipalAxis(m).vector;
            output_original_reversed.R(:,m) = -output_original.R(:,m);
        end
    end
    if reversal_flag
        if abs(det(output_original_reversed.R) - 1) < orthogonal_tol
            disp("Axis reversal successful");
        end
    end
end

function PCT_performance = evaluate_PCT(output_eps_star, inertial_properties_offset)
    PCT_performance = struct;
    PCT_performance.MOI_norm = abs(output_eps_star.offset.PrincipalMOINorm) < abs(inertial_properties_offset.PrincipalMOINorm);
    PCT_performance.CM_distance.boolVal = norm(output_eps_star.offset.CMVector) < norm(inertial_properties_offset.CMVector);
    PCT_performance.CM_distance.oldNorm = norm(inertial_properties_offset.CMVector);
    PCT_performance.CM_distance.newNorm = norm(output_eps_star.offset.CMVector);
    angleList = false(3,1);
    angleOldList = nan(3,1);
    angleNewList = nan(3,1);
    angleDiffList = nan(3,1);
    for m=1:3
        angleOld = inertial_properties_offset.angle(m);
        angleNew = output_eps_star.offset.angle(m);

        if angleOld > 90
            angleOld = 180 - inertial_properties_offset.angle(m);
        end

        if angleNew > 90
            angleNew = 180 - output_eps_star.offset.angle(m);
        end

        angleOldList(m) = angleOld;
        angleNewList(m) = angleNew;
        angleDiffList(m) = angleOld - angleNew;
        if angleNew < angleOld
            angleList(m) = true;
        end
    end
    PCT_performance.angle.boolVal = all(angleList);
    PCT_performance.angle.oldList = angleOldList;
    PCT_performance.angle.newList = angleNewList;
end

function MOI_EpsilonOffset = compute_EpsilonOffset(eps, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset(eps, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    MOI_EpsilonOffset = inertial_properties_EpsilonOffset.PrincipalMOINorm;
    First_MOI_EpsilonOffset = inertial_properties_EpsilonOffset.PrincipalMOI(1);
end

function First_MOI_EpsilonOffset = compute_First_EpsilonOffset(eps, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset(eps, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    First_MOI_EpsilonOffset = inertial_properties_EpsilonOffset.PrincipalMOI(1);
end

function MOI_EpsilonOffset_2D = compute_EpsilonOffset_2D(eps_vector, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset_2D(eps_vector, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    MOI_EpsilonOffset_2D = inertial_properties_EpsilonOffset.PrincipalMOINorm;
    First_MOI_EpsilonOffset_2D = inertial_properties_EpsilonOffset.PrincipalMOI(1);
end

function MOI_EpsilonOffset_2D_modified = compute_EpsilonOffset_2D_modified(eps_vector, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset_2D(eps_vector, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    MOI_EpsilonOffset_2D_modified = (inertial_properties_EpsilonOffset.PrincipalMOINorm)^2;
    First_MOI_EpsilonOffset_2D = inertial_properties_EpsilonOffset.PrincipalMOI(1);
end

function [c,ceq] = compute_First_EpsilonOffset_2D_modified(eps_vector, visualize_input, First_MOI_EpsilonOffset)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset_2D(eps_vector, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    First_MOI_EpsilonOffset_2D = inertial_properties_EpsilonOffset.PrincipalMOI(1);
    c = (First_MOI_EpsilonOffset_2D)^2 - (First_MOI_EpsilonOffset)^2;
    ceq = [];
end

function [c,ceq] = compute_Second_EpsilonOffset_2D_modified(eps_vector, visualize_input, Second_MOI_EpsilonOffset)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset_2D(eps_vector, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    Second_MOI_EpsilonOffset_2D = inertial_properties_EpsilonOffset.PrincipalMOI(2);
    c = (Second_MOI_EpsilonOffset_2D)^2 - (Second_MOI_EpsilonOffset)^2;
    ceq = [];
end

function [c,ceq] = compute_Third_EpsilonOffset_2D_modified(eps_vector, visualize_input, Third_MOI_EpsilonOffset)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset_2D(eps_vector, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    Third_MOI_EpsilonOffset_2D = inertial_properties_EpsilonOffset.PrincipalMOI(3);
    c = (Third_MOI_EpsilonOffset_2D)^2 - (Third_MOI_EpsilonOffset)^2;
    ceq = [];
end

function MOI_EpsilonOffset_3D = compute_EpsilonOffset_3D(eps_vector, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    eps_output = compute_total_EpsilonOffset_3D(eps_vector, visualize_input);
    inertial_properties_EpsilonOffset = eps_output.offset;
    MOI_EpsilonOffset_3D = inertial_properties_EpsilonOffset.PrincipalMOINorm;
    First_MOI_EpsilonOffset_3D = inertial_properties_EpsilonOffset.PrincipalMOI(1);
end

function eps_output = compute_total_EpsilonOffset(eps, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    local_offset = visualize_input.local_offset;
    max_local_offset = max(local_offset);
    
    inertial_properties_reference = visualize_input.inertial_properties_reference;
    mass_distribution_eps = visualize_input.mass_distribution_noisy;
    N = visualize_input.N;
    for k=1:N
        mass_distribution_eps(k).mass = eps - (local_offset(k)/max_local_offset);
    end
    [inertial_properties_eps, mass_distribution_eps] = compute_inertial_properties(mass_distribution_eps, N);
    inertial_properties_EpsilonOffset = calculate_offset(inertial_properties_reference, inertial_properties_eps);
    if axis_reversal_bool
        inertial_properties_eps_reversed = axis_reversal_original(inertial_properties_eps, inertial_properties_EpsilonOffset);
        inertial_properties_eps = inertial_properties_eps_reversed;
        inertial_properties_EpsilonOffset = calculate_offset(inertial_properties_reference, inertial_properties_eps);
    end
    
    eps_output = struct;
    eps_output.mass_distribution = mass_distribution_eps;
    eps_output.inertial_properties = inertial_properties_eps;
    eps_output.offset = inertial_properties_EpsilonOffset;
end

function eps_output = compute_total_EpsilonOffset_2D(eps_vector, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    local_offset = visualize_input.local_offset;
    local_angular_offset = visualize_input.local_angular_offset;
    max_local_offset = max(local_offset);
    max_local_angular_offset = max(local_angular_offset);
    
    eps_x = eps_vector(1);
    eps_y = eps_vector(2);
    
    local_offset_vector = visualize_input.local_offset_vector;
    local_offset_x = local_offset_vector(:,1);
    local_offset_y = local_offset_vector(:,2);
    local_offset_z = local_offset_vector(:,3);
    
    max_local_offset_x = max(local_offset_x);
    max_local_offset_y = max(local_offset_y);
    
    
    mass_redistribution = 'Angular 2D';
    
    inertial_properties_reference = visualize_input.inertial_properties_reference;
    mass_distribution_eps = visualize_input.mass_distribution_noisy;
    N = visualize_input.N;
    for k=1:N
        switch mass_redistribution
            case 'Symmetric 2D'
                mass_x = eps_x - (local_offset(k)/max_local_offset);
                mass_y = eps_y - (local_offset(k)/max_local_offset);
                mass_z = 0;
            case 'Standard 2D'
                mass_x = eps_x - (local_offset_x(k)/max_local_offset_x);
                mass_y = eps_y - (local_offset_y(k)/max_local_offset_y);
                mass_z = 0;
            case 'Absolute Standard 2D'
                mass_x = eps_x - abs(local_offset_x(k)/max_local_offset_x);
                mass_y = eps_y - abs(local_offset_y(k)/max_local_offset_y);
                mass_z = 0;
            case 'Product 2D'
                mass_x = eps_y*(eps_x - (local_offset(k)/max_local_offset));
                mass_y = 0;
                mass_z = 0;
            case 'Angular 2D'
                mass_x = eps_x - (local_offset(k)/max_local_offset);
                mass_y = eps_y - (local_angular_offset(k)/max_local_angular_offset);
                mass_z = 0;
            case 'Absolute Angular 2D'
                mass_x = eps_x - (local_offset(k)/max_local_offset);
                mass_y = eps_y - abs(local_angular_offset(k)/max_local_angular_offset);
                mass_z = 0;
        end
        mass_distribution_eps(k).mass = sqrt(mass_x^2 + mass_y^2 + mass_z^2);
        
    end
    [inertial_properties_eps, mass_distribution_eps] = compute_inertial_properties(mass_distribution_eps, N);
    inertial_properties_EpsilonOffset = calculate_offset(inertial_properties_reference, inertial_properties_eps);
    if axis_reversal_bool
        inertial_properties_eps_reversed = axis_reversal_original(inertial_properties_eps, inertial_properties_EpsilonOffset);
        inertial_properties_eps = inertial_properties_eps_reversed;
        inertial_properties_EpsilonOffset = calculate_offset(inertial_properties_reference, inertial_properties_eps);
    end
    
    eps_output = struct;
    eps_output.mass_distribution = mass_distribution_eps;
    eps_output.inertial_properties = inertial_properties_eps;
    eps_output.offset = inertial_properties_EpsilonOffset;
end

function eps_output = compute_total_EpsilonOffset_3D(eps_vector, visualize_input)
    axis_reversal_bool = visualize_input.optimized_axis_reversal_bool;
    local_offset = visualize_input.local_offset;
    max_local_offset = max(local_offset);
    
    eps_x = eps_vector(1);
    eps_y = eps_vector(2);
    eps_z = eps_vector(3);
    
    local_offset_vector = visualize_input.local_offset_vector;
    local_offset_x = local_offset_vector(:,1);
    local_offset_y = local_offset_vector(:,2);
    local_offset_z = local_offset_vector(:,3);
    
    max_local_offset_x = max(local_offset_x);
    max_local_offset_y = max(local_offset_y);
    max_local_offset_z = max(local_offset_z);
    
    
    mass_redistribution = 'Symmetric 3D';
    
    inertial_properties_reference = visualize_input.inertial_properties_reference;
    mass_distribution_eps = visualize_input.mass_distribution_noisy;
    N = visualize_input.N;
    for k=1:N
        switch mass_redistribution
            case 'Symmetric 3D'
                mass_x = eps_x - (local_offset(k)/max_local_offset);
                mass_y = eps_y - (local_offset(k)/max_local_offset);
                mass_z = eps_z - (local_offset(k)/max_local_offset);
            case 'Standard 3D'
                mass_x = eps_x - (local_offset_x(k)/max_local_offset_x);
                mass_y = eps_y - (local_offset_y(k)/max_local_offset_y);
                mass_z = eps_z - (local_offset_z(k)/max_local_offset_z);
            case 'Absolute Standard 3D'
                mass_x = eps_x - abs(local_offset_x(k)/max_local_offset_x);
                mass_y = eps_y - abs(local_offset_y(k)/max_local_offset_y);
                mass_z = eps_z - abs(local_offset_z(k)/max_local_offset_z);
        end
        mass_distribution_eps(k).mass = sqrt(mass_x^2 + mass_y^2 + mass_z^2);
        
    end
    [inertial_properties_eps, mass_distribution_eps] = compute_inertial_properties(mass_distribution_eps, N);
    inertial_properties_EpsilonOffset = calculate_offset(inertial_properties_reference, inertial_properties_eps);
    if axis_reversal_bool
        inertial_properties_eps_reversed = axis_reversal_original(inertial_properties_eps, inertial_properties_EpsilonOffset);
        inertial_properties_eps = inertial_properties_eps_reversed;
        inertial_properties_EpsilonOffset = calculate_offset(inertial_properties_reference, inertial_properties_eps);
    end
    
    eps_output = struct;
    eps_output.mass_distribution = mass_distribution_eps;
    eps_output.inertial_properties = inertial_properties_eps;
    eps_output.offset = inertial_properties_EpsilonOffset;
end

function [visualize_input, local_offset, local_angular_offset, local_offset_vector] = compute_localOffset(visualize_input)
    mass_distribution_reference = visualize_input.mass_distribution_reference;
    mass_distribution_noisy = visualize_input.mass_distribution_noisy;
    inertial_properties_reference = visualize_input.inertial_properties_reference;
    inertial_properties_noisy = visualize_input.inertial_properties_noisy;
    T_reference = inertial_properties_reference.T;
    R_reference = inertial_properties_reference.R;
    T_noisy = inertial_properties_noisy.T;
    R_noisy = inertial_properties_noisy.R;
    N = visualize_input.N;
    local_offset = nan(N,1);
    local_angular_offset = nan(N,1);
    local_offset_vector = nan(N,3);
    for k=1:N
        position_global_reference = mass_distribution_reference(k).vector;
        mass_distribution_reference(k).local_displacement = global2local(T_reference, R_reference, position_global_reference);
        position_global_noisy = mass_distribution_noisy(k).vector;
        mass_distribution_noisy(k).local_displacement = global2local(T_noisy, R_noisy, position_global_noisy);
        local_offset_vector(k,:) = mass_distribution_noisy(k).local_displacement - mass_distribution_reference(k).local_displacement;
        local_offset(k) = norm(mass_distribution_noisy(k).local_displacement - mass_distribution_reference(k).local_displacement);
        local_angular_offset(k) = rad2deg(acos(dot(mass_distribution_reference(k).local_displacement, mass_distribution_noisy(k).local_displacement)/(norm(mass_distribution_reference(k).local_displacement)*norm(mass_distribution_noisy(k).local_displacement))));
        mass_distribution_noisy(k).local_offset = local_offset(k);
        mass_distribution_noisy(k).local_angular_offset = local_angular_offset(k);
    end
    visualize_input.mass_distribution_reference = mass_distribution_reference;
    visualize_input.mass_distribution_noisy = mass_distribution_noisy;
end

function cylinder_data = generate_cylinder_data
    t = 0:0.1:50;
    n = length(t);
    np = 100;
    r = 1 + 0.05*t;
    [X, Y, Z] = cylinder(r, np-1);
    Zscaled = 50*Z;

    cylinder_data = struct;
    cylinder_data.X = X;
    cylinder_data.Y = Y;
    cylinder_data.Z = Zscaled;
    cylinder_data.np = np;
    cylinder_data.n = n;
end

function mass_distribution = generate_random_distribution(cylinder_data, N)
    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;
    np = cylinder_data.np;
    n = cylinder_data.n;
    
    dZ = 0.1;
    deltaZ_max = 6;
    deltaZ_random = 0.5;
    Nmax = deltaZ_max/dZ;
    Nrandom = deltaZ_random/dZ;
    Nsubtract = Nmax+Nrandom+1;

    plane_idx = randi([1, np], 1, N);
    height_idx = randi([Nsubtract, n-Nsubtract], 1, N);

    mass_distribution = struct;
    for k=1:N
        current_plane_idx = plane_idx(k);
        current_height_idx = height_idx(k);
        mass_distribution(k).mass = 1;
        mass_distribution(k).plane_idx = current_plane_idx;
        mass_distribution(k).height_idx = current_height_idx;
        mass_distribution(k).X = X(current_height_idx, current_plane_idx);
        mass_distribution(k).Y = Y(current_height_idx, current_plane_idx);
        mass_distribution(k).Z = Z(current_height_idx, current_plane_idx);
        mass_distribution(k).vector = [mass_distribution(k).X; mass_distribution(k).Y; mass_distribution(k).Z];
    end

end

function noisy_distribution = generate_noisy_ditribution(mass_distribution, N, cylinder_data)
    
    dZ = 0.1;
    deltaZ_max = 6;
    deltaZ_random = 0.5;
    Nmax = deltaZ_max/dZ;
    N_random = deltaZ_random/dZ;
    
    noiseModel = 'PCTmodel';
    
    switch noiseModel
        case 'RandomLevel'
            noise_idx = randi([1, Nmax], 1, N);
        case 'PCTmodel'
            systematic_noise_idx = Nmax*ones(1,N);
            random_noise_idx = randi([-N_random, N_random], 1, N);
            noise_idx = nan(1,N);
            noise_idx(1:2:end) = random_noise_idx(1:2:end);
            noise_idx(2:2:end) = systematic_noise_idx(2:2:end) + random_noise_idx(2:2:end);
        case 'PCTCapozzoHybridmodel'
            systematic_noise_idx = Nmax*ones(1,N);
            random_noise_idx = randi([-N_random, N_random], 1, N);
            noise_idx = nan(1,N);
            noise_idx(1:4) = random_noise_idx(1:4);
            noise_idx(5:6) = systematic_noise_idx(5:6) + random_noise_idx(5:6);
        case 'HighLevel'
            noise_idx = randi([round(0.7*Nmax), Nmax], 1, N);
        case 'HighLevelAlternate'
            noise_idx = randi([round(0.7*Nmax), Nmax], 1, N);
            noise_idx(1:2:end) = -noise_idx(1:2:end);
        case 'MaxLevelAlternate'
            noise_idx = Nmax*ones(1,N);
            noise_idx(1:2:end) = -noise_idx(1:2:end);
    end
    
    noisy_distribution = struct;
    
    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;
    
    for k=1:N
        current_plane_idx = mass_distribution(k).plane_idx;
        current_height_idx = mass_distribution(k).height_idx;
        current_noise_idx = noise_idx(k);
        
        noisy_height_idx = current_height_idx + current_noise_idx;
        
        noisy_distribution(k).mass = mass_distribution(k).mass;
        noisy_distribution(k).X = X(noisy_height_idx, current_plane_idx);
        noisy_distribution(k).Y = Y(noisy_height_idx, current_plane_idx);
        noisy_distribution(k).Z = Z(noisy_height_idx, current_plane_idx);
        noisy_distribution(k).vector = [noisy_distribution(k).X; noisy_distribution(k).Y; noisy_distribution(k).Z];
        noisy_distribution(k).noise_idx = current_noise_idx;
    end

end

function [inertial_properties, mass_distribution] = compute_inertial_properties(mass_distribution, N)
    localOrigin = 'cm';
    
    inertial_properties = struct;
    M_sum = 0;
    CM_X_sum = 0;
    CM_Y_sum = 0;
    CM_Z_sum = 0;
    for k=1:N
        switch localOrigin
            case 'centroidal'
                mass_multiplier = 1;
            case 'cm'
                mass_multiplier = mass_distribution(k).mass;
        end
        M_sum = M_sum + mass_multiplier;
        CM_X_sum = CM_X_sum + mass_multiplier*mass_distribution(k).X;
        CM_Y_sum = CM_Y_sum + mass_multiplier*mass_distribution(k).Y;
        CM_Z_sum = CM_Z_sum + mass_multiplier*mass_distribution(k).Z;
    end
    CM_X = CM_X_sum / M_sum;
    CM_Y = CM_Y_sum / M_sum;
    CM_Z = CM_Z_sum / M_sum;

    I_XX = 0;
    I_YY = 0;
    I_ZZ = 0;
    I_XY = 0;
    I_XZ = 0;
    I_YZ = 0;
    for k=1:N
        mass_distribution(k).X_rel = mass_distribution(k).X - CM_X;
        mass_distribution(k).Y_rel = mass_distribution(k).Y - CM_Y;
        mass_distribution(k).Z_rel = mass_distribution(k).Z - CM_Z;
        I_XX = I_XX + mass_distribution(k).mass*((mass_distribution(k).Y_rel)^2 + (mass_distribution(k).Z_rel)^2);
        I_YY = I_YY + mass_distribution(k).mass*((mass_distribution(k).X_rel)^2 + (mass_distribution(k).Z_rel)^2);
        I_ZZ = I_ZZ + mass_distribution(k).mass*((mass_distribution(k).X_rel)^2 + (mass_distribution(k).Y_rel)^2);
        I_XY = I_XY + ((-1) * mass_distribution(k).mass * mass_distribution(k).X_rel * mass_distribution(k).Y_rel);
        I_XZ = I_XZ + ((-1) * mass_distribution(k).mass * mass_distribution(k).X_rel * mass_distribution(k).Z_rel);
        I_YZ = I_YZ + ((-1) * mass_distribution(k).mass * mass_distribution(k).Y_rel * mass_distribution(k).Z_rel);
    end

    tol_orthogonal = 10^-10;
    I = [I_XX I_XY I_XZ; I_XY I_YY I_YZ; I_XZ I_YZ I_ZZ];
    [V, D] = eig(I);
    diag_MOI = diag(D);
    third_axis = V(:,3);
    
    if dot(cross(V(:,1), V(:,2)), V(:,3)) + 1 < tol_orthogonal
        V(:,3) = -third_axis;
    end
    
    inertial_properties.CM.x = CM_X;
    inertial_properties.CM.y = CM_Y;
    inertial_properties.CM.z = CM_Z;
    inertial_properties.CM.vector = [CM_X; CM_Y; CM_Z];
    inertial_properties.MOI = I;
    inertial_properties.PrincipalAxis = struct;
    inertial_properties.PrincipalAxis(1).vector = V(:,1);
    inertial_properties.PrincipalAxis(2).vector = V(:,2);
    inertial_properties.PrincipalAxis(3).vector = V(:,3);
    inertial_properties.PrincipalMOI = struct;
    inertial_properties.PrincipalMOI(1).value = diag_MOI(1);
    inertial_properties.PrincipalMOI(2).value = diag_MOI(2);
    inertial_properties.PrincipalMOI(3).value = diag_MOI(3);
    inertial_properties.PrincipalMOINorm = sqrt((diag_MOI(1)^2) + (diag_MOI(2)^2) + (diag_MOI(3)^2));
    inertial_properties.R = V;
    inertial_properties.T = inertial_properties.CM.vector;
end

function inertial_properties_reflected = compute_inertial_properties_reflected(cm_reflect_bool, cm_reflected_vector, inertial_properties)
    
    inertial_properties_reflected = inertial_properties;
    
    if cm_reflect_bool
        inertial_properties_reflected.CM.x = cm_reflected_vector(1);
        inertial_properties_reflected.CM.y = cm_reflected_vector(2);
        inertial_properties_reflected.CM.z = cm_reflected_vector(3);
        inertial_properties_reflected.CM.vector = [cm_reflected_vector(1); cm_reflected_vector(2); cm_reflected_vector(3)];
        inertial_properties_reflected.T = inertial_properties_reflected.CM.vector;
    end
end

function inertial_properties_added = add_inertial_properties_perturbation(vectors_final, inertiaTensor_final, inertial_properties_old)

    inertial_properties_added = inertial_properties_old;

    inertial_properties_added.MOI = inertiaTensor_final;
    inertial_properties_added.PrincipalAxis = struct;
    inertial_properties_added.PrincipalAxis(1).vector = vectors_final(:,1);
    inertial_properties_added.PrincipalAxis(2).vector = vectors_final(:,2);
    inertial_properties_added.PrincipalAxis(3).vector = vectors_final(:,3);
    inertial_properties_added.R = vectors_final;
end

function inertial_properties_hybrid = add_inertial_properties_hybrid(inertial_properties_optimized, inertial_properties_perturbed)
    inertial_properties_hybrid = inertial_properties_optimized;
    
    inertial_properties_hybrid.CM.x = inertial_properties_perturbed.CM.x;
    inertial_properties_hybrid.CM.y = inertial_properties_perturbed.CM.y;
    inertial_properties_hybrid.CM.z = inertial_properties_perturbed.CM.z;
    inertial_properties_hybrid.CM.vector = [inertial_properties_hybrid.CM.x; inertial_properties_hybrid.CM.y; inertial_properties_hybrid.CM.z];
    inertial_properties_hybrid.T = inertial_properties_hybrid.CM.vector;
    
end

function sanity_check_bool = sanity_check_axes(mass_distribution, inertial_properties, N)
    if N==3
        tol_orthogonal = 10^-10;
        plane_vector = mass_distribution(1).vector - mass_distribution(2).vector;
        if abs(dot(inertial_properties.PrincipalAxis(3).vector, plane_vector)) < tol_orthogonal
            disp("For N=3, Last principal axis is orthogonal to plane of triangle");
            axes_orthogonal = true;
        else
            disp("For N=3, Orthogonality test failed");
            axes_orthogonal = false;
        end
        
        [triangle_inequality_satisfied, triangle_symmetric] = apply_triangle_inequality(mass_distribution);
        triangle_conditions = triangle_inequality_satisfied && ~triangle_symmetric;
        
        if axes_orthogonal && triangle_conditions
            sanity_check_bool = true;
        else
            sanity_check_bool = false;
        end
    else
        sanity_check_bool = true;
    end
end

function [triangle_inequality_satisfied, triangle_symmetric] = apply_triangle_inequality(mass_distribution)
    triangle_inequality_satisfied = false;
    triangle_symmetric = false;
    tol = 1.0e-4;
    mass_distribution_matrix = horzcat(mass_distribution(1).vector, mass_distribution(2).vector, mass_distribution(3).vector);
    if ~any(any(isnan(mass_distribution_matrix)))
        r1 = mass_distribution(1).vector;
        r2 = mass_distribution(2).vector;
        r3 = mass_distribution(3).vector;
        r_12 = norm(r2-r1);
        r_13 = norm(r3-r1);
        r_23 = norm(r3-r2);
        sum_of_sides = r_12 + r_23;
        third_side = r_13;
        diff_side = abs(sum_of_sides - third_side);
        if diff_side > tol
            triangle_inequality_satisfied = true;
        end
        if or(or((r_12 == r_13), (r_12 == r_23)), (r_23 == r_13))
           triangle_symmetric = true;
           disp("Triangle has symmetry, PCT symmetry check failed");
        else
            disp("PCT symmetry check success");
        end
    end
end

function inertial_properties_offset = calculate_offset(inertial_properties_reference, inertial_properties_noisy)
    inertial_properties_offset = struct;
    for j=1:3
        inertial_properties_offset.angle(j) = rad2deg(acos(dot(inertial_properties_reference.PrincipalAxis(j).vector, inertial_properties_noisy.PrincipalAxis(j).vector)));
        inertial_properties_offset.PrincipalMOI(j) = inertial_properties_reference.PrincipalMOI(j).value - inertial_properties_noisy.PrincipalMOI(j).value;
    end
    inertial_properties_offset.PrincipalMOINorm = inertial_properties_reference.PrincipalMOINorm - inertial_properties_noisy.PrincipalMOINorm;
    inertial_properties_offset.CMVector = inertial_properties_reference.CM.vector - inertial_properties_noisy.CM.vector;
end

function visualize_distribution(visualize_input)
    visualize_flag = visualize_input.visualize_flag;
    
    if visualize_flag
        N = visualize_input.N;
        mass_distribution_reference = visualize_input.mass_distribution_reference;
        mass_distribution_noisy = visualize_input.mass_distribution_noisy;
        mass_distribution_optimized = visualize_input.output_eps_star.mass_distribution;
        inertial_properties_reference = visualize_input.inertial_properties_reference;
        inertial_properties_noisy = visualize_input.inertial_properties_noisy;
        inertial_properties_optimized = visualize_input.output_eps_star.inertial_properties;
        cylinder_data = visualize_input.cylinder_data;

        X = cylinder_data.X;
        Y = cylinder_data.Y;
        Z = cylinder_data.Z;

        color_reference_list = ["r", "g", "b"];
        color_noisy_list = ["r--", "g--", "b--"];
        color_optimized_list = ["r:", "g:", "b:"];

        figure('units','normalized','outerposition',[0 0 1 1]);
        surf(X,Y,Z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', '0.1');
        hold on;
        for k=1:N
            scatter3(mass_distribution_reference(k).X, mass_distribution_reference(k).Y, mass_distribution_reference(k).Z, 'm', 'filled');
            scatter3(mass_distribution_noisy(k).X, mass_distribution_noisy(k).Y, mass_distribution_noisy(k).Z, 'm');
            scatter3(mass_distribution_optimized(k).X, mass_distribution_optimized(k).Y, mass_distribution_optimized(k).Z, 'cyan', 'LineWidth', mass_distribution_optimized(k).mass);
        end
        X_list.reference = nan(N+1,1);
        Y_list.reference = nan(N+1,1);
        Z_list.reference = nan(N+1,1);

        X_list.noisy = nan(N+1,1);
        Y_list.noisy = nan(N+1,1);
        Z_list.noisy = nan(N+1,1);

        for k=1:N
            X_list.reference(k) = mass_distribution_reference(k).X;
            Y_list.reference(k) = mass_distribution_reference(k).Y;
            Z_list.reference(k) = mass_distribution_reference(k).Z;

            X_list.noisy(k) = mass_distribution_noisy(k).X;
            Y_list.noisy(k) = mass_distribution_noisy(k).Y;
            Z_list.noisy(k) = mass_distribution_noisy(k).Z;
        end

        X_list.reference(N+1) = mass_distribution_reference(1).X;
        Y_list.reference(N+1) = mass_distribution_reference(1).Y;
        Z_list.reference(N+1) = mass_distribution_reference(1).Z;

        X_list.noisy(N+1) = mass_distribution_noisy(1).X;
        Y_list.noisy(N+1) = mass_distribution_noisy(1).Y;
        Z_list.noisy(N+1) = mass_distribution_noisy(1).Z;

        if N==3
            plot3(X_list.reference, Y_list.reference, Z_list.reference, 'black')
            plot3(X_list.noisy, Y_list.noisy, Z_list.noisy, 'black--')
        end
        scatter3(inertial_properties_reference.CM.x, inertial_properties_reference.CM.y, inertial_properties_reference.CM.z, 'b', 'filled');
        scatter3(inertial_properties_noisy.CM.x, inertial_properties_noisy.CM.y, inertial_properties_noisy.CM.z, 'b');
        scatter3(inertial_properties_optimized.CM.x, inertial_properties_optimized.CM.y, inertial_properties_optimized.CM.z, 'cyan');

        for j=1:3
            current_axis_vector = struct;

            current_axis_vector.reference.start = inertial_properties_reference.CM.vector;
            current_axis_vector.reference.end = inertial_properties_reference.CM.vector + inertial_properties_reference.PrincipalAxis(j).vector;
            current_axis_vector.reference.x = [current_axis_vector.reference.start(1), current_axis_vector.reference.end(1)];
            current_axis_vector.reference.y = [current_axis_vector.reference.start(2), current_axis_vector.reference.end(2)];
            current_axis_vector.reference.z = [current_axis_vector.reference.start(3), current_axis_vector.reference.end(3)];

            current_axis_vector.noisy.start = inertial_properties_noisy.CM.vector;
            current_axis_vector.noisy.end = inertial_properties_noisy.CM.vector + inertial_properties_noisy.PrincipalAxis(j).vector;
            current_axis_vector.noisy.x = [current_axis_vector.noisy.start(1), current_axis_vector.noisy.end(1)];
            current_axis_vector.noisy.y = [current_axis_vector.noisy.start(2), current_axis_vector.noisy.end(2)];
            current_axis_vector.noisy.z = [current_axis_vector.noisy.start(3), current_axis_vector.noisy.end(3)];

            current_axis_vector.optimized.start = inertial_properties_optimized.CM.vector;
            current_axis_vector.optimized.end = inertial_properties_optimized.CM.vector + inertial_properties_optimized.PrincipalAxis(j).vector;
            current_axis_vector.optimized.x = [current_axis_vector.optimized.start(1), current_axis_vector.optimized.end(1)];
            current_axis_vector.optimized.y = [current_axis_vector.optimized.start(2), current_axis_vector.optimized.end(2)];
            current_axis_vector.optimized.z = [current_axis_vector.optimized.start(3), current_axis_vector.optimized.end(3)];

            plot3(current_axis_vector.reference.x, current_axis_vector.reference.y, current_axis_vector.reference.z, color_reference_list(j), 'LineWidth', 1);
            plot3(current_axis_vector.noisy.x, current_axis_vector.noisy.y, current_axis_vector.noisy.z, color_noisy_list(j), 'LineWidth', 1);
            plot3(current_axis_vector.optimized.x, current_axis_vector.optimized.y, current_axis_vector.optimized.z, color_optimized_list(j), 'LineWidth', 1);

        end
        description_axes = zeros(13, 1);
        description_axes(1) = plot(NaN,NaN,'m');
        description_axes(2) = plot(NaN,NaN,'-m');
        description_axes(3) = plot(NaN,NaN,'black');
        description_axes(4) = plot(NaN,NaN,'-black');
        description_axes(5) = plot(NaN,NaN,'b');
        description_axes(6) = plot(NaN,NaN,'-b');
        description_axes(7) = plot(NaN,NaN,'r');
        description_axes(8) = plot(NaN,NaN,'-r');
        description_axes(9) = plot(NaN,NaN,'g');
        description_axes(10) = plot(NaN,NaN,'-g');
        description_axes(11) = plot(NaN,NaN,'b');
        description_axes(12) = plot(NaN,NaN,'-b');
        description_axes(13) = plot(NaN,NaN,'-cyan');
        legend(description_axes, 'reference markers: magenta solid', 'noisy markers: magenta hollow', 'reference boundary: black solid', 'noisy boundary: black dashed', 'reference CM: blue solid', 'noisy CM: blue hollow', 'reference axis one: red solid', 'noisy axis one: red dashed', 'reference axis two: green solid', 'noisy axis two: green dashed', 'reference axis three: blue solid', 'noisy axis three: blue dashed', 'deformed cylinder surface: transparent blue');
        zlim([-1,51])
        fileName_img = "MarkerConfigurationImages\MarkerDistribution_" + string(visualize_input.TrialNumber) + ".png";
        saveas(gcf,fileName_img);
        close(gcf);
    end
end

function visualize_EpsilonTransform(visualize_input)
    visualize_flag = visualize_input.visualize_flag;
   
    inertial_properties_reference = visualize_input.inertial_properties_reference;
    inertial_properties_noisy = visualize_input.inertial_properties_noisy;
    output_eps_star = visualize_input.output_eps_star;
    
    if visualize_flag
        r = inertial_properties_reference.PrincipalMOINorm;
        [X, Y, Z] = sphere;
        X2 = X*r;
        Y2 = Y*r;
        Z2 = Z*r;
        figure;
        surf(X2,Y2,Z2, 'FaceColor', 'b', 'FaceAlpha', '0.1');
        hold on;
        plot3([0, inertial_properties_reference.PrincipalMOI(1).value], [0, inertial_properties_reference.PrincipalMOI(2).value], [0, inertial_properties_reference.PrincipalMOI(3).value], 'r');
        plot3([0, inertial_properties_noisy.PrincipalMOI(1).value], [0, inertial_properties_noisy.PrincipalMOI(2).value], [0, inertial_properties_noisy.PrincipalMOI(3).value], 'g');
        plot3([0, output_eps_star.inertial_properties.PrincipalMOI(1).value], [0, output_eps_star.inertial_properties.PrincipalMOI(2).value], [0, output_eps_star.inertial_properties.PrincipalMOI(3).value], 'b');
        scatter3(inertial_properties_reference.PrincipalMOI(1).value, inertial_properties_reference.PrincipalMOI(2).value, inertial_properties_reference.PrincipalMOI(3).value, 'r');
        scatter3(inertial_properties_noisy.PrincipalMOI(1).value, inertial_properties_noisy.PrincipalMOI(2).value, inertial_properties_noisy.PrincipalMOI(3).value, 'g');
        scatter3(output_eps_star.inertial_properties.PrincipalMOI(1).value, output_eps_star.inertial_properties.PrincipalMOI(2).value, output_eps_star.inertial_properties.PrincipalMOI(3).value, 'b');
    end

    vectorNoisy = [inertial_properties_noisy.PrincipalMOI(1).value; inertial_properties_noisy.PrincipalMOI(2).value; inertial_properties_noisy.PrincipalMOI(3).value];
    vectorEpsilon = [output_eps_star.inertial_properties.PrincipalMOI(1).value; output_eps_star.inertial_properties.PrincipalMOI(2).value; output_eps_star.inertial_properties.PrincipalMOI(3).value];
    vectorRotationAxis = cross(vectorNoisy, vectorEpsilon);
    vectorRotationAxisUnit = vectorRotationAxis/norm(vectorRotationAxis);
    disp("Components of equivalent rotation axis in eigenvector space")
    disp(vectorRotationAxisUnit);
end

function position_local = global2local(T, R, position_global)
    position_local = nan(3,1);
    tol = 1.0e-10;
    if ~any(any(isnan(R)))
        same_inverse_transpose_bool = all(all(abs(inv(R) - R') < tol));
        det_R_bool = abs(det(R)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            position_local = R'*(position_global - T);
        end
    end
end

function position_local = global2local_perturbed(T, R, position_global)
    position_local = nan(3,1);
    tol = 1.0e-3;
    if ~any(any(isnan(R)))
        same_inverse_transpose_bool = all(all(abs(inv(R) - R') < tol));
        det_R_bool = abs(det(R)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            position_local = R'*(position_global - T);
        end
    end
end

function position_global = local2global(T, R, position_local)
    position_global = nan(3,1);
    tol = 1.0e-10;
    if ~any(any(isnan(R)))
        same_inverse_transpose_bool = all(all(abs(inv(R) - R') < tol));
        det_R_bool = abs(det(R)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            position_global = T + R*position_local;
        end
    end
end

function position_global = local2global_perturbed(T, R, position_local)
    position_global = nan(3,1);
    tol = 1.0e-3;
    if ~any(any(isnan(R)))
        same_inverse_transpose_bool = all(all(abs(inv(R) - R') < tol));
        det_R_bool = abs(det(R)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            position_global = T + R*position_local;
        end
    end
end

function [vectors_final, lambda_final, inertiaTensor_final, iterate_data] = iterative_perturbation(vectors_start, lambda_start, lambda_target)

    s=1;
    s_min = 1/(2^13);
    tol = 10^-7;

    vectors_current = vectors_start;
    lambda_current = lambda_start;
    Dlambda = s*(lambda_target - lambda_current);
    K_total = (vectors_start*diag(lambda_current))/vectors_start;

    lambda_iterates = [];
    vectorOne_iterates = [];
    vectorTwo_iterates = [];
    vectorThree_iterates = [];
    s_iterates = [];

    while s>s_min
        disp(s);
        lambda_iterates = [lambda_iterates, lambda_current];
        vectorOne_iterates = [vectorOne_iterates, vectors_current(:,1)];
        vectorTwo_iterates = [vectorTwo_iterates, vectors_current(:,2)];
        vectorThree_iterates = [vectorThree_iterates, vectors_current(:,3)];
        s_iterates = [s_iterates, s];
        A = generate_current_perturbation_lhs(vectors_current);
        p = lsqminnorm(A,Dlambda);
        dK = [p(1) p(2) p(3); p(2) p(4) p(5); p(3) p(5) p(6)];
        [dLambda, dVec] = compute_perturbation(vectors_current, lambda_current, dK);
        geometric_conditions = evaluate_geometric_conditions(vectors_current, dVec);
        if all(all(geometric_conditions))
            lambda_current = lambda_current + dLambda;
            vectors_current = vectors_current + dVec;
            K_total = K_total + dK;
            vectors_current(:,1) = vectors_current(:,1)/norm(vectors_current(:,1));
            vectors_current(:,2) = vectors_current(:,2)/norm(vectors_current(:,2));
            vectors_current(:,3) = vectors_current(:,3)/norm(vectors_current(:,3));
            disp(lambda_current);
        else
            s=s/2;
            Dlambda = Dlambda/2;
        end

        if norm(lambda_current - lambda_target)<tol
            break
        end
    end
    lambda_iterates = [lambda_iterates, lambda_current];
    s_iterates = [s_iterates, s];
    vectorOne_iterates = [vectorOne_iterates, vectors_current(:,1)];
    vectorTwo_iterates = [vectorTwo_iterates, vectors_current(:,2)];
    vectorThree_iterates = [vectorThree_iterates, vectors_current(:,3)];

    vectors_final = vectors_current;
    final_orthogonality_bool = evaluate_final_orthogonality(vectors_final);

    inertiaTensor_final = K_total;

    if all(abs(lambda_target - lambda_current) < tol)
        lambda_final = lambda_current;
    else
        lambda_final = nan(3,1);
        vectors_final = nan(3,3);
        inertiaTensor_final = nan(3,3);
    end

    if ~final_orthogonality_bool
        lambda_final = nan(3,1);
        vectors_final = nan(3,3);
        inertiaTensor_final = nan(3,3);
    end
    
    iterate_data = struct;
    iterate_data.s_iterates = s_iterates;
    iterate_data.lambda_iterates = lambda_iterates;
    iterate_data.vectorOne_iterates = vectorOne_iterates;
    iterate_data.vectorTwo_iterates = vectorTwo_iterates;
    iterate_data.vectorThree_iterates = vectorThree_iterates;

end

function F = nLsysMD(m, mass_distribution, inertiaTensor_final, N)
    x_sum = 0;
    y_sum = 0;
    z_sum = 0;
    m_sum = 0;
    for i=1:N
        x_sum = x_sum + m(i)*mass_distribution(i).X;
        y_sum = y_sum + m(i)*mass_distribution(i).Y;
        z_sum = z_sum + m(i)*mass_distribution(i).Z;
        m_sum = m_sum + m(i);
    end
    x_cm = x_sum/m_sum;
    y_cm = y_sum/m_sum;
    z_cm = z_sum/m_sum;
    
    sumVal = zeros(6,1);
    for i=1:N
        x_rel = mass_distribution(i).X - x_cm;
        y_rel = mass_distribution(i).Y - y_cm;
        z_rel = mass_distribution(i).Z - z_cm;
        sumVal(1) = sumVal(1) + m(i)*(y_rel^2 + z_rel^2);
        sumVal(2) = sumVal(2) - m(i)*x_rel*y_rel;
        sumVal(3) = sumVal(3) - m(i)*x_rel*z_rel;
        sumVal(4) = sumVal(4) + m(i)*(x_rel^2 + z_rel^2);
        sumVal(5) = sumVal(5) - m(i)*y_rel*z_rel;
        sumVal(6) = sumVal(6) + m(i)*(x_rel^2 + y_rel^2);
    end
    
    F(1) = sumVal(1) - inertiaTensor_final(1,1);
    F(2) = sumVal(2) - inertiaTensor_final(1,2);
    F(3) = sumVal(3) - inertiaTensor_final(1,3);
    F(4) = sumVal(4) - inertiaTensor_final(2,2);
    F(5) = sumVal(5) - inertiaTensor_final(2,3);
    F(6) = sumVal(6) - inertiaTensor_final(3,3);
end

function G = nLsysMDsum(m, mass_distribution, inertiaTensor_final, N)
    x_sum = 0;
    y_sum = 0;
    z_sum = 0;
    m_sum = 0;
    for i=1:N
        x_sum = x_sum + m(i)*mass_distribution(i).X;
        y_sum = y_sum + m(i)*mass_distribution(i).Y;
        z_sum = z_sum + m(i)*mass_distribution(i).Z;
        m_sum = m_sum + m(i);
    end
    x_cm = x_sum/m_sum;
    y_cm = y_sum/m_sum;
    z_cm = z_sum/m_sum;
    
    sumVal = zeros(6,1);
    for i=1:N
        x_rel = mass_distribution(i).X - x_cm;
        y_rel = mass_distribution(i).Y - y_cm;
        z_rel = mass_distribution(i).Z - z_cm;
        sumVal(1) = sumVal(1) + m(i)*(y_rel^2 + z_rel^2);
        sumVal(2) = sumVal(2) - m(i)*x_rel*y_rel;
        sumVal(3) = sumVal(3) - m(i)*x_rel*z_rel;
        sumVal(4) = sumVal(4) + m(i)*(x_rel^2 + z_rel^2);
        sumVal(5) = sumVal(5) - m(i)*y_rel*z_rel;
        sumVal(6) = sumVal(6) + m(i)*(x_rel^2 + y_rel^2);
    end
    
    F = nan(6,1);
    F(1) = sumVal(1) - inertiaTensor_final(1,1);
    F(2) = sumVal(2) - inertiaTensor_final(1,2);
    F(3) = sumVal(3) - inertiaTensor_final(1,3);
    F(4) = sumVal(4) - inertiaTensor_final(2,2);
    F(5) = sumVal(5) - inertiaTensor_final(2,3);
    F(6) = sumVal(6) - inertiaTensor_final(3,3);
    G = norm(F)^2;
end

function [c, ceq] = cm_constraint(m, mass_distribution, inertial_properties, N)
    x_sum = 0;
    y_sum = 0;
    z_sum = 0;
    m_sum = 0;
    cm_dist_threshold = 4;
    for i=1:N
        x_sum = x_sum + m(i)*mass_distribution(i).X;
        y_sum = y_sum + m(i)*mass_distribution(i).Y;
        z_sum = z_sum + m(i)*mass_distribution(i).Z;
        m_sum = m_sum + m(i);
    end
    x_cm = x_sum/m_sum;
    y_cm = y_sum/m_sum;
    z_cm = z_sum/m_sum;
    vector_cm = [x_cm; y_cm; z_cm];
    c = norm(inertial_properties.CM.vector - vector_cm) - cm_dist_threshold;
    ceq = [];
end

function vector_cm_reflected = cm_reflect(m, mass_distribution, inertial_properties, N)
    x_sum = 0;
    y_sum = 0;
    z_sum = 0;
    m_sum = 0;
    for i=1:N
        x_sum = x_sum + m(i)*mass_distribution(i).X;
        y_sum = y_sum + m(i)*mass_distribution(i).Y;
        z_sum = z_sum + m(i)*mass_distribution(i).Z;
        m_sum = m_sum + m(i);
    end
    x_cm = x_sum/m_sum;
    y_cm = y_sum/m_sum;
    z_cm = z_sum/m_sum;
    vector_cm = [x_cm; y_cm; z_cm];
    vector_cm_relative = vector_cm - inertial_properties.CM.vector;
    vector_cm_reflected = inertial_properties.CM.vector - vector_cm_relative;
end

function PCT_perturbation_performance = evaluate_PCT_perturbation(inertial_properties_perturbed_offset, inertial_properties_offset)
    PCT_perturbation_performance = struct;
    PCT_perturbation_performance.MOI_norm = abs(inertial_properties_perturbed_offset.PrincipalMOINorm) < abs(inertial_properties_offset.PrincipalMOINorm);
    PCT_perturbation_performance.CM_distance.boolVal = norm(inertial_properties_perturbed_offset.CMVector) < norm(inertial_properties_offset.CMVector);
    PCT_perturbation_performance.CM_distance.oldNorm = norm(inertial_properties_offset.CMVector);
    PCT_perturbation_performance.CM_distance.newNorm = norm(inertial_properties_perturbed_offset.CMVector);
    angleList = false(3,1);
    angleOldList = nan(3,1);
    angleNewList = nan(3,1);
    angleDiffList = nan(3,1);
    for m=1:3
        angleOld = inertial_properties_offset.angle(m);
        angleNew = inertial_properties_perturbed_offset.angle(m);

        if angleOld > 90
            angleOld = 180 - inertial_properties_offset.angle(m);
        end

        if angleNew > 90
            angleNew = 180 - inertial_properties_perturbed_offset.angle(m);
        end

        angleOldList(m) = angleOld;
        angleNewList(m) = angleNew;
        angleDiffList(m) = angleOld - angleNew;
        if angleNew < angleOld
            angleList(m) = true;
        end
    end
    PCT_perturbation_performance.angle.boolVal = all(angleList);
    PCT_perturbation_performance.angle.oldList = angleOldList;
    PCT_perturbation_performance.angle.newList = angleNewList;
end

function PCT_hybrid_performance = evaluate_PCT_hybrid(inertial_properties_hybrid_offset, inertial_properties_offset)
    PCT_hybrid_performance = struct;
    PCT_hybrid_performance.MOI_norm = abs(inertial_properties_hybrid_offset.PrincipalMOINorm) < abs(inertial_properties_offset.PrincipalMOINorm);
    PCT_hybrid_performance.CM_distance.boolVal = norm(inertial_properties_hybrid_offset.CMVector) < norm(inertial_properties_offset.CMVector);
    PCT_hybrid_performance.CM_distance.oldNorm = norm(inertial_properties_offset.CMVector);
    PCT_hybrid_performance.CM_distance.newNorm = norm(inertial_properties_hybrid_offset.CMVector);
    angleList = false(3,1);
    angleOldList = nan(3,1);
    angleNewList = nan(3,1);
    angleDiffList = nan(3,1);
    for m=1:3
        angleOld = inertial_properties_offset.angle(m);
        angleNew = inertial_properties_hybrid_offset.angle(m);

        if angleOld > 90
            angleOld = 180 - inertial_properties_offset.angle(m);
        end

        if angleNew > 90
            angleNew = 180 - inertial_properties_hybrid_offset.angle(m);
        end

        angleOldList(m) = angleOld;
        angleNewList(m) = angleNew;
        angleDiffList(m) = angleOld - angleNew;
        if angleNew < angleOld
            angleList(m) = true;
        end
    end
    PCT_hybrid_performance.angle.boolVal = all(angleList);
    PCT_hybrid_performance.angle.oldList = angleOldList;
    PCT_hybrid_performance.angle.newList = angleNewList;
end

function visualize_input = compute_Final_localOffset_perturbed(visualize_input)
    T = visualize_input.inertial_properties_perturbed.T;
    R = visualize_input.inertial_properties_perturbed.R;
    N = visualize_input.N;
    for k=1:N
        current_global_vector = visualize_input.mass_distribution_perturbed(k).vector;
        current_ref_vector = visualize_input.mass_distribution_reference(k).local_displacement;
        original_global_vector = visualize_input.mass_distribution_reference(k).vector;
        current_local_vector = global2local_perturbed(T, R, current_global_vector);
        global_vector_reconstructed = local2global_perturbed(T, R, current_ref_vector);
        global_offset = norm(global_vector_reconstructed - original_global_vector);
        visualize_input.mass_distribution_perturbed(k).local_displacement = current_local_vector;
        visualize_input.mass_distribution_perturbed(k).local_offset = norm(current_local_vector - current_ref_vector);
        visualize_input.mass_distribution_perturbed(k).global_reconstructed = global_vector_reconstructed;
        visualize_input.mass_distribution_perturbed(k).global_offset = global_offset;
    end
    
    inertial_properties_noisy_reversed = axis_reversal_original(visualize_input.inertial_properties_noisy, visualize_input.inertial_properties_offset);
    visualize_input.inertial_properties_noisy = inertial_properties_noisy_reversed;
    visualize_input.inertial_properties_offset = calculate_offset(visualize_input.inertial_properties_reference, visualize_input.inertial_properties_noisy);
    T_noisy = visualize_input.inertial_properties_noisy.T;
    R_noisy = visualize_input.inertial_properties_noisy.R;
    for k=1:N
        current_ref_vector = visualize_input.mass_distribution_reference(k).local_displacement;
        original_global_vector = visualize_input.mass_distribution_reference(k).vector;
        global_vector_reconstructed_noisy = local2global(T_noisy, R_noisy, current_ref_vector);
        global_offset_noisy = norm(global_vector_reconstructed_noisy - original_global_vector);
        visualize_input.mass_distribution_noisy(k).global_reconstructed = global_vector_reconstructed_noisy;
        visualize_input.mass_distribution_noisy(k).global_offset = global_offset_noisy;
    end
end

function visualize_input = compute_Final_localOffset_hybrid(visualize_input)
    T = visualize_input.inertial_properties_hybrid.T;
    R = visualize_input.inertial_properties_hybrid.R;
    N = visualize_input.N;
    for k=1:N
        current_global_vector = visualize_input.mass_distribution_hybrid(k).vector;
        current_ref_vector = visualize_input.mass_distribution_reference(k).local_displacement;
        original_global_vector = visualize_input.mass_distribution_reference(k).vector;
        current_local_vector = global2local(T, R, current_global_vector);
        global_vector_reconstructed = local2global(T, R, current_ref_vector);
        global_offset = norm(global_vector_reconstructed - original_global_vector);
        visualize_input.mass_distribution_hybrid(k).local_displacement = current_local_vector;
        visualize_input.mass_distribution_hybrid(k).local_offset = norm(current_local_vector - current_ref_vector);
        visualize_input.mass_distribution_hybrid(k).global_reconstructed = global_vector_reconstructed;
        visualize_input.mass_distribution_hybrid(k).global_offset = global_offset;
    end
    
    inertial_properties_noisy_reversed = axis_reversal_original(visualize_input.inertial_properties_noisy, visualize_input.inertial_properties_offset);
    visualize_input.inertial_properties_noisy = inertial_properties_noisy_reversed;
    visualize_input.inertial_properties_offset = calculate_offset(visualize_input.inertial_properties_reference, visualize_input.inertial_properties_noisy);
    T_noisy = visualize_input.inertial_properties_noisy.T;
    R_noisy = visualize_input.inertial_properties_noisy.R;
    for k=1:N
        current_ref_vector = visualize_input.mass_distribution_reference(k).local_displacement;
        original_global_vector = visualize_input.mass_distribution_reference(k).vector;
        global_vector_reconstructed_noisy = local2global(T_noisy, R_noisy, current_ref_vector);
        global_offset_noisy = norm(global_vector_reconstructed_noisy - original_global_vector);
        visualize_input.mass_distribution_noisy(k).global_reconstructed = global_vector_reconstructed_noisy;
        visualize_input.mass_distribution_noisy(k).global_offset = global_offset_noisy;
    end
end

function A = generate_current_perturbation_lhs(vectors_current)
    N = 3;
    A = nan(3,6);
    for j=1:N
        Vj = vectors_current(:,j);
        Vx = Vj(1);
        Vy = Vj(2);
        Vz = Vj(3);
        A(j,1) = Vx^2;
        A(j,2) = 2*Vx*Vy;
        A(j,3) = 2*Vx*Vz;
        A(j,4) = Vy^2;
        A(j,5) = 2*Vy*Vz;
        A(j,6) = Vz^2;
    end
end

function [dLambda, dVec] = compute_perturbation(vectors_current, lambda_current, dK)
    N=3;
    V = vectors_current;
    dLambda = nan(N,1);
    dVec = nan(N,N);
    D_perturbed_calc = zeros(3,3);
    for i=1:N
        dLambda(i) = transpose(V(:,i))*dK*V(:,i);
        dVecSum = zeros(N,1);
        for j=1:N
            if j~=i
                dVecSum = dVecSum + ((transpose(V(:,j))*dK*V(:,i))/(lambda_current(i)-lambda_current(j)))*V(:,j);
            end
        end
        dVec(:,i) = dVecSum;
    end
end

function geometric_conditions = evaluate_geometric_conditions(vectors_current, dVec)
    N=3;
    tol = 10^-7;
    norm_tol = 10^-4;
    vectors_perturbed = vectors_current + dVec;
    geometric_conditions = false(2,3);
    for j=1:N
        if abs(dot(dVec(:,j), vectors_current(:,j))) < tol
            geometric_conditions(1,j) = true;
        end
        if abs(norm(vectors_perturbed(:,j))-1) < norm_tol
            geometric_conditions(2,j) = true;
        end
    end
    
end

function final_orthogonality_check = evaluate_final_orthogonality(vectors_final)
    final_orthogonality_check = false;
    final_orthogonality_array = false(3,1);
    
    tol = 10^-3;
    
    v1 = vectors_final(:,1);
    v2 = vectors_final(:,2);
    v3 = vectors_final(:,3);
    
    if abs(dot(v1,v2)) < tol
        final_orthogonality_array(1) = true;
    end
    
    if abs(dot(v2,v3)) < tol
        final_orthogonality_array(2) = true;
    end
    
    if abs(dot(v3,v1)) < tol
        final_orthogonality_array(3) = true;
    end
    
    if all(final_orthogonality_array)
        final_orthogonality_check = true;
    end
end

function anatomical_landmarks_thigh_data = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose)
    anatomical_landmarks_thigh_data = struct;
    
    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;
    n = cylinder_data.n;
    np = cylinder_data.np;
    
    T_ref = reference_pose.T;
    R_ref = reference_pose.R;
    
    T_opt = optimized_pose.T;
    R_opt = optimized_pose.R;
    
    name_list = ["Lateral Epicondyle";
                 "Medial Epicondyle";
                 "Greater Trochanter";
                 "Femoral Head"];
             
    height_idx_list = [1, 1, n, n];
    plane_idx_list = [floor(np/4), floor(3*np/4), floor(np/4), 1];
    
    for j=1:4
        current_height_idx = height_idx_list(j);
        current_plane_idx = plane_idx_list(j);
        
        anatomical_landmarks_thigh_data(j).name = name_list(j);
        if j==4
            anatomical_landmarks_thigh_data(j).x_global_reference = 0;
        else
            anatomical_landmarks_thigh_data(j).x_global_reference = X(current_height_idx, current_plane_idx);
        end
        anatomical_landmarks_thigh_data(j).y_global_reference = Y(current_height_idx, current_plane_idx);
        anatomical_landmarks_thigh_data(j).z_global_reference = Z(current_height_idx, current_plane_idx);
        anatomical_landmarks_thigh_data(j).vector_global_reference = [anatomical_landmarks_thigh_data(j).x_global_reference;
                                                                      anatomical_landmarks_thigh_data(j).y_global_reference;
                                                                      anatomical_landmarks_thigh_data(j).z_global_reference];
        anatomical_landmarks_thigh_data(j).vector_local_reference = global2local_perturbed(T_ref, R_ref, anatomical_landmarks_thigh_data(j).vector_global_reference);
        anatomical_landmarks_thigh_data(j).vector_global_reconstructed = local2global_perturbed(T_opt, R_opt, anatomical_landmarks_thigh_data(j).vector_local_reference);
        anatomical_landmarks_thigh_data(j).global_offset = norm(anatomical_landmarks_thigh_data(j).vector_global_reconstructed - anatomical_landmarks_thigh_data(j).vector_global_reference);
    end
    
end