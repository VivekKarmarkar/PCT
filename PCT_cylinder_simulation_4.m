N = 4;
NumConfigurations = 1000;

GTD_filename = "ground_truth_data_thigh.mat";
load(GTD_filename);

STA_filename = "box_markers_gait_" + string(N) + ".mat";
load(STA_filename);
input_size = generate_input_size_struct(N, NumConfigurations, box_markers_gait);

input_settings = generate_input_settings;

cylinder_data = generate_cylinder_data;

final_output = struct;

% animation_flag = false(NumConfigurations, 1);
% animation_flag(1) = true;
% fig_scale_flag = false(NumConfigurations, 1);
% fig_scale_flag(1) = true;

visualization_flags = struct;
% visualization_flags.animation = false(NumConfigurations, 1);
% visualization_flags.fig_scale = false(NumConfigurations, 1);
% visualization_flags.animation(1) = true;
% visualization_flags.animation(15) = true;
% visualization_flags.fig_scale(1) = true;

visualization_flags.animation = true(NumConfigurations, 1);
visualization_flags.fig_scale = true(NumConfigurations, 1);

for j=1:NumConfigurations  
    configuration_idx = j;
    folderName_settings = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\MarkerConfigurationConstrainedSettings_" + string(N);
    fileName_configuration = folderName_settings + "\FixedMarkerDistribution_" + string(configuration_idx);
    load(fileName_configuration);

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

    output_containers = generate_output_containers(input_size);
    output_containers = process_all_algorithms(input_struct, output_containers);
    output_containers = post_processing_discontinuity(output_containers);
    
    final_output(configuration_idx).output_containers = output_containers;
    fig_scale = generate_fig_scale_info_PCT(visualization_flags.fig_scale(configuration_idx), input_struct, output_containers);
    input_struct.visualization_settings = struct;
    input_struct.visualization_settings.fig_scale = fig_scale;
    input_struct.visualization_settings.animation_flag = visualization_flags.animation(configuration_idx);
    input_struct.visualization_settings.movingView_flag_val = false;
    create_animation_RF(input_struct, output_containers, "TF");
    create_animation_RF(input_struct, output_containers, "AF");
    generate_output_plots(input_struct, output_containers);
end

final_output_mean = generate_ouput_stats(final_output, input_size);
save('STA_output', 'final_output', '-v7.3');

function cylinder_data_pose = generate_cylinder_data_pose(cylinder_data, current_pose)
    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;
    
    n = cylinder_data.n;
    np = cylinder_data.np;
    
    T = current_pose.T;
    R = current_pose.R;
    
    X_pose = nan(n, np);
    Y_pose = nan(n, np);
    Z_pose = nan(n, np);
    
    for k=1:n
        for j=1:np
            vector_standard = [X(k,j); Y(k,j); Z(k,j)];
            vector_STA = rotx(-90)*vector_standard;
            vector_pose = local2global_perturbed(T, R, vector_STA); 
            X_pose(k,j) = vector_pose(1);
            Y_pose(k,j) = vector_pose(2);
            Z_pose(k,j) = vector_pose(3);
        end
    end
    
    cylinder_data_pose = struct;
    cylinder_data_pose.X = X_pose;
    cylinder_data_pose.Y = Y_pose;
    cylinder_data_pose.Z = Z_pose;
    cylinder_data_pose.np = np;
    cylinder_data_pose.n = n;
end

function mass_distribution_pose = apply_current_pose(mass_distribution_standard_pose, current_pose, N)
    mass_distribution_pose = struct;
    
    T = current_pose.T;
    R = current_pose.R;
    
    for j=1:N
        vector_standard = mass_distribution_standard_pose(j).vector;
        vector_STA = rotx(-90)*vector_standard;
        vector_pose = T + R*vector_STA;
        mass_distribution_pose(j).mass = 1;
        mass_distribution_pose(j).X = vector_pose(1);
        mass_distribution_pose(j).Y = vector_pose(2);
        mass_distribution_pose(j).Z = vector_pose(3);
        mass_distribution_pose(j).vector = vector_pose;
    end
    
end

function final_output_mean = generate_ouput_stats(final_output, input_size)

    final_output_mean = struct;
    
    final_output_mean.T_estimated = struct;
    final_output_mean.T_estimated.SVDLS = nan(input_size.NumFrames, 3, 1);
    final_output_mean.T_estimated.PTUR = nan(input_size.NumFrames, 3, 1);

    final_output_mean.CM_offset = struct;
    final_output_mean.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_mean.CM_offset.PTUR = nan(input_size.NumFrames, 1);

    final_output_mean.cluster_offset = struct;
    final_output_mean.cluster_offset.SVDLS = nan(input_size.NumFrames, input_size.N);
    final_output_mean.cluster_offset.PTUR = nan(input_size.NumFrames, input_size.N);

    final_output_mean.AL_offset = struct;
    final_output_mean.AL_offset.SVDLS = nan(input_size.NumFrames, input_size.Num_AL_Thigh);
    final_output_mean.AL_offset.PTUR = nan(input_size.NumFrames, input_size.Num_AL_Thigh);

    NumEigenvectors = 3;
    final_output_mean.eigenvector_offset = struct;
    final_output_mean.eigenvector_offset.SVDLS = nan(input_size.NumFrames, NumEigenvectors);
    final_output_mean.eigenvector_offset.PTUR = nan(input_size.NumFrames, NumEigenvectors);
    
    NumConfigurations = min([length(final_output)-1, input_size.NumConfigurations]);
    
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    ToeOff = 0.6;

    for t=1:input_size.NumFrames
        disp(t);
        
        T_estimated_current = struct;
        T_estimated_current.SVDLS = nan(NumConfigurations, 3, 1);
        T_estimated_current.PTUR = nan(NumConfigurations, 3, 1);
        
        CM_offset_current = struct;
        CM_offset_current.SVDLS = nan(NumConfigurations, 1);
        CM_offset_current.PTUR = nan(NumConfigurations, 1);

        cluster_offset_current = struct;
        cluster_offset_current.SVDLS = nan(NumConfigurations, input_size.N);
        cluster_offset_current.PTUR = nan(NumConfigurations, input_size.N);

        AL_offset_current = struct;
        AL_offset_current.SVDLS = nan(NumConfigurations, input_size.Num_AL_Thigh);
        AL_offset_current.PTUR = nan(NumConfigurations, input_size.Num_AL_Thigh);

        eigenvector_offset_current = struct;
        eigenvector_offset_current.SVDLS = nan(NumConfigurations, NumEigenvectors);
        eigenvector_offset_current.PTUR = nan(NumConfigurations, NumEigenvectors);

        for k=1:NumConfigurations
            disp(k);
            T_estimated_current.SVDLS(k,:,1) = final_output(k).output_containers.T_estimated.SVDLS(t);
            T_estimated_current.PTUR(k,:,1) = final_output(k).output_containers.T_estimated.PTUR(t);
            CM_offset_current.SVDLS(k) = final_output(k).output_containers.CM_offset.SVDLS(t);
            CM_offset_current.PTUR(k) = final_output(k).output_containers.CM_offset.PTUR(t);
            for n=1:input_size.N
                cluster_offset_current.SVDLS(k,n) = final_output(k).output_containers.cluster_offset.SVDLS(t,n);
                cluster_offset_current.PTUR(k,n) = final_output(k).output_containers.cluster_offset.PTUR(t,n);
            end

            for n=1:input_size.Num_AL_Thigh
                AL_offset_current.SVDLS(k,n) = final_output(k).output_containers.AL_offset.SVDLS(t,n);
                AL_offset_current.PTUR(k,n) = final_output(k).output_containers.AL_offset.PTUR(t,n);
            end

            for n=1:NumEigenvectors
                eigenvector_offset_current.SVDLS(k,n) = final_output(k).output_containers.eigenvector_offset.SVDLS(t,n);
                eigenvector_offset_current.PTUR(k,n) = final_output(k).output_containers.eigenvector_offset.PTUR(t,n);
            end

        end
        final_output_mean.CM_offset.SVDLS(t) = nanmean(CM_offset_current.SVDLS);
        final_output_mean.CM_offset.PTUR(t) = nanmean(CM_offset_current.PTUR);

        for n=1:input_size.N
            final_output_mean.cluster_offset.SVDLS(t,n) = nanmean(cluster_offset_current.SVDLS(:,n));
            final_output_mean.cluster_offset.PTUR(t,n) = nanmean(cluster_offset_current.PTUR(:,n));
        end

        for n=1:input_size.Num_AL_Thigh
            final_output_mean.AL_offset.SVDLS(t,n) = nanmean(AL_offset_current.SVDLS(:,n));
            final_output_mean.AL_offset.PTUR(t,n) = nanmean(AL_offset_current.PTUR(:,n));
        end

        for n=1:NumEigenvectors
            final_output_mean.eigenvector_offset.SVDLS(t,n) = nanmean(eigenvector_offset_current.SVDLS(:,n));
            final_output_mean.eigenvector_offset.PTUR(t,n) = nanmean(eigenvector_offset_current.PTUR(:,n));
        end
    end

    figure
    xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    plot(fraction_completed_gait_cycle, final_output_mean.CM_offset.SVDLS, 'k', 'DisplayName', 'SVDLS');
    plot(fraction_completed_gait_cycle, final_output_mean.CM_offset.PTUR, 'b', 'DisplayName', 'PTUR');
    grid;
    legend;
    title('Average CM offset during gait cycle');
    xlabel('Frame Number');
    ylabel('Error (cm)');

    figure
    for n=1:input_size.N
        subplot(2, input_size.N/2, n)
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_mean.cluster_offset.SVDLS(:,n), 'k', 'DisplayName', 'SVDLS');
        plot(fraction_completed_gait_cycle, final_output_mean.cluster_offset.PTUR(:,n), 'b', 'DisplayName', 'PTUR');
        grid;
        if n==3
            legend;
        end
        if n>3
            xlabel("Fraction completed gait cycle");
        end
        if or(n==1, n==4)
            ylabel("Error (cm)");
        end
        title("Avg Error Cluster Point IDX = " + string(n));
    end

    AL_name_list = ["Lateral Epicondyle";
                     "Medial Epicondyle";
                     "Greater Trochanter";
                     "Femoral Head"];

    figure
    for n=1:input_size.Num_AL_Thigh
        subplot(2, input_size.Num_AL_Thigh/2, n)
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_mean.AL_offset.SVDLS(:,n), 'k', 'DisplayName', 'SVDLS');
        plot(fraction_completed_gait_cycle, final_output_mean.AL_offset.PTUR(:,n), 'b', 'DisplayName', 'PTUR');
        grid;
        if n==2
            legend;
        end
        if n>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(n,2) == 0)
            ylabel("Error (cm)");
        end
        title("Avg Error for " + AL_name_list(n));
    end

    figure
    for n=1:NumEigenvectors
        subplot(NumEigenvectors, 1, n)
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_mean.eigenvector_offset.SVDLS(:,n), 'k', 'DisplayName', 'SVDLS');
        plot(fraction_completed_gait_cycle, final_output_mean.eigenvector_offset.PTUR(:,n), 'b', 'DisplayName', 'PTUR');
        grid;
        if n==1
            legend;
        end
        if n==3
            xlabel("Fraction completed gait cycle");
        end
        ylabel('Error (degress)');
        title("Avg Error for Eigenvector IDX = " + string(n));
    end

end

function output_containers = post_processing_discontinuity(output_containers)
    roots_info = find_roots_location(output_containers.deltaY_offset.PTNR);
    roots_idx = roots_info.roots_idx;
    roots_num = length(roots_idx);
    if roots_num > 0
        first_root = roots_idx(1);
        first_root_switch_negative_bool = roots_info.switch_negative_idx(first_root);
        if mod(roots_num,2) == 0
            for j=1:roots_num-1
                current_root = roots_idx(j);
                if first_root_switch_negative_bool
                    if roots_info.switch_negative_idx(current_root)
                        output_containers.CM_offset.PTUR(roots_idx(j): roots_idx(j+1)) = output_containers.CM_offset.PTLR(roots_idx(j): roots_idx(j+1));
                        output_containers.cluster_offset.PTUR(roots_idx(j): roots_idx(j+1), :) = output_containers.cluster_offset.PTLR(roots_idx(j): roots_idx(j+1), :);
                        output_containers.AL_offset.PTUR(roots_idx(j): roots_idx(j+1), :) = output_containers.AL_offset.PTLR(roots_idx(j): roots_idx(j+1), :);
                    end
                else
                    if roots_info.switch_positive_idx(current_root)
                        output_containers.CM_offset.PTUR(roots_idx(j): roots_idx(j+1)) = output_containers.CM_offset.PTNR(roots_idx(j): roots_idx(j+1));
                        output_containers.cluster_offset.PTUR(roots_idx(j): roots_idx(j+1), :) = output_containers.cluster_offset.PTNR(roots_idx(j): roots_idx(j+1), :);
                        output_containers.AL_offset.PTUR(roots_idx(j): roots_idx(j+1), :) = output_containers.AL_offset.PTNR(roots_idx(j): roots_idx(j+1), :);
                    end
                end
            end
        end
    end
    
    output_containers.eigenvalueNorm_offset.PTNR(1:2) = nan;
    output_containers.eigenvector_offset.PTNR(1:2, :) = nan;
    output_containers.mass_redistribution.PTNR(1:2, :) = nan;
    output_containers.CM_offset.PTNR(1:2) = nan;
    output_containers.cluster_offset.PTNR(1:2, :) = nan;
    output_containers.AL_offset.PTNR(1:2, :) = nan;
    output_containers.deltaY_offset.PTNR(1:2) = nan;
    output_containers.T_estimated.PTNR(1:2, :, 1) = nan(2,3);
    output_containers.R_estimated.PTNR(1, :, :) = nan(3,3);
    output_containers.R_estimated.PTNR(2, :, :) = nan(3,3);
    output_containers.AL_estimated.PTNR(1:2, :, 1) = nan(2,4);
    output_containers.AL_reference.PTNR(1:2, :, 1) = nan(2,4);
    output_containers.AF_Translation_estimated.PTNR(1:2, :, 1) = nan(2,3);
    output_containers.AF_Rotation_estimated.PTNR(1, :, :) = nan(3,3);
    output_containers.AF_Rotation_estimated.PTNR(2, :, :) = nan(3,3);
    
    output_containers.eigenvalueNorm_offset.PTUR(1:2) = nan;
    output_containers.eigenvector_offset.PTUR(1:2, :) = nan;
    output_containers.mass_redistribution.PTUR(1:2, :) = nan;
    output_containers.CM_offset.PTUR(1:2) = nan;
    output_containers.cluster_offset.PTUR(1:2, :) = nan;
    output_containers.AL_offset.PTUR(1:2, :) = nan;
    output_containers.deltaY_offset.PTUR(1:2) = nan;
    output_containers.T_estimated.PTUR(1:2, :, 1) = nan(2,3);
    output_containers.R_estimated.PTUR(1, :, :) = nan(3,3);
    output_containers.R_estimated.PTUR(2, :, :) = nan(3,3);
    output_containers.AL_estimated.PTUR(1:2, :, 1) = nan(2,4);
    output_containers.AL_reference.PTUR(1:2, :, 1) = nan(2,4);
    output_containers.AF_Translation_estimated.PTUR(1:2, :, 1) = nan(2,3);
    output_containers.AF_Rotation_estimated.PTUR(1, :, :) = nan(3,3);
    output_containers.AF_Rotation_estimated.PTUR(2, :, :) = nan(3,3);
    
    output_containers.eigenvalueNorm_offset.PTLR(1:2) = nan;
    output_containers.eigenvector_offset.PTLR(1:2, :) = nan;
    output_containers.mass_redistribution.PTLR(1:2, :) = nan;
    output_containers.CM_offset.PTLR(1:2) = nan;
    output_containers.cluster_offset.PTLR(1:2, :) = nan;
    output_containers.AL_offset.PTLR(1:2, :) = nan;
    output_containers.deltaY_offset.PTLR(1:2) = nan;
    output_containers.T_estimated.PTLR(1:2, :, 1) = nan(2,3);
    output_containers.R_estimated.PTLR(1, :, :) = nan(3,3);
    output_containers.R_estimated.PTLR(2, :, :) = nan(3,3);
    output_containers.AL_estimated.PTLR(1:2, :, 1) = nan(2,4);
    output_containers.AL_reference.PTLR(1:2, :, 1) = nan(2,4);
    output_containers.AF_Translation_estimated.PTLR(1:2, :, 1) = nan(2,3);
    output_containers.AF_Rotation_estimated.PTLR(1, :, :) = nan(3,3);
    output_containers.AF_Rotation_estimated.PTLR(2, :, :) = nan(3,3);
    
end

function roots_info = find_roots_location(data)
    NumPoints = length(data);
    roots_idx = [];
    switch_negative_idx = false(NumPoints, 1);
    switch_positive_idx = false(NumPoints, 1);
    for j=1:NumPoints-1
        if (data(j)>0) && (data(j+1)<0)
            roots_idx = [roots_idx; j+1];
            switch_negative_idx(j+1) = true;
        elseif (data(j)<0) && (data(j+1)>0)
            roots_idx = [roots_idx; j];
            switch_positive_idx(j) = true;
        end
    end
    
    roots_info = struct;
    roots_info.roots_idx = roots_idx;
    roots_info.switch_negative_idx = switch_negative_idx;
    roots_info.switch_positive_idx = switch_positive_idx;
end

function output_containers = process_all_algorithms(input_struct, output_containers)
    input_size = input_struct.input_size;
    input_settings = input_struct.input_settings;
    input_data = input_struct.input_data;
    
    N = input_size.N;
    NumFrames = input_size.NumFrames;
    Num_AL_Thigh = input_size.Num_AL_Thigh;
    NumAngles = input_size.NumAngles;
    
    initial_axes_reversal = input_settings.initial_axes_reversal;
    initial_axes_perturbed_reversal = input_settings.initial_axes_perturbed_reversal;
    optimized_axes_reversal = input_settings.optimized_axes_reversal;
    m0_switch_idx = input_settings.m0_switch_idx;
    
    cylinder_data = input_data.cylinder_data;
    mass_distribution_reference_standard_pose = input_data.mass_distribution_reference;
    box_markers_gait = input_data.box_markers_gait;
    pose_gait = input_data.pose_gait;

    for j=1:NumFrames
        STA_frame_idx = j;
        fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle(STA_frame_idx);
        cm_dist_threshold = compute_cm_dist_threshold_STA(fraction_completed_gait_cycle);

        mass_distribution_noisy_standard_pose = generate_noisy_distribution_STA(mass_distribution_reference_standard_pose, STA_frame_idx, box_markers_gait, N);
        
        current_pose = pose_gait(STA_frame_idx);
        mass_distribution_reference = apply_current_pose(mass_distribution_reference_standard_pose, current_pose, N);
        mass_distribution_noisy = apply_current_pose(mass_distribution_noisy_standard_pose, current_pose, N);
        
        [inertial_properties_reference, mass_distribution_reference] = compute_inertial_properties(mass_distribution_reference, N);
        [inertial_properties_noisy, mass_distribution_noisy] = compute_inertial_properties(mass_distribution_noisy, N);
        inertial_properties_offset = calculate_offset(inertial_properties_reference, inertial_properties_noisy);
        if initial_axes_reversal
            inertial_properties_noisy_reversed = axis_reversal_original(inertial_properties_noisy, inertial_properties_offset);
            inertial_properties_noisy = inertial_properties_noisy_reversed;
            inertial_properties_offset = calculate_offset(inertial_properties_reference, inertial_properties_noisy);
        end

        inertial_properties_noisy_copy = inertial_properties_noisy;
        inertial_properties_offset_copy = calculate_offset(inertial_properties_reference, inertial_properties_noisy_copy);
        disp("Original offset = " + string(inertial_properties_offset_copy.angle));

        if initial_axes_perturbed_reversal
            inertial_properties_noisy_reversed = axis_reversal_original(inertial_properties_noisy_copy, inertial_properties_offset_copy);
            inertial_properties_noisy_copy = inertial_properties_noisy_reversed;
            inertial_properties_offset_copy = calculate_offset(inertial_properties_reference, inertial_properties_noisy_copy);
        end
        disp("New offset = " + string(inertial_properties_offset_copy.angle));

        algorithm_container = struct;

        algorithm_container.PCTUO = struct;
        algorithm_container.PCTUO.TrialNumber = STA_frame_idx;
        algorithm_container.PCTUO.mass_distribution_reference = mass_distribution_reference;
        algorithm_container.PCTUO.mass_distribution_noisy = mass_distribution_noisy;
        algorithm_container.PCTUO.inertial_properties_reference = inertial_properties_reference;
        algorithm_container.PCTUO.inertial_properties_noisy = inertial_properties_noisy;
        algorithm_container.PCTUO.inertial_properties_offset = inertial_properties_offset;
        algorithm_container.PCTUO.N = N;
        algorithm_container.PCTUO.cylinder_data = cylinder_data;
        algorithm_container.PCTUO.current_pose = current_pose;
        algorithm_container.PCTUO.optimized_axis_reversal_bool = optimized_axes_reversal;

        [algorithm_container.PCTUO, local_offset, local_angular_offset, local_offset_vector] = compute_localOffset(algorithm_container.PCTUO);
        algorithm_container.PCTUO.local_offset = local_offset;
        algorithm_container.PCTUO.local_angular_offset = local_angular_offset;
        algorithm_container.PCTUO.local_offset_vector = local_offset_vector;

        algorithm_container.SVDLS = algorithm_container.PCTUO;
        algorithm_container.SVDLS.mass_distribution_SVDLS = mass_distribution_noisy;
        algorithm_container.SVDLS = SVDLS_processing(algorithm_container.SVDLS);

        algorithm_container.PTNR = algorithm_container.PCTUO;
        algorithm_container.PTNR.inertial_properties_noisy_original = inertial_properties_noisy;
        algorithm_container.PTNR.inertial_properties_noisy = inertial_properties_noisy_copy;
        algorithm_container.PTNR.inertial_properties_offset = inertial_properties_offset_copy;
        [algorithm_container.PTNR, local_offset_copy_NR, local_angular_offset_copy_NR, local_offset_vector_copy_NR] = compute_localOffset(algorithm_container.PTNR);
        algorithm_container.PTNR.local_offset = local_offset_copy_NR;
        algorithm_container.PTNR.local_angular_offset = local_angular_offset_copy_NR;
        algorithm_container.PTNR.local_offset_vector = local_offset_vector_copy_NR;
        algorithm_container.PTNR.cm_dist_threshold = cm_dist_threshold;
        algorithm_container.PTNR.fraction_completed_gait_cycle = fraction_completed_gait_cycle;

        if STA_frame_idx < m0_switch_idx
            m0 = ones(N,1);
        else
            m0 = output_containers.mass_redistribution.PTNR(STA_frame_idx - 1, :);
        end

        algorithm_container.PTNR.m0 = m0;
        algorithm_container.PTNR.reflection_setting = "NR";
        algorithm_container.PTNR = Perturbation_Theory_processing(algorithm_container.PTNR);
        
        algorithm_container.PTUR = algorithm_container.PCTUO;
        algorithm_container.PTUR.inertial_properties_noisy_original = inertial_properties_noisy;
        algorithm_container.PTUR.inertial_properties_noisy = inertial_properties_noisy_copy;
        algorithm_container.PTUR.inertial_properties_offset = inertial_properties_offset_copy;
        [algorithm_container.PTUR, local_offset_copy_UR, local_angular_offset_copy_UR, local_offset_vector_copy_UR] = compute_localOffset(algorithm_container.PTUR);
        algorithm_container.PTUR.local_offset = local_offset_copy_UR;
        algorithm_container.PTUR.local_angular_offset = local_angular_offset_copy_UR;
        algorithm_container.PTUR.local_offset_vector = local_offset_vector_copy_UR;
        algorithm_container.PTUR.cm_dist_threshold = cm_dist_threshold;
        algorithm_container.PTUR.fraction_completed_gait_cycle = fraction_completed_gait_cycle;

        if STA_frame_idx < m0_switch_idx
            m0 = ones(N,1);
        else
            m0 = output_containers.mass_redistribution.PTUR(STA_frame_idx - 1, :);
        end

        algorithm_container.PTUR.m0 = m0;
        algorithm_container.PTUR.reflection_setting = "UR";
        algorithm_container.PTUR = Perturbation_Theory_processing(algorithm_container.PTUR);
        
        algorithm_container.PTLR = algorithm_container.PCTUO;
        algorithm_container.PTLR.inertial_properties_noisy_original = inertial_properties_noisy;
        algorithm_container.PTLR.inertial_properties_noisy = inertial_properties_noisy_copy;
        algorithm_container.PTLR.inertial_properties_offset = inertial_properties_offset_copy;
        [algorithm_container.PTLR, local_offset_copy_LR, local_angular_offset_copy_LR, local_offset_vector_copy_LR] = compute_localOffset(algorithm_container.PTLR);
        algorithm_container.PTLR.local_offset = local_offset_copy_LR;
        algorithm_container.PTLR.local_angular_offset = local_angular_offset_copy_LR;
        algorithm_container.PTLR.local_offset_vector = local_offset_vector_copy_LR;
        algorithm_container.PTLR.cm_dist_threshold = cm_dist_threshold;
        algorithm_container.PTLR.fraction_completed_gait_cycle = fraction_completed_gait_cycle;

        if STA_frame_idx < m0_switch_idx
            m0 = ones(N,1);
        else
            m0 = output_containers.mass_redistribution.PTLR(STA_frame_idx - 1, :);
        end

        algorithm_container.PTLR.m0 = m0;
        algorithm_container.PTLR.reflection_setting = "LR";
        algorithm_container.PTLR = Perturbation_Theory_processing(algorithm_container.PTLR);

        for k=1:N
            output_containers.mass_redistribution.SVDLS(j, k, 1) = algorithm_container.SVDLS.mass_distribution_SVDLS(k).mass;
            output_containers.mass_redistribution.PTNR(j, k, 1) = algorithm_container.PTNR.mass_distribution_perturbed(k).mass;
            output_containers.mass_redistribution.PTUR(j, k, 1) = algorithm_container.PTUR.mass_distribution_perturbed(k).mass;
            output_containers.mass_redistribution.PTLR(j, k, 1) = algorithm_container.PTLR.mass_distribution_perturbed(k).mass;
            output_containers.cluster_offset.SVDLS(j, k, 1) = algorithm_container.SVDLS.mass_distribution_SVDLS(k).global_offset;
            output_containers.cluster_offset.PTNR(j, k, 1) = algorithm_container.PTNR.mass_distribution_perturbed(k).global_offset;
            output_containers.cluster_offset.PTUR(j, k, 1) = algorithm_container.PTUR.mass_distribution_perturbed(k).global_offset;
            output_containers.cluster_offset.PTLR(j, k, 1) = algorithm_container.PTLR.mass_distribution_perturbed(k).global_offset;
        end

        for k=1:Num_AL_Thigh
            output_containers.AL_offset.SVDLS(j, k, 1) = algorithm_container.SVDLS.anatomical_landmarks_thigh_data(k).global_offset;
            output_containers.AL_offset.PTNR(j, k, 1) = algorithm_container.PTNR.anatomical_landmarks_thigh_data(k).global_offset;
            output_containers.AL_offset.PTUR(j, k, 1) = algorithm_container.PTUR.anatomical_landmarks_thigh_data(k).global_offset;
            output_containers.AL_offset.PTLR(j, k, 1) = algorithm_container.PTLR.anatomical_landmarks_thigh_data(k).global_offset;
        end

        current_PCT_performanceSVDLS = evaluate_PCT_SVDLS(algorithm_container.SVDLS.inertial_properties_SVDLS_offset, inertial_properties_offset);
        current_PCT_performancePerturbedNR = evaluate_PCT_perturbation(algorithm_container.PTNR.inertial_properties_perturbed_offset, inertial_properties_offset_copy);
        current_PCT_performancePerturbedUR = evaluate_PCT_perturbation(algorithm_container.PTUR.inertial_properties_perturbed_offset, inertial_properties_offset_copy);
        current_PCT_performancePerturbedLR = evaluate_PCT_perturbation(algorithm_container.PTLR.inertial_properties_perturbed_offset, inertial_properties_offset_copy);
        for k=1:NumAngles
            output_containers.eigenvector_offset.SVDLS(j, k, 1) = current_PCT_performanceSVDLS.angle.newList(k);
            output_containers.eigenvector_offset.PTNR(j, k, 1) = current_PCT_performancePerturbedNR.angle.newList(k);
            output_containers.eigenvector_offset.PTUR(j, k, 1) = current_PCT_performancePerturbedUR.angle.newList(k);
            output_containers.eigenvector_offset.PTLR(j, k, 1) = current_PCT_performancePerturbedLR.angle.newList(k);
        end

        output_containers.CM_offset.SVDLS(j) = current_PCT_performanceSVDLS.CM_distance.newNorm;
        output_containers.CM_offset.PTNR(j) = current_PCT_performancePerturbedNR.CM_distance.newNorm;
        output_containers.CM_offset.PTUR(j) = current_PCT_performancePerturbedUR.CM_distance.newNorm;
        output_containers.CM_offset.PTLR(j) = current_PCT_performancePerturbedLR.CM_distance.newNorm;

        output_containers.eigenvalueNorm_offset.PTNR(j) = (abs(algorithm_container.PTNR.inertial_properties_perturbed.PrincipalMOINorm - inertial_properties_reference.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
        output_containers.eigenvalueNorm_offset.PTUR(j) = (abs(algorithm_container.PTUR.inertial_properties_perturbed.PrincipalMOINorm - inertial_properties_reference.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
        output_containers.eigenvalueNorm_offset.PTLR(j) = (abs(algorithm_container.PTLR.inertial_properties_perturbed.PrincipalMOINorm - inertial_properties_reference.PrincipalMOINorm)/inertial_properties_reference.PrincipalMOINorm)*100;
        
        output_containers.deltaY_offset.PTNR(j) = algorithm_container.PTNR.deltaY_offset;
        output_containers.deltaY_offset.PTUR(j) = algorithm_container.PTUR.deltaY_offset;
        output_containers.deltaY_offset.PTLR(j) = algorithm_container.PTLR.deltaY_offset;
        
        output_containers.T_estimated.SVDLS(j, :, 1) = algorithm_container.SVDLS.inertial_properties_SVDLS.T;
        output_containers.T_estimated.PTNR(j, :, 1) = algorithm_container.PTNR.inertial_properties_perturbed.T;
        output_containers.T_estimated.PTUR(j, :, 1) = algorithm_container.PTUR.inertial_properties_perturbed.T;
        output_containers.T_estimated.PTLR(j, :, 1) = algorithm_container.PTLR.inertial_properties_perturbed.T;
        
        output_containers.R_estimated.SVDLS(j, :, :) = algorithm_container.SVDLS.inertial_properties_SVDLS.R;
        output_containers.R_estimated.PTNR(j, :, :) = algorithm_container.PTNR.inertial_properties_perturbed.R;
        output_containers.R_estimated.PTUR(j, :, :) = algorithm_container.PTUR.inertial_properties_perturbed.R;
        output_containers.R_estimated.PTLR(j, :, :) = algorithm_container.PTLR.inertial_properties_perturbed.R;
        
        output_containers.T_reference.SVDLS(j, :, 1) = algorithm_container.SVDLS.inertial_properties_reference.T;
        output_containers.R_reference.SVDLS(j, :, :) = algorithm_container.SVDLS.inertial_properties_reference.R;
        
        for k=1:Num_AL_Thigh
            output_containers.AL_estimated.SVDLS(j, k, :) = algorithm_container.SVDLS.anatomical_landmarks_thigh_data(k).vector_global_reconstructed;
            output_containers.AL_estimated.PTNR(j, k, :) = algorithm_container.PTNR.anatomical_landmarks_thigh_data(k).vector_global_reconstructed;
            output_containers.AL_estimated.PTUR(j, k, :) = algorithm_container.PTUR.anatomical_landmarks_thigh_data(k).vector_global_reconstructed;
            output_containers.AL_estimated.PTLR(j, k, :) = algorithm_container.PTLR.anatomical_landmarks_thigh_data(k).vector_global_reconstructed;

            output_containers.AL_reference.SVDLS(j, k, :) = algorithm_container.SVDLS.anatomical_landmarks_thigh_data(k).vector_global_reference;
            output_containers.AL_reference.PTNR(j, k, :) = algorithm_container.PTNR.anatomical_landmarks_thigh_data(k).vector_global_reference;
            output_containers.AL_reference.PTUR(j, k, :) = algorithm_container.PTUR.anatomical_landmarks_thigh_data(k).vector_global_reference;
            output_containers.AL_reference.PTLR(j, k, :) = algorithm_container.PTLR.anatomical_landmarks_thigh_data(k).vector_global_reference;
        end
        
        output_containers.AF_Translation_estimated.SVDLS(j, :, 1) = algorithm_container.SVDLS.anatomical_frame.T;
        output_containers.AF_Translation_estimated.PTNR(j, :, 1) = algorithm_container.PTNR.anatomical_frame.T;
        output_containers.AF_Translation_estimated.PTUR(j, :, 1) = algorithm_container.PTUR.anatomical_frame.T;
        output_containers.AF_Translation_estimated.PTLR(j, :, 1) = algorithm_container.PTLR.anatomical_frame.T;
        
        output_containers.AF_Rotation_estimated.SVDLS(j, :, :) = algorithm_container.SVDLS.anatomical_frame.R;
        output_containers.AF_Rotation_estimated.PTNR(j, :, :) = algorithm_container.PTNR.anatomical_frame.R;
        output_containers.AF_Rotation_estimated.PTUR(j, :, :) = algorithm_container.PTUR.anatomical_frame.R;
        output_containers.AF_Rotation_estimated.PTLR(j, :, :) = algorithm_container.PTLR.anatomical_frame.R;
        
        output_containers.AF_Translation_reference.SVDLS(j, :, 1) = algorithm_container.SVDLS.anatomical_frame.T_ref;
        output_containers.AF_Rotation_reference.SVDLS(j, :, :) = algorithm_container.SVDLS.anatomical_frame.R_ref;
        
    end
end

function generate_output_plots(input_struct, output_containers)
    input_size = input_struct.input_size;
    input_data = input_struct.input_data;

    configuration_idx = input_data.configuration_idx;
    N = input_size.N;
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    ToeOff = 0.6;
    
    color_list = ['k', 'b', 'r', 'm'];
    alg_list = ["SVDLS", "PTUR"];

    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\STA_time_series_pose\N=" + string(N) + "\";
    
    path_CM = path_img + "CM_plots\";
    figure('units','normalized','outerposition',[0 0 1 1])
    xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    for n=1:length(alg_list)
        current_alg = alg_list(n);
        current_data = output_containers.CM_offset.(current_alg);
        if current_alg == "PTNR"
            scatter(fraction_completed_gait_cycle, current_data, color_list(n), 'DisplayName', current_alg);
        else
            plot(fraction_completed_gait_cycle, current_data, color_list(n), 'DisplayName', current_alg);
        end
    end
    grid;
    legend;
    xlabel("Fraction completed gait cycle");
    ylabel("Error (cm)");
    title("Center of Mass offset , IDX = " + string(configuration_idx));
    filename_CM = path_CM + "CM_plot_" + string(configuration_idx) + ".png";
    saveas(gcf, filename_CM);
    close(gcf);

    path_Eigenvector = path_img + "Eigenvector_plots\";
    figure('units','normalized','outerposition',[0 0 1 1])
    for j=1:3
        subplot(3,1,j);
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_data = output_containers.eigenvector_offset.(current_alg);
            plot(fraction_completed_gait_cycle, current_data(:,j), color_list(n), 'DisplayName', current_alg);
        end
        grid;
        if j==1
            legend;
        end
        if j==3
            xlabel("Fraction completed gait cycle");
        end
        ylabel("Error (degrees)");
        title("Eigenvector " + string(j));
    end
    suptitle("Eigenvector Offset , IDX = " + string(configuration_idx));
    filename_Eigenvector = path_Eigenvector + "Eigenvector_plot_" + string(configuration_idx) + ".png";
    saveas(gcf, filename_Eigenvector);
    close(gcf);
    
    AL_name_list = ["Lateral Epicondyle";
                 "Medial Epicondyle";
                 "Greater Trochanter";
                 "Femoral Head"];
             
    path_AL = path_img + "AL_plots\";
    figure('units','normalized','outerposition',[0 0 1 1])
    for j=1:4
        subplot(2,2,j);
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_data = output_containers.AL_offset.(current_alg);
            if current_alg == "PTNR"
                scatter(fraction_completed_gait_cycle, current_data(:,j), color_list(n), 'DisplayName', current_alg);
            else
                plot(fraction_completed_gait_cycle, current_data(:,j), color_list(n), 'DisplayName', current_alg);
            end
        end
        grid;
        if j==2
            legend;
        end
        if j>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(j,2) == 0)
            ylabel("Error (cm)");
        end
        title(AL_name_list(j));
    end
    suptitle("Anatomical Landmark Offset , IDX = " + string(configuration_idx))
    filename_AL = path_AL + "AL_plot_" + string(configuration_idx) + ".png";
    saveas(gcf, filename_AL);
    close(gcf);
    
    path_ClusterPt = path_img + "ClusterPt_plots\";
    figure('units','normalized','outerposition',[0 0 1 1])
    for j=1:N
        subplot(2, N/2, j)
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_data = output_containers.cluster_offset.(current_alg);
            if current_alg == "PTNR"
                scatter(fraction_completed_gait_cycle, current_data(:,j), color_list(n), 'DisplayName', current_alg);
            else
                plot(fraction_completed_gait_cycle, current_data(:,j), color_list(n), 'DisplayName', current_alg);
            end
        end
        grid;
        if j==3
            legend;
        end
        if j>3
            xlabel("Fraction completed gait cycle");
        end
        if or(j==1, j==4)
            ylabel("Error (cm)");
        end
        title("Cluster Point " + string(j));
    end
    suptitle("Cluster Point Global Position Reconstruction Error , IDX = " + string(configuration_idx));
    filename_ClusterPt = path_ClusterPt + "ClusterPt_plot_" + string(configuration_idx) + ".png";
    saveas(gcf, filename_ClusterPt);
    close(gcf);

end

function algorithm_containerPT = Perturbation_Theory_processing(algorithm_containerPT)
    cylinder_data = algorithm_containerPT.cylinder_data;
    current_pose = algorithm_containerPT.current_pose;
    N = algorithm_containerPT.N;
    inertial_properties_noisy = algorithm_containerPT.inertial_properties_noisy_original;
    inertial_properties_noisy_copy = algorithm_containerPT.inertial_properties_noisy;
    inertial_properties_reference = algorithm_containerPT.inertial_properties_reference;
    mass_distribution_noisy = algorithm_containerPT.mass_distribution_noisy;
    mass_distribution_reference = algorithm_containerPT.mass_distribution_reference;
    cm_dist_threshold = algorithm_containerPT.cm_dist_threshold;
    reflection_setting = algorithm_containerPT.reflection_setting;

    lambda_start = nan(3,1);
    lambda_target = nan(3,1);
    vectors_start = nan(3,3);
    for j=1:3
        lambda_start(j) = inertial_properties_noisy_copy.PrincipalMOI(j).value;
        lambda_target(j) = inertial_properties_reference.PrincipalMOI(j).value;
        vectors_start(:,j) = inertial_properties_noisy_copy.PrincipalAxis(j).vector;
    end
    [vectors_final, ~, inertiaTensor_final, ~] = iterative_perturbation(vectors_start, lambda_start, lambda_target);
    m0 = algorithm_containerPT.m0;
    options_nlcon = optimoptions('fmincon','Display','iter','Algorithm','sqp');
    func = @(m)nLsysMD(m, mass_distribution_noisy, inertiaTensor_final, N);
    funcRef = @(m)nLsysMD(m, mass_distribution_reference, inertial_properties_reference.MOI, N);
    funcSum = @(m)nLsysMDsum(m, mass_distribution_noisy, inertiaTensor_final, N);
    cm_nonlcon = @(m)cm_constraint(m, mass_distribution_noisy, inertial_properties_noisy, N, cm_dist_threshold);
    
    mass_tol = 10^-10;
    cm_reflect_bool = false;
    
    if all(abs(fsolve(funcRef, m0) - m0) < mass_tol*ones(N,1))
        try
            m_star = fsolve(func, m0);
            cm_constraint_val = cm_constraint(m_star, mass_distribution_noisy, inertial_properties_noisy, N, cm_dist_threshold);
            if cm_constraint_val > 0
                m_star_cm = fmincon(funcSum, m_star, [], [], [], [], zeros(N,1), [], cm_nonlcon, options_nlcon);
                cm_constraint_val_new = cm_constraint(m_star_cm, mass_distribution_noisy, inertial_properties_noisy, N, cm_dist_threshold);
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
    
    cm_vector_noisy = inertial_properties_noisy.CM.vector;
    cm_vector_reflected = cm_reflect(m_star, mass_distribution_noisy, inertial_properties_noisy, N);

    cm_vector_relative_global = cm_vector_noisy - cm_vector_reflected;
    deltaY_offset = dot(cm_vector_relative_global , vectors_final(:,1));
    algorithm_containerPT.deltaY_offset = deltaY_offset;
    if reflection_setting == "UR"
        if deltaY_offset > 0 
            cm_reflect_bool = true;
        end
    elseif reflection_setting == "LR"
        if deltaY_offset < 0
            cm_reflect_bool = true;
        end
    end
    
    disp(m_star);
    mass_distribution_perturbed = mass_distribution_noisy;
    for i=1:N
        mass_distribution_perturbed(i).mass = m_star(i);
    end
    
    [inertial_properties_perturbed, mass_distribution_perturbed] = compute_inertial_properties(mass_distribution_perturbed, N);
    inertial_properties_perturbed = compute_inertial_properties_reflected(cm_reflect_bool, cm_vector_reflected, inertial_properties_perturbed);
    inertial_properties_perturbed = add_inertial_properties_perturbation(vectors_final, inertiaTensor_final, inertial_properties_perturbed);
    inertial_properties_perturbed_offset = calculate_offset(inertial_properties_reference, inertial_properties_perturbed);
    algorithm_containerPT.inertial_properties_perturbed_offset = inertial_properties_perturbed_offset;
    
    algorithm_containerPT.mass_distribution_perturbed = mass_distribution_perturbed;
    algorithm_containerPT.inertial_properties_perturbed = inertial_properties_perturbed;
    algorithm_containerPT = compute_Final_localOffset_perturbed(algorithm_containerPT);
    
    reference_pose = struct;
    reference_pose.T = inertial_properties_reference.T;
    reference_pose.R = inertial_properties_reference.R;
    
    optimized_pose_perturbed = struct;
    optimized_pose_perturbed.T = inertial_properties_perturbed.T;
    optimized_pose_perturbed.R = inertial_properties_perturbed.R;
    
    [anatomical_landmarks_thigh_data_perturbed, anatomical_frame_perturbed] = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose_perturbed, current_pose);
    algorithm_containerPT.anatomical_landmarks_thigh_data = anatomical_landmarks_thigh_data_perturbed;
    algorithm_containerPT.anatomical_frame = anatomical_frame_perturbed;
end

function algorithm_containerSVDLS = SVDLS_processing(algorithm_containerSVDLS)
    cylinder_data = algorithm_containerSVDLS.cylinder_data;
    inertial_properties_reference = algorithm_containerSVDLS.inertial_properties_reference;

    [T_SVDLS, R_SVDLS] = pose_estimate_SVD_LS(algorithm_containerSVDLS);
    inertial_properties_SVDLS = add_inertial_properties_SVDLS(T_SVDLS, R_SVDLS);
    inertial_properties_SVDLS_offset = calculate_SVDLS_offset(inertial_properties_reference, inertial_properties_SVDLS);
    algorithm_containerSVDLS.inertial_properties_SVDLS = inertial_properties_SVDLS;
    algorithm_containerSVDLS.inertial_properties_SVDLS_offset = inertial_properties_SVDLS_offset;
    algorithm_containerSVDLS = compute_Final_localOffset_SVDLS(algorithm_containerSVDLS);
    
    reference_pose = struct;
    reference_pose.T = inertial_properties_reference.T;
    reference_pose.R = inertial_properties_reference.R;
    
    optimized_pose_SVDLS = struct;
    optimized_pose_SVDLS.T = inertial_properties_SVDLS.T;
    optimized_pose_SVDLS.R = inertial_properties_SVDLS.R;
    
    current_pose = algorithm_containerSVDLS.current_pose;
    
    [anatomical_landmarks_thigh_data_SVDLS, anatomical_frame_SVDLS] = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose_SVDLS, current_pose);
    algorithm_containerSVDLS.anatomical_landmarks_thigh_data = anatomical_landmarks_thigh_data_SVDLS;
    algorithm_containerSVDLS.anatomical_frame = anatomical_frame_SVDLS;
    
end

function cm_dist_threshold = compute_cm_dist_threshold_STA(fraction_completed_gait_cycle)
    if fraction_completed_gait_cycle <= 0.5
        cm_dist_threshold = 0.7 + 2.2*fraction_completed_gait_cycle;
    else
        cm_dist_threshold = 2.9 - 2.2*fraction_completed_gait_cycle;
    end
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

function output_containers = generate_output_containers(input_size)
    output_containers = struct;
    
    algorithm_str_list = ["SVDLS", "PCTUO", "PTNR", "PTUR", "PTLR"];
    NumAlgorithms = length(algorithm_str_list);
    
    NumFrames = input_size.NumFrames;
    Num_AL_Thigh = input_size.Num_AL_Thigh;
    N = input_size.N;
    
    NumEigenvectors = 3;
    
    eigenvalueNorm_offset = struct;
    eigenvector_offset = struct;
    eps_star = struct;
    mass_redistribution = struct;
    CM_offset = struct;
    cluster_offset = struct;
    AL_offset = struct;
    deltaY_offset = struct;
    
    T_estimated = struct;
    R_estimated = struct;
    
    T_reference = struct;
    R_reference = struct;
    
    AL_estimated = struct;
    AL_reference = struct;
    
    AF_Translation_estimated = struct;
    AF_Rotation_estimated = struct;
    
    AF_Translation_reference = struct;
    AF_Rotation_reference = struct;
    
    for n=1:NumAlgorithms
        algorithm_str_current = algorithm_str_list(n);
        eigenvalueNorm_offset.(algorithm_str_current) = nan(NumFrames, 1);
        eigenvector_offset.(algorithm_str_current) = nan(NumFrames, NumEigenvectors, 1);
        eps_star.(algorithm_str_current) = nan(NumFrames, 1);
        mass_redistribution.(algorithm_str_current) = nan(NumFrames, N, 1);
        CM_offset.(algorithm_str_current) = nan(NumFrames, 1);
        cluster_offset.(algorithm_str_current) = nan(NumFrames, N, 1);
        AL_offset.(algorithm_str_current) = nan(NumFrames, Num_AL_Thigh, 1);
        deltaY_offset.(algorithm_str_current) = nan(NumFrames, 1);
        
        T_estimated.(algorithm_str_current) = nan(NumFrames, 3, 1);
        R_estimated.(algorithm_str_current) = nan(NumFrames, 3, 3);
        
        T_reference.(algorithm_str_current) = nan(NumFrames, 3, 1);
        R_reference.(algorithm_str_current) = nan(NumFrames, 3, 3);
        
        AL_estimated.(algorithm_str_current) = nan(NumFrames, Num_AL_Thigh, 3);
        AL_reference.(algorithm_str_current) = nan(NumFrames, Num_AL_Thigh, 3);
        
        AF_Translation_estimated.(algorithm_str_current) = nan(NumFrames, 3, 1);
        AF_Rotation_estimated.(algorithm_str_current) = nan(NumFrames, 3, 3);
        
        AF_Translation_reference.(algorithm_str_current) = nan(NumFrames, 3, 1);
        AF_Rotation_reference.(algorithm_str_current) = nan(NumFrames, 3, 3);
    end
    
    output_containers.eigenvalueNorm_offset = eigenvalueNorm_offset;
    output_containers.eigenvector_offset = eigenvector_offset;
    output_containers.eps_star = eps_star;
    output_containers.mass_redistribution = mass_redistribution;
    output_containers.CM_offset = CM_offset;
    output_containers.cluster_offset = cluster_offset;
    output_containers.AL_offset = AL_offset;
    output_containers.deltaY_offset = deltaY_offset;
    
    output_containers.T_estimated = T_estimated;
    output_containers.R_estimated = R_estimated;
    
    output_containers.T_reference = T_reference;
    output_containers.R_reference = R_reference;
    
    output_containers.AL_estimated = AL_estimated;
    output_containers.AL_reference = AL_reference;
    
    output_containers.AF_Translation_estimated = AF_Translation_estimated;
    output_containers.AF_Rotation_estimated = AF_Rotation_estimated;
    
    output_containers.AF_Translation_reference = AF_Translation_reference;
    output_containers.AF_Rotation_reference = AF_Rotation_reference;
    
end

function input_settings = generate_input_settings
    input_settings = struct;
    
    initial_axes_reversal = false;
    initial_axes_perturbed_reversal = true;
    optimized_axes_reversal = true;
    
    input_settings.initial_axes_reversal = initial_axes_reversal;
    input_settings.initial_axes_perturbed_reversal = initial_axes_perturbed_reversal;
    input_settings.optimized_axes_reversal = optimized_axes_reversal;
    
    input_settings.m0_switch_idx = 2;
    
end

function input_size = generate_input_size_struct(N, NumConfigurations, box_markers_gait)
    input_size = struct;
    
    fraction_completed_gait_cycle = box_markers_gait(1).fraction_completed_gait_cycle;
    NumFrames = length(fraction_completed_gait_cycle);
    Num_AL_Thigh = 4;
    
    input_size.N = N;
    input_size.fraction_completed_gait_cycle = fraction_completed_gait_cycle;
    input_size.NumFrames = NumFrames;
    input_size.Num_AL_Thigh = Num_AL_Thigh;
    input_size.NumAngles = 3;
    input_size.NumConfigurations = NumConfigurations;
end

function noisy_distribution_STA = generate_noisy_distribution_STA(mass_distribution_reference, frame_idx, box_markers_gait, N)
    noisy_distribution_STA = struct;
    
    for j=1:N
        STA_data = box_markers_gait(j).vector(frame_idx, :)';
        noise_mm = rotx(90)*STA_data;
        noise = noise_mm/10;
        noisy_distribution_STA(j).mass = 1;
        noisy_vector = mass_distribution_reference(j).vector + noise;
        noisy_distribution_STA(j).X = noisy_vector(1);
        noisy_distribution_STA(j).Y = noisy_vector(2);
        noisy_distribution_STA(j).Z = noisy_vector(3);
        noisy_distribution_STA(j).vector = noisy_vector;
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
    try
        [V, D] = eig(I);
    catch
        V = nan(3,3);
        D = diag(nan(3,1));
    end
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

function inertial_properties_offset = calculate_offset(inertial_properties_reference, inertial_properties_noisy)
    inertial_properties_offset = struct;
    for j=1:3
        inertial_properties_offset.angle(j) = rad2deg(acos(dot(inertial_properties_reference.PrincipalAxis(j).vector, inertial_properties_noisy.PrincipalAxis(j).vector)));
        inertial_properties_offset.PrincipalMOI(j) = inertial_properties_reference.PrincipalMOI(j).value - inertial_properties_noisy.PrincipalMOI(j).value;
    end
    inertial_properties_offset.PrincipalMOINorm = inertial_properties_reference.PrincipalMOINorm - inertial_properties_noisy.PrincipalMOINorm;
    inertial_properties_offset.CMVector = inertial_properties_reference.CM.vector - inertial_properties_noisy.CM.vector;
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

function [algorithm_container, local_offset, local_angular_offset, local_offset_vector] = compute_localOffset(algorithm_container)
    mass_distribution_reference = algorithm_container.mass_distribution_reference;
    mass_distribution_noisy = algorithm_container.mass_distribution_noisy;
    inertial_properties_reference = algorithm_container.inertial_properties_reference;
    inertial_properties_noisy = algorithm_container.inertial_properties_noisy;
    T_reference = inertial_properties_reference.T;
    R_reference = inertial_properties_reference.R;
    T_noisy = inertial_properties_noisy.T;
    R_noisy = inertial_properties_noisy.R;
    N = algorithm_container.N;
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
    algorithm_container.mass_distribution_reference = mass_distribution_reference;
    algorithm_container.mass_distribution_noisy = mass_distribution_noisy;
end

function [T_val, R_val] = pose_estimate_SVD_LS(algorithm_container)
    N = algorithm_container.N;
    
    global_vectors = struct;
    local_vectors = struct;
    for c=1:N
        global_vectors(c).vec = algorithm_container.mass_distribution_noisy(c).vector;
        global_vectors(c).x = global_vectors(c).vec(1);
        global_vectors(c).y = global_vectors(c).vec(2);
        global_vectors(c).z = global_vectors(c).vec(3);

        local_vectors(c).vec = algorithm_container.mass_distribution_reference(c).local_displacement;
        local_vectors(c).x = local_vectors(c).vec(1);
        local_vectors(c).y = local_vectors(c).vec(2);
        local_vectors(c).z = local_vectors(c).vec(3);
    end
    
    global_vectors_x_list = nan(N,1);
    global_vectors_y_list = nan(N,1);
    global_vectors_z_list = nan(N,1);
    
    local_vectors_x_list = nan(N,1);
    local_vectors_y_list = nan(N,1);
    local_vectors_z_list = nan(N,1);
    
    for c=1:N
        global_vectors_x_list(c) = global_vectors(c).x;
        global_vectors_y_list(c) = global_vectors(c).y;
        global_vectors_z_list(c) = global_vectors(c).z;
        
        local_vectors_x_list(c) = local_vectors(c).x;
        local_vectors_y_list(c) = local_vectors(c).y;
        local_vectors_z_list(c) = local_vectors(c).z;
    end
    
    global_vectors_mean = struct;
    global_vectors_mean.x = mean(global_vectors_x_list, 'omitnan');
    global_vectors_mean.y = mean(global_vectors_y_list, 'omitnan');
    global_vectors_mean.z = mean(global_vectors_z_list, 'omitnan');
    global_vectors_mean.vec = [global_vectors_mean.x; global_vectors_mean.y; global_vectors_mean.z];
    
    local_vectors_mean = struct;
    local_vectors_mean.x = mean(local_vectors_x_list, 'omitnan');
    local_vectors_mean.y = mean(local_vectors_y_list, 'omitnan');
    local_vectors_mean.z = mean(local_vectors_z_list, 'omitnan');
    local_vectors_mean.vec = [local_vectors_mean.x; local_vectors_mean.y; local_vectors_mean.z];

    q_global = struct;
    q_local = struct;
    h_matrix = zeros(3,3);
    for c=1:N
        q_global(c).x = global_vectors(c).x - global_vectors_mean.x;
        q_global(c).y = global_vectors(c).y - global_vectors_mean.y;
        q_global(c).z = global_vectors(c).z - global_vectors_mean.z;
        q_global(c).vec = [q_global(c).x; q_global(c).y; q_global(c).z];
        
        q_local(c).x = local_vectors(c).x - local_vectors_mean.x;
        q_local(c).y = local_vectors(c).y - local_vectors_mean.y;
        q_local(c).z = local_vectors(c).z - local_vectors_mean.z;
        q_local(c).vec = [q_local(c).x; q_local(c).y; q_local(c).z];

        h_matrix = h_matrix + (q_local(c).vec * transpose(q_global(c).vec));
    end

    sanityCheck_tol = 10^-8;
    X = nan(3,3);
    if ~any(any(isnan(h_matrix)))
        [U, ~, V] = svd(h_matrix);
        X = V*U';
    end

    if ~all(all(isnan(X)))
        if abs(det(X) - 1) > sanityCheck_tol
            if abs(det(X) + 1) < sanityCheck_tol
                V_new = horzcat(V(:,1:2), -V(:,3));
                X = V_new*U';
            else
                X = nan(3,3);
            end
        end
    end

    R_val = X;
    T_val = global_vectors_mean.vec - R_val*local_vectors_mean.vec;
   
end

function inertial_properties_SVDLS = add_inertial_properties_SVDLS(T_SVDLS, R_SVDLS)

    inertial_properties_SVDLS = struct;
    
    inertial_properties_SVDLS.CM = struct;
    inertial_properties_SVDLS.CM.x = T_SVDLS(1);
    inertial_properties_SVDLS.CM.y = T_SVDLS(2);
    inertial_properties_SVDLS.CM.z = T_SVDLS(3);
    inertial_properties_SVDLS.CM.vector = [inertial_properties_SVDLS.CM.x; inertial_properties_SVDLS.CM.y; inertial_properties_SVDLS.CM.z];
    inertial_properties_SVDLS.T = T_SVDLS;
    
    inertial_properties_SVDLS.PrincipalAxis = struct;
    inertial_properties_SVDLS.PrincipalAxis(1).vector = R_SVDLS(:,1);
    inertial_properties_SVDLS.PrincipalAxis(2).vector = R_SVDLS(:,2);
    inertial_properties_SVDLS.PrincipalAxis(3).vector = R_SVDLS(:,3);
    inertial_properties_SVDLS.R = R_SVDLS;
end

function inertial_properties_offset = calculate_SVDLS_offset(inertial_properties_reference, inertial_properties_SVDLS)
    inertial_properties_offset = struct;
    for j=1:3
        inertial_properties_offset.angle(j) = rad2deg(acos(dot(inertial_properties_reference.PrincipalAxis(j).vector, inertial_properties_SVDLS.PrincipalAxis(j).vector)));
    end
    inertial_properties_offset.CMVector = inertial_properties_reference.CM.vector - inertial_properties_SVDLS.CM.vector;
end

function algorithm_container = compute_Final_localOffset_SVDLS(algorithm_container)
    T = algorithm_container.inertial_properties_SVDLS.T;
    R = algorithm_container.inertial_properties_SVDLS.R;
    N = algorithm_container.N;
    for k=1:N
        current_global_vector = algorithm_container.mass_distribution_SVDLS(k).vector;
        current_ref_vector = algorithm_container.mass_distribution_reference(k).local_displacement;
        original_global_vector = algorithm_container.mass_distribution_reference(k).vector;
        current_local_vector = global2local(T, R, current_global_vector);
        global_vector_reconstructed = local2global(T, R, current_ref_vector);
        global_offset = norm(global_vector_reconstructed - original_global_vector);
        algorithm_container.mass_distribution_SVDLS(k).local_displacement = current_local_vector;
        algorithm_container.mass_distribution_SVDLS(k).local_offset = norm(current_local_vector - current_ref_vector);
        algorithm_container.mass_distribution_SVDLS(k).global_reconstructed = global_vector_reconstructed;
        algorithm_container.mass_distribution_SVDLS(k).global_offset = global_offset;
    end
    
    inertial_properties_noisy_reversed = axis_reversal_original(algorithm_container.inertial_properties_noisy, algorithm_container.inertial_properties_offset);
    algorithm_container.inertial_properties_noisy = inertial_properties_noisy_reversed;
    algorithm_container.inertial_properties_offset = calculate_offset(algorithm_container.inertial_properties_reference, algorithm_container.inertial_properties_noisy);
    T_noisy = algorithm_container.inertial_properties_noisy.T;
    R_noisy = algorithm_container.inertial_properties_noisy.R;
    for k=1:N
        current_ref_vector = algorithm_container.mass_distribution_reference(k).local_displacement;
        original_global_vector = algorithm_container.mass_distribution_reference(k).vector;
        global_vector_reconstructed_noisy = local2global(T_noisy, R_noisy, current_ref_vector);
        global_offset_noisy = norm(global_vector_reconstructed_noisy - original_global_vector);
        algorithm_container.mass_distribution_noisy(k).global_reconstructed = global_vector_reconstructed_noisy;
        algorithm_container.mass_distribution_noisy(k).global_offset = global_offset_noisy;
    end
end

function [anatomical_landmarks_thigh_data, anatomical_frame] = compute_anatomical_landmarks_thigh_data(cylinder_data, reference_pose, optimized_pose, current_pose)
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
    
    T_current = current_pose.T;
    R_current = current_pose.R;
    
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
            x_global_reference = 0;
        else
            x_global_reference = X(current_height_idx, current_plane_idx);
        end
        y_global_reference = Y(current_height_idx, current_plane_idx);
        z_global_reference = Z(current_height_idx, current_plane_idx);
        vector_standard_reference = [x_global_reference; y_global_reference; z_global_reference];
        vector_STA_reference = rotx(-90)*vector_standard_reference;
        vector_global_reference = T_current + R_current*vector_STA_reference;
        anatomical_landmarks_thigh_data(j).vector_global_reference = vector_global_reference;
        anatomical_landmarks_thigh_data(j).vector_local_reference = global2local_perturbed(T_ref, R_ref, anatomical_landmarks_thigh_data(j).vector_global_reference);
        anatomical_landmarks_thigh_data(j).vector_global_reconstructed = local2global_perturbed(T_opt, R_opt, anatomical_landmarks_thigh_data(j).vector_local_reference);
        anatomical_landmarks_thigh_data(j).global_offset = norm(anatomical_landmarks_thigh_data(j).vector_global_reconstructed - anatomical_landmarks_thigh_data(j).vector_global_reference);
    end
    
    current_reconstructed_position_val = struct;
    current_reconstructed_position_val.("left_lat_femoral_epicondyle") = anatomical_landmarks_thigh_data(1).vector_global_reconstructed;
    current_reconstructed_position_val.("left_med_femoral_epicondyle") = anatomical_landmarks_thigh_data(2).vector_global_reconstructed;
    current_reconstructed_position_val.("left_femoral_head") = anatomical_landmarks_thigh_data(4).vector_global_reconstructed;
    [T_AF, R_AF, unit_vectors_AF] = pose_estimate_femoral_frame(current_reconstructed_position_val);
    
    current_reference_position_val = struct;
    current_reference_position_val.("left_lat_femoral_epicondyle") = anatomical_landmarks_thigh_data(1).vector_global_reference;
    current_reference_position_val.("left_med_femoral_epicondyle") = anatomical_landmarks_thigh_data(2).vector_global_reference;
    current_reference_position_val.("left_femoral_head") = anatomical_landmarks_thigh_data(4).vector_global_reference;
    [T_AF_ref, R_AF_ref, unit_vectors_AF_ref] = pose_estimate_femoral_frame(current_reference_position_val);
    
    anatomical_frame = struct;
    anatomical_frame.T = T_AF;
    anatomical_frame.R = R_AF;
    anatomical_frame.unit_vectors = unit_vectors_AF;
    
    anatomical_frame.T_ref = T_AF_ref;
    anatomical_frame.R_ref = R_AF_ref;
    anatomical_frame.unit_vectors_ref = unit_vectors_AF_ref;
end

function PCT_SVDLS_performance = evaluate_PCT_SVDLS(inertial_properties_SVDLS_offset, inertial_properties_offset)
    PCT_SVDLS_performance = struct;
    
    PCT_SVDLS_performance.CM_distance.boolVal = norm(inertial_properties_SVDLS_offset.CMVector) < norm(inertial_properties_offset.CMVector);
    PCT_SVDLS_performance.CM_distance.oldNorm = norm(inertial_properties_offset.CMVector);
    PCT_SVDLS_performance.CM_distance.newNorm = norm(inertial_properties_SVDLS_offset.CMVector);
    
    angleList = false(3,1);
    angleOldList = nan(3,1);
    angleNewList = nan(3,1);
    angleDiffList = nan(3,1);
    for m=1:3
        angleOld = inertial_properties_offset.angle(m);
        angleNew = inertial_properties_SVDLS_offset.angle(m);

        if angleOld > 90
            angleOld = 180 - inertial_properties_offset.angle(m);
        end

        if angleNew > 90
            angleNew = 180 - inertial_properties_SVDLS_offset.angle(m);
        end

        angleOldList(m) = angleOld;
        angleNewList(m) = angleNew;
        angleDiffList(m) = angleOld - angleNew;
        if angleNew < angleOld
            angleList(m) = true;
        end
    end
    PCT_SVDLS_performance.angle.boolVal = all(angleList);
    PCT_SVDLS_performance.angle.oldList = angleOldList;
    PCT_SVDLS_performance.angle.newList = angleNewList;
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

function [c, ceq] = cm_constraint(m, mass_distribution, inertial_properties, N, cm_dist_threshold)
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

function algorithm_containerPT = compute_Final_localOffset_perturbed(algorithm_containerPT)
    T = algorithm_containerPT.inertial_properties_perturbed.T;
    R = algorithm_containerPT.inertial_properties_perturbed.R;
    N = algorithm_containerPT.N;
    for k=1:N
        current_global_vector = algorithm_containerPT.mass_distribution_perturbed(k).vector;
        current_ref_vector = algorithm_containerPT.mass_distribution_reference(k).local_displacement;
        original_global_vector = algorithm_containerPT.mass_distribution_reference(k).vector;
        current_local_vector = global2local_perturbed(T, R, current_global_vector);
        global_vector_reconstructed = local2global_perturbed(T, R, current_ref_vector);
        global_offset = norm(global_vector_reconstructed - original_global_vector);
        algorithm_containerPT.mass_distribution_perturbed(k).local_displacement = current_local_vector;
        algorithm_containerPT.mass_distribution_perturbed(k).local_offset = norm(current_local_vector - current_ref_vector);
        algorithm_containerPT.mass_distribution_perturbed(k).global_reconstructed = global_vector_reconstructed;
        algorithm_containerPT.mass_distribution_perturbed(k).global_offset = global_offset;
    end
    
    inertial_properties_noisy_reversed = axis_reversal_original(algorithm_containerPT.inertial_properties_noisy, algorithm_containerPT.inertial_properties_offset);
    algorithm_containerPT.inertial_properties_noisy = inertial_properties_noisy_reversed;
    algorithm_containerPT.inertial_properties_offset = calculate_offset(algorithm_containerPT.inertial_properties_reference, algorithm_containerPT.inertial_properties_noisy);
    T_noisy = algorithm_containerPT.inertial_properties_noisy.T;
    R_noisy = algorithm_containerPT.inertial_properties_noisy.R;
    for k=1:N
        current_ref_vector = algorithm_containerPT.mass_distribution_reference(k).local_displacement;
        original_global_vector = algorithm_containerPT.mass_distribution_reference(k).vector;
        global_vector_reconstructed_noisy = local2global(T_noisy, R_noisy, current_ref_vector);
        global_offset_noisy = norm(global_vector_reconstructed_noisy - original_global_vector);
        algorithm_containerPT.mass_distribution_noisy(k).global_reconstructed = global_vector_reconstructed_noisy;
        algorithm_containerPT.mass_distribution_noisy(k).global_offset = global_offset_noisy;
    end
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

function create_animation_RF(input_struct, output_containers, reference_frame_str)
    
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\STA_time_series_pose\N=" + string(input_struct.input_size.N) + "\";
    
    movingView_flag = input_struct.visualization_settings.movingView_flag_val;
    animation_flag = input_struct.visualization_settings.animation_flag;
    fig_scale = input_struct.visualization_settings.fig_scale;
    
    substruct_str = struct;
    
    if reference_frame_str == "TF"
        
        substruct_str.R_estimated =  "R_estimated";
        substruct_str.T_estimated =  "T_estimated";
        
        substruct_str.R_reference =  "R_reference";
        substruct_str.T_reference =  "T_reference";
        
        folder_name_start = "Movement_Videos_PCT_Thigh_TF\movement_subNo_";
        if movingView_flag
            folder_name_start = "Movement_MovingPtView_Videos_PCT_Thigh_TF\movement_subNo_";
        end
        folder_name_start = path_img + folder_name_start;
        
    elseif reference_frame_str == "AF"
        
        substruct_str.R_estimated = "AF_Rotation_estimated";
        substruct_str.T_estimated = "AF_Translation_estimated";
        
        substruct_str.R_reference = "AF_Rotation_reference";
        substruct_str.T_reference = "AF_Translation_reference";
        
        folder_name_start = "Movement_Videos_PCT_Thigh_AF\movement_subNo_";
        if movingView_flag
            folder_name_start = "Movement_MovingPtView_Videos_PCT_Thigh_AF\movement_subNo_";
        end
        folder_name_start = path_img + folder_name_start;
        
    end
    
    if animation_flag
        
        input_data = input_struct.input_data;
        input_size = input_struct.input_size;
        
        configuration_idx = input_data.configuration_idx;
        pose_gait = input_data.pose_gait;
        cylinder_data = input_data.cylinder_data;
        mass_distribution_reference_standard_pose = input_data.mass_distribution_reference;
        box_markers_gait = input_data.box_markers_gait;
        
        num_measurements = input_size.NumFrames;
        N = input_size.N;
        
        original_view_vector = fig_scale.view_vector;
        original_azimuth = original_view_vector(1);
        original_elevation = original_view_vector(2);
        dTheta = 0;
        if ~isnan(num_measurements)
            if num_measurements > 1
                dTheta = rad2deg((2*pi)/(num_measurements-1));
            end
        end
        
        my_writer = VideoWriter(folder_name_start + string(configuration_idx), 'MPEG-4');
        my_writer.FrameRate = 50;
        open(my_writer);
        
        figh = figure;
        figh.WindowState = 'maximized';

        for k=1:num_measurements
            unit_vector_local = struct;
            
            unit_vector_local.REF = generate_unit_vectors_local(output_containers.(substruct_str.R_reference).SVDLS(k, :, :));
            unit_vector_local.SVDLS = generate_unit_vectors_local(output_containers.(substruct_str.R_estimated).SVDLS(k, :, :));
            unit_vector_local.PTUR = generate_unit_vectors_local(output_containers.(substruct_str.R_estimated).PTUR(k, :, :));
            
            current_pose = pose_gait(k);
            cylinder_data_pose = generate_cylinder_data_pose(cylinder_data, current_pose);
            X = cylinder_data_pose.X;
            Y = cylinder_data_pose.Y;
            Z = cylinder_data_pose.Z;
            
            mass_distribution_noisy_standard_pose = generate_noisy_distribution_STA(mass_distribution_reference_standard_pose, k, box_markers_gait, N);
            mass_distribution_reference = apply_current_pose(mass_distribution_reference_standard_pose, current_pose, N);
            mass_distribution_noisy = apply_current_pose(mass_distribution_noisy_standard_pose, current_pose, N);
            
            clf
            global_frame_plotter;
            surf(X, Y, Z, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', '0.1');
            hold on;
            technical_frame_plotter(output_containers.(substruct_str.T_reference).SVDLS(k, :, 1), unit_vector_local.REF, "REF");
            technical_frame_plotter(output_containers.(substruct_str.T_estimated).SVDLS(k, :, 1), unit_vector_local.SVDLS, "SVDLS");
            technical_frame_plotter(output_containers.(substruct_str.T_estimated).PTUR(k, :, 1), unit_vector_local.PTUR, "PTUR");
            for j=1:N
                scatter3(mass_distribution_reference(j).X, mass_distribution_reference(j).Y, mass_distribution_reference(j).Z, 'k', 'filled');
                scatter3(mass_distribution_noisy(j).X, mass_distribution_noisy(j).Y, mass_distribution_noisy(j).Z, 'm');
            end
            pretty_plotter(fig_scale);
           
            if movingView_flag
                new_azimuth = original_azimuth + (k-1)*dTheta;
                fig_scale.view_vector = [new_azimuth original_elevation];
            end
           
            current_frame = getframe(figh);
            writeVideo(my_writer, current_frame);
        end
        close(my_writer);
        close(figh);

    end

end

function global_frame_plotter
    axis_size = 10;
    i = [axis_size;0;0];
    j = [0;axis_size;0];
    k = [0;0;axis_size];
    plot3([0 i(1)],[0 i(2)],[0 i(3)], 'r', 'LineWidth',2)
    hold on
    plot3([0 j(1)],[0 j(2)],[0 j(3)], 'g', 'LineWidth',2)
    plot3([0 k(1)],[0 k(2)],[0 k(3)], 'b', 'LineWidth',2)
end

function technical_frame_plotter(origin_local, unit_vector_local, algorithm_str)
    tf_axis_size = 50;
    T = origin_local;
    
    i_hat = tf_axis_size*unit_vector_local.x;
    j_hat = tf_axis_size*unit_vector_local.y;
    k_hat = tf_axis_size*unit_vector_local.z;
    
    if algorithm_str == "SVDLS"
        plot3([T(1) T(1)+i_hat(1)],[T(2) T(2)+i_hat(2)],[T(3) T(3)+i_hat(3)], 'r--', 'LineWidth', 1)
        plot3([T(1) T(1)+j_hat(1)],[T(2) T(2)+j_hat(2)],[T(3) T(3)+j_hat(3)], 'g--', 'LineWidth', 1)
        plot3([T(1) T(1)+k_hat(1)],[T(2) T(2)+k_hat(2)],[T(3) T(3)+k_hat(3)], 'b--', 'LineWidth', 1)
    elseif algorithm_str == "PTUR"
        plot3([T(1) T(1)+i_hat(1)],[T(2) T(2)+i_hat(2)],[T(3) T(3)+i_hat(3)], 'r', 'LineWidth', 1)
        plot3([T(1) T(1)+j_hat(1)],[T(2) T(2)+j_hat(2)],[T(3) T(3)+j_hat(3)], 'g', 'LineWidth', 1)
        plot3([T(1) T(1)+k_hat(1)],[T(2) T(2)+k_hat(2)],[T(3) T(3)+k_hat(3)], 'b', 'LineWidth', 1)
    elseif algorithm_str == "REF"
        plot3([T(1) T(1)+i_hat(1)],[T(2) T(2)+i_hat(2)],[T(3) T(3)+i_hat(3)], 'r:', 'LineWidth', 1)
        plot3([T(1) T(1)+j_hat(1)],[T(2) T(2)+j_hat(2)],[T(3) T(3)+j_hat(3)], 'g:', 'LineWidth', 1)
        plot3([T(1) T(1)+k_hat(1)],[T(2) T(2)+k_hat(2)],[T(3) T(3)+k_hat(3)], 'b:', 'LineWidth', 1)
    end
end

function unit_vectors_local = generate_unit_vectors_local(R_input)
    R = reshape(R_input, [3,3]);
    
    unit_vectors_local = struct;
    unit_vectors_local.x = R(:,1);
    unit_vectors_local.y = R(:,2);
    unit_vectors_local.z = R(:,3);
end

function fig_scale_info = generate_fig_scale_info_PCT(fig_scale_flag, input_struct, output_containers)
    if fig_scale_flag
        plot_edge_width = 50;
        view_vector = [-37.5 30];

        input_size = input_struct.input_size;
        num_measurements = input_size.NumFrames;
        cylinder_x_mean_list = nan(num_measurements,1);
        cylinder_y_mean_list = nan(num_measurements,1);
        cylinder_z_mean_list = nan(num_measurements,1);
        for n=1:num_measurements
            cylinder_x_mean_list(n) = output_containers.T_estimated.SVDLS(n,1);
            cylinder_y_mean_list(n) = output_containers.T_estimated.SVDLS(n,2);
            cylinder_z_mean_list(n) = output_containers.T_estimated.SVDLS(n,3);

        end

        x_max_limit = max(cylinder_x_mean_list) + plot_edge_width;
        y_max_limit = max(cylinder_y_mean_list) + plot_edge_width;
        z_max_limit = max(cylinder_z_mean_list) + plot_edge_width;

        x_min_limit = min(cylinder_x_mean_list) - plot_edge_width;
        y_min_limit = min(cylinder_y_mean_list) - plot_edge_width;
        z_min_limit = min(cylinder_z_mean_list) - plot_edge_width;

        fig_scale_info = struct;

        fig_scale_info.x_max = x_max_limit;
        fig_scale_info.y_max = y_max_limit;
        fig_scale_info.z_max = z_max_limit;

        fig_scale_info.x_min = x_min_limit;
        fig_scale_info.y_min = y_min_limit;
        fig_scale_info.z_min = z_min_limit;

        fig_scale_info.view_vector = view_vector;
    end
end

function pretty_plotter(fig_scale)
    grid

    description_axes = zeros(11, 1);
    description_axes(1) = plot(NaN,NaN,'r--');
    description_axes(2) = plot(NaN,NaN,'g--');
    description_axes(3) = plot(NaN,NaN,'b--');
    description_axes(4) = plot(NaN,NaN,'r');
    description_axes(5) = plot(NaN,NaN,'g');
    description_axes(6) = plot(NaN,NaN,'b');
    description_axes(7) = plot(NaN,NaN,'r:');
    description_axes(8) = plot(NaN,NaN,'g:');
    description_axes(9) = plot(NaN,NaN,'b:');
    description_axes(10) = plot(NaN,NaN,'m');
    description_axes(11) = plot(NaN,NaN,'black');
    legend(description_axes, 'x: red dashed (SVDLS)','y: green dashed (SVDLS)','z: blue dashed (SVDLS)', 'x: red (PTUR)','y: green (PTUR)','z: blue (PTUR)', 'x: red dotted (REF)','y: green dotted (REF)','z: blue dotted (REF)', 'noisy markers: hollow magenta', 'reference markers: solid black');
    if ~isnan(fig_scale.x_min) && ~isnan(fig_scale.x_max)
        xlim([fig_scale.x_min, fig_scale.x_max])
    end
    if ~isnan(fig_scale.y_min) && ~isnan(fig_scale.y_max)
        ylim([fig_scale.y_min, fig_scale.y_max])
    end
    if ~isnan(fig_scale.z_min) && ~isnan(fig_scale.z_max)
        zlim([fig_scale.z_min, fig_scale.z_max])
    end
    
    view(fig_scale.view_vector);
    
end

function [T,R,unit_vectors] = pose_estimate_femoral_frame(current_reconstructed_position_val)
    
    current_side = "left";

    r_lat_femoral_epicondyle = current_reconstructed_position_val.(current_side + "_lat_femoral_epicondyle");
    r_med_femoral_epicondyle = current_reconstructed_position_val.(current_side + "_med_femoral_epicondyle");
    r_femoral_head = current_reconstructed_position_val.(current_side + "_femoral_head");
    r_origin = 0.5*(r_lat_femoral_epicondyle + r_med_femoral_epicondyle);
    
    ptCloud = pointCloud(transpose(horzcat(r_lat_femoral_epicondyle, r_med_femoral_epicondyle, r_femoral_head)));
    pts_collinear = true;
    if ~any(any(isnan(ptCloud.Location)))
        [~,plane_inlierIndices,~] = pcfitplane(ptCloud,1);
        if ~isempty(plane_inlierIndices)
            pts_collinear = false;
        end
    end

    i = [1;0;0];
    j = [0;1;0];
    k = [0;0;1];

    T = r_origin;
    R = nan(3,3);
    R_calc = nan(3,3);
    
    i_hat = nan(3,1);
    j_hat = nan(3,1);
    k_hat = nan(3,1);

    if ~pts_collinear
        
        ry = r_femoral_head - r_origin;
        
        r_plane = r_lat_femoral_epicondyle - r_med_femoral_epicondyle;
        r_normal = cross(r_plane, ry);
        
        syms u_x u_y u_z;
        eqn_1 = r_normal(1)*u_x + r_normal(2)*u_y + r_normal(3)*u_z == 0;
        eqn_2 = ry(1)*u_x + ry(2)*u_y + ry(3)*u_z == 0;
        eqn_3 = u_x^2 + u_y^2 + u_z^2 == 1;
        [sol_x, sol_y, sol_z] = solve(eqn_1, eqn_2, eqn_3);
        rz_1 = [double(sol_x(1,1)); double(sol_y(1,1)); double(sol_z(1,1))];
        rz_2 = [double(sol_x(2,1)); double(sol_y(2,1)); double(sol_z(2,1))];
        
        if dot(rz_1, r_plane) < 0
            rz_minus = rz_1;
            rz_plus = rz_2;
        else
            rz_minus = rz_2;
            rz_plus = rz_1;
        end
        
        if current_side == "left"
            rz = rz_minus;
        else
            rz = rz_plus;
        end
        
        rx = cross(ry, rz);
        
        i_hat_dir = rx;
        j_hat_dir = ry;
        k_hat_dir = rz;
        
        i_hat = i_hat_dir/norm(i_hat_dir);
        j_hat = j_hat_dir/norm(j_hat_dir);
        k_hat = k_hat_dir/norm(k_hat_dir);

        R_calc = [dot(i_hat,i) dot(j_hat,i) dot(k_hat,i)
             dot(i_hat,j) dot(j_hat,j) dot(k_hat,j)
             dot(i_hat,k) dot(j_hat,k) dot(k_hat,k)];
    end
    
    tol = 1.0e-10;
    if ~any(any(isnan(R_calc)))
        same_inverse_transpose_bool = all(all(abs(inv(R_calc) - R_calc') < tol ));
        det_R_bool = abs(det(R_calc)-1.0) < tol;
        rotation_matrix_properties_satisfied_bool = same_inverse_transpose_bool && det_R_bool;
        if rotation_matrix_properties_satisfied_bool
            R = R_calc;
        end
    end
    
    unit_vectors.x = i_hat;
    unit_vectors.y = j_hat;
    unit_vectors.z = k_hat;
end