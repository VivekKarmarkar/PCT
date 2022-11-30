display_flag = struct;
display_flag.legend_bool = false;
display_flag.title_bool = false;

color_info = assign_color_info(false);
display_flag.color_info = color_info;

N = 4;
filename_thigh = "STA_output_thigh_updated_PCT_" + string(N) + ".mat";
filename_shank = "STA_output_shank_updated_PCT_" + string(N) + ".mat";
load(filename_thigh);
load(filename_shank);
NumConfigurations = length(final_output_shank);

GTD_filename_thigh = "ground_truth_data_thigh_updated.mat";
load(GTD_filename_thigh);
GTD_filename_shank = "ground_truth_data_shank_updated.mat";
load(GTD_filename_shank);

STA_filename_thigh = "box_markers_gait_" + string(N) + ".mat";
load(STA_filename_thigh);
STA_filename_shank = "box_markers_gait_shank_" + string(N) + ".mat";
load(STA_filename_shank);

input_size = generate_input_size_struct(N, NumConfigurations, box_markers_gait_shank);
input_settings = generate_input_settings;

cylinder_data_thigh = generate_cylinder_data_thigh;
cylinder_data_shank = generate_cylinder_data_shank;

kneeAngle_output = struct;
AFoffset_output = struct;

discontinuity_flag_PCT = struct;
discontinuity_flag_PCT.thigh = false(NumConfigurations,1);
discontinuity_flag_PCT.shank = false(NumConfigurations,1);
for j=1:NumConfigurations
    discontinuity_flag_PCT.thigh(j) = final_output(j).discontinuity_flag_PCT;
    discontinuity_flag_PCT.shank(j) = final_output_shank(j).discontinuity_flag_PCT;
end
discontinuity_list_PCT = struct;
discontinuity_list_PCT.thigh = find(~discontinuity_flag_PCT.thigh);
discontinuity_list_PCT.shank = find(~discontinuity_flag_PCT.shank);

continuity_list_PCT = struct;
if length(discontinuity_list_PCT.shank) < length(discontinuity_list_PCT.thigh)
    continuity_list_PCT.thigh = setdiff(1:NumConfigurations, discontinuity_list_PCT.thigh(1:length(discontinuity_list_PCT.shank)));
    continuity_list_PCT.shank = setdiff(1:NumConfigurations, discontinuity_list_PCT.shank);
    input_size.NumPCT = length(discontinuity_list_PCT.shank);
else
    continuity_list_PCT.thigh = setdiff(1:NumConfigurations, discontinuity_list_PCT.thigh);
    continuity_list_PCT.shank = setdiff(1:NumConfigurations, discontinuity_list_PCT.shank(1:length(discontinuity_list_PCT.thigh)));
    input_size.NumPCT = length(discontinuity_list_PCT.thigh);
end

for j=1:NumConfigurations
    configuration_idx = j;
    
    folderName_settings_thigh = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\MarkerConfigurationConstrainedSettingsUpdated_" + string(N);
    fileName_configuration_thigh = folderName_settings_thigh + "\FixedMarkerDistribution_" + string(configuration_idx);
    load(fileName_configuration_thigh);
    
    folderName_settings_shank = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\MarkerConfigurationConstrainedSettingsShank_" + string(N);
    fileName_configuration_shank = folderName_settings_shank + "\FixedMarkerDistribution_" + string(configuration_idx);
    mass_distribution_reference_shank = load(fileName_configuration_shank);
    
    input_data_thigh = struct;
    input_data_thigh.cylinder_data_thigh = cylinder_data_thigh;
    input_data_thigh.configuration_idx = configuration_idx;
    input_data_thigh.mass_distribution_reference_thigh = mass_distribution_reference;
    input_data_thigh.box_markers_gait_thigh = box_markers_gait;
    input_data_thigh.pose_gait_thigh = pose_gait;
    
    input_data_shank = struct;
    input_data_shank.cylinder_data_shank = cylinder_data_shank;
    input_data_shank.configuration_idx = configuration_idx;
    input_data_shank.mass_distribution_reference_shank = mass_distribution_reference_shank.mass_distribution_reference;
    input_data_shank.box_markers_gait_shank = box_markers_gait_shank;
    input_data_shank.pose_gait_shank = pose_gait_shank;

    input_struct = struct;
    input_struct.input_size = input_size;
    input_struct.input_settings = input_settings;
    input_struct.input_data_thigh = input_data_thigh;
    input_struct.input_data_shank = input_data_shank;
    input_struct.color_info = color_info;
    
    output_container_current_PCT = struct;
    knee_angle_current = struct;
    knee_angle_current.estimated = struct;
    
    if j<=input_size.NumPCT
        current_PCT_thigh_idx = discontinuity_list_PCT.thigh(j);
        current_PCT_shank_idx = discontinuity_list_PCT.shank(j);
        output_container_current_PCT.thigh = final_output(current_PCT_thigh_idx).output_containers;
        output_container_current_PCT.shank = final_output_shank(current_PCT_shank_idx).output_containers;
    else
        current_PCT_thigh_idx = continuity_list_PCT.thigh(j-input_size.NumPCT);
        current_PCT_shank_idx = continuity_list_PCT.shank(j-input_size.NumPCT);
        output_container_current_PCT.thigh = final_output(current_PCT_thigh_idx).output_containers;
        output_container_current_PCT.shank = final_output_shank(current_PCT_shank_idx).output_containers;
    end
    
    knee_angle_current.reference = estimate_knee_angle_gs(output_container_current_PCT, input_size, "reference", "SVDLS");
    knee_angle_current.estimated.PCTO = estimate_knee_angle_gs(output_container_current_PCT, input_size, "estimated", "PCTO");
    knee_angle_current.estimated.SVDLS = estimate_knee_angle_gs(output_container_current_PCT, input_size, "estimated", "SVDLS");
    knee_angle_current.estimated.PTUR = estimate_knee_angle_gs(output_container_current_PCT, input_size, "estimated", "PTUR");
  
    generate_output_plots_thigh(input_struct, output_container_current_PCT.thigh);
    generate_output_plots_shank(input_struct, output_container_current_PCT.shank);
   
    generate_knee_angle_plot(knee_angle_current, input_struct);
    AFoffset_current = generate_AL_origin_offset_plot(output_container_current_PCT, input_struct);
    
    current_fig_scale_flag = true;
    current_animation_flag = false;
    fig_scale = generate_fig_scale_info_PCT(current_fig_scale_flag, input_struct, output_container_current_PCT);
    input_struct.visualization_settings = struct;
    input_struct.visualization_settings.fig_scale = fig_scale;
    input_struct.visualization_settings.animation_flag = current_animation_flag;
    input_struct.visualization_settings.movingView_flag_val = false;
    create_animation_RF(input_struct, output_container_current_PCT, "TF");
    create_animation_RF(input_struct, output_container_current_PCT, "AF");
    
    kneeAngle_output(j).REF = knee_angle_current.reference;
    kneeAngle_output(j).SVDLS = knee_angle_current.estimated.SVDLS;
    kneeAngle_output(j).PTUR = knee_angle_current.estimated.PTUR;
    kneeAngle_output(j).PCTO = knee_angle_current.estimated.PCTO;
    
    AFoffset_output(j).thigh = AFoffset_current.thigh;
    AFoffset_output(j).shank = AFoffset_current.shank;
    
    disp(j);
end

final_output_thigh_mean = generate_ouput_stats_thigh(final_output, input_size, display_flag);
final_output_shank_mean = generate_ouput_stats_shank(final_output_shank, input_size, display_flag);
final_output_mean = generate_ouput_stats(AFoffset_output, kneeAngle_output, input_size, NumConfigurations, display_flag);

plot_PCT_discontinuous_markerRatio(input_size, final_output, display_flag);
plot_PCT_discontinuous_epsilon(input_size, final_output, display_flag);

function color_info = assign_color_info(default_style_bool)
    color_info = struct;
    color_info.reference = 'k';
    
    if default_style_bool
        color_info.SVDLS = [1 0 1];
        color_info.PTUR = [0 0 1];
        color_info.PCTO = [1 0 0];
    else
        color_info.SVDLS = [6/255 199/255 151/255];
        color_info.PTUR = [0 114/255 189/255];
        color_info.PCTO = [162/255 20/255 47/255];
    end
end

function plot_PCT_discontinuous_epsilon(input_size, final_output, display_flag)
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    N = input_size.N;
    
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N);
    path_PCTdiscontinuous = path_img + "\PCTdiscontinuous_plots\";
    
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\PCT paper drafts\Images\";
    path_PCTdiscontinuous_final = path_img_final;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(3,1, 'TileSpacing', 'compact');
    
    nexttile;
    plot(fraction_completed_gait_cycle, final_output(3).output_containers.eps_star.PCTO, 'k', 'LineWidth', 1, 'DisplayName', 'Epsilon');
    hold on;
    yL = ylim;
    fill([0.17 0.22 0.22 0.17], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    fill([0.37 0.57 0.57 0.37], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    fill([0.67 0.75 0.75 0.67], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    title('Optimization parameter (epsilon)');
    ax = gca;
    ax.FontSize = 20;
    
    nexttile;
    plot(fraction_completed_gait_cycle, final_output(3).output_containers.mass_redistribution.PCTO(:,1), 'Color', [162/255 20/255 47/255], 'LineWidth', 1, 'DisplayName', 'm_1');
    hold on;
    plot(fraction_completed_gait_cycle, final_output(3).output_containers.mass_redistribution.PCTO(:,2), 'Color', [0 114/255 189/255], 'LineWidth', 2, 'DisplayName', 'm_2');
    plot(fraction_completed_gait_cycle, final_output(3).output_containers.mass_redistribution.PCTO(:,3), 'Color', [6/255 199/255 151/255], 'LineWidth', 1, 'DisplayName', 'm_3');
    plot(fraction_completed_gait_cycle, final_output(3).output_containers.mass_redistribution.PCTO(:,4), 'k', 'LineWidth', 1, 'DisplayName', 'm_4');
    title('Marker mass');
    ylim([0,2]);
    yL = ylim;
    fill([0.17 0.22 0.22 0.17], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    fill([0.37 0.57 0.57 0.37], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    fill([0.67 0.75 0.75 0.67], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ax = gca;
    ax.FontSize = 20;
    
    nexttile;
    plot(input_size.fraction_completed_gait_cycle, final_output(3).output_containers.CM_offset.PCTO, 'k', 'LineWidth', 1, 'DisplayName', 'CM');
    hold on;
    yL = ylim;
    fill([0.17 0.22 0.22 0.17], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    fill([0.37 0.57 0.57 0.37], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    fill([0.67 0.75 0.75 0.67], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    xlabel('Fraction completed gait cycle');
    ax = gca;
    ax.FontSize = 20;
    title('CM offset (cm)');
    if display_flag.legend_bool
        legend('Location', 'best');
    end
    grid;
    if display_flag.title_bool
        title('PCT discontinuity');
    end
    
    filename_PCTdiscontinuous = path_PCTdiscontinuous + "PCTdiscontinuous_Epsilon_plot.png";
    exportgraphics(gcf, filename_PCTdiscontinuous);
    filename_PCTdiscontinuous_final = path_PCTdiscontinuous_final + "PCTdiscontinuous_Epsilon.png";
    exportgraphics(gcf, filename_PCTdiscontinuous_final);
    close(gcf);
end

function plot_PCT_discontinuous_markerRatio(input_size, final_output, display_flag)
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    N = input_size.N;
    
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N);
    path_PCTdiscontinuous = path_img + "\PCTdiscontinuous_plots\";
    
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\PCT paper drafts\Images\";
    path_PCTdiscontinuous_final = path_img_final;
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    plot(fraction_completed_gait_cycle, final_output(5).output_containers.eps_star.PCTO, 'k', 'LineWidth', 1, 'DisplayName', 'Epsilon');
    hold on;
    scatter(fraction_completed_gait_cycle, final_output(5).output_containers.eps_star.PCTO - 1, [], [0 114/255 189/255], 'LineWidth', 0.1, 'DisplayName', 'Epsilon - 1');
    plot(fraction_completed_gait_cycle, final_output(5).output_containers.mass_redistribution.PCTO(:,4), 'Color', [162/255 20/255 47/255], 'LineWidth', 1, 'DisplayName', 'm_4');
    plot(fraction_completed_gait_cycle, final_output(5).output_containers.mass_redistribution.PCTO(:,2), 'Color', [6/255 199/255 151/255], 'LineWidth', 1, 'DisplayName', 'm_2');
    plot(input_size.fraction_completed_gait_cycle, final_output(5).output_containers.CM_offset.PCTO, 'Color', [237/255 241/255 19/255], 'LineWidth', 1.5, 'DisplayName', 'CM');
    yL = ylim;
    fill([0.35 0.45 0.45 0.35], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    if display_flag.legend_bool
        legend('Location', 'best');
    end
    grid;
    xlabel('Fraction completed gait cycle');
    if display_flag.title_bool
        title('PCT discontinuity');
    end
    ax = gca;
    ax.FontSize = 20;
    
    filename_PCTdiscontinuous = path_PCTdiscontinuous + "PCTdiscontinuous_MarkerSwitching_plot.png";
    exportgraphics(gcf, filename_PCTdiscontinuous);
    filename_PCTdiscontinuous_final = path_PCTdiscontinuous_final + "PCTdiscontinuous_MarkerSwitching.png";
    exportgraphics(gcf, filename_PCTdiscontinuous_final);
    close(gcf);
end

function final_output_mean = generate_ouput_stats(AFoffset_output, kneeAngle_output, input_size, NumConfigurations, display_flag)

    final_output_mean = struct;
    
    color_info = display_flag.color_info;

    kneeAngle_str_list = ["flexion", "adduction" , "externalRotation"];
    
    NumKneeAngles = length(kneeAngle_str_list);
    
    final_output_mean.kneeAngle = struct;
    final_output_mean.kneeAngle.REF = nan(input_size.NumFrames, NumKneeAngles);
    final_output_mean.kneeAngle.SVDLS = nan(input_size.NumFrames, NumKneeAngles);
    final_output_mean.kneeAngle.PTUR = nan(input_size.NumFrames, NumKneeAngles);
    final_output_mean.kneeAngle.PCTO = nan(input_size.NumFrames, NumKneeAngles);
    
    final_output_std.kneeAngle = struct;
    final_output_std.kneeAngle.REF = nan(input_size.NumFrames, NumKneeAngles);
    final_output_std.kneeAngle.SVDLS = nan(input_size.NumFrames, NumKneeAngles);
    final_output_std.kneeAngle.PTUR = nan(input_size.NumFrames, NumKneeAngles);
    final_output_std.kneeAngle.PCTO = nan(input_size.NumFrames, NumKneeAngles);
    
    final_output_ub.kneeAngle = struct;
    final_output_ub.kneeAngle.REF = nan(input_size.NumFrames, NumKneeAngles);
    final_output_ub.kneeAngle.SVDLS = nan(input_size.NumFrames, NumKneeAngles);
    final_output_ub.kneeAngle.PTUR = nan(input_size.NumFrames, NumKneeAngles);
    final_output_ub.kneeAngle.PCTO = nan(input_size.NumFrames, NumKneeAngles);
    
    final_output_lb.kneeAngle = struct;
    final_output_lb.kneeAngle.REF = nan(input_size.NumFrames, NumKneeAngles);
    final_output_lb.kneeAngle.SVDLS = nan(input_size.NumFrames, NumKneeAngles);
    final_output_lb.kneeAngle.PTUR = nan(input_size.NumFrames, NumKneeAngles);
    final_output_lb.kneeAngle.PCTO = nan(input_size.NumFrames, NumKneeAngles);
    
    final_output_mean.AFoffset = struct;
    
    final_output_mean.AFoffset.thigh = struct;
    final_output_mean.AFoffset.thigh.SVDLS = nan(input_size.NumFrames, 1);
    final_output_mean.AFoffset.thigh.PTUR = nan(input_size.NumFrames, 1);
    final_output_mean.AFoffset.thigh.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_std.AFoffset.thigh = struct;
    final_output_std.AFoffset.thigh.SVDLS = nan(input_size.NumFrames, 1);
    final_output_std.AFoffset.thigh.PTUR = nan(input_size.NumFrames, 1);
    final_output_std.AFoffset.thigh.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_ub.AFoffset.thigh = struct;
    final_output_ub.AFoffset.thigh.SVDLS = nan(input_size.NumFrames, 1);
    final_output_ub.AFoffset.thigh.PTUR = nan(input_size.NumFrames, 1);
    final_output_ub.AFoffset.thigh.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_lb.AFoffset.thigh = struct;
    final_output_lb.AFoffset.thigh.SVDLS = nan(input_size.NumFrames, 1);
    final_output_lb.AFoffset.thigh.PTUR = nan(input_size.NumFrames, 1);
    final_output_lb.AFoffset.thigh.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_mean.AFoffset.shank = struct;
    final_output_mean.AFoffset.shank.SVDLS = nan(input_size.NumFrames, 1);
    final_output_mean.AFoffset.shank.PTUR = nan(input_size.NumFrames, 1);
    final_output_mean.AFoffset.shank.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_std.AFoffset.shank = struct;
    final_output_std.AFoffset.shank.SVDLS = nan(input_size.NumFrames, 1);
    final_output_std.AFoffset.shank.PTUR = nan(input_size.NumFrames, 1);
    final_output_std.AFoffset.shank.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_ub.AFoffset.shank = struct;
    final_output_ub.AFoffset.shank.SVDLS = nan(input_size.NumFrames, 1);
    final_output_ub.AFoffset.shank.PTUR = nan(input_size.NumFrames, 1);
    final_output_ub.AFoffset.shank.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_lb.AFoffset.shank = struct;
    final_output_lb.AFoffset.shank.SVDLS = nan(input_size.NumFrames, 1);
    final_output_lb.AFoffset.shank.PTUR = nan(input_size.NumFrames, 1);
    final_output_lb.AFoffset.shank.PCTO = nan(input_size.NumFrames, 1);
    
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    x2 = vertcat(fraction_completed_gait_cycle', flipud(fraction_completed_gait_cycle'));
    ToeOff = 0.6;

    for t=1:input_size.NumFrames
        disp(t);

        kneeAngle_current = struct;
        kneeAngle_current.REF = nan(NumConfigurations, NumKneeAngles);
        kneeAngle_current.SVDLS = nan(NumConfigurations, NumKneeAngles);
        kneeAngle_current.PTUR = nan(NumConfigurations, NumKneeAngles);
        kneeAngle_current.PCTO = nan(NumConfigurations, NumKneeAngles);
        
        AFoffset_current = struct;
        
        AFoffset_current.thigh = struct;
        AFoffset_current.thigh.SVDLS = nan(NumConfigurations, 1);
        AFoffset_current.thigh.PTUR = nan(NumConfigurations, 1);
        AFoffset_current.thigh.PCTO = nan(NumConfigurations, 1);
        
        AFoffset_current.shank = struct;
        AFoffset_current.shank.SVDLS = nan(NumConfigurations, 1);
        AFoffset_current.shank.PTUR = nan(NumConfigurations, 1);
        AFoffset_current.shank.PCTO = nan(NumConfigurations, 1);

        for k=1:NumConfigurations
            disp(k);
            
            if k <= input_size.NumPCT
                for n=1:NumKneeAngles
                    kneeAngle_current_data.PCTO = kneeAngle_output(k).PCTO.(kneeAngle_str_list(n));
                    kneeAngle_current.PCTO(k,n) = kneeAngle_current_data.PCTO(t);
                end
                AFoffset_current.thigh.PCTO(k) = AFoffset_output(k).thigh.PCTO(t);
                AFoffset_current.shank.PCTO(k) = AFoffset_output(k).shank.PCTO(t);
            end
            
            for n=1:NumKneeAngles
                
                kneeAngle_current_data = struct;
                
                kneeAngle_current_data.REF = kneeAngle_output(k).REF.(kneeAngle_str_list(n));
                kneeAngle_current_data.SVDLS = kneeAngle_output(k).SVDLS.(kneeAngle_str_list(n));
                kneeAngle_current_data.PTUR = kneeAngle_output(k).PTUR.(kneeAngle_str_list(n));
                
                kneeAngle_current.REF(k,n) = kneeAngle_current_data.REF(t);
                kneeAngle_current.SVDLS(k,n) = kneeAngle_current_data.SVDLS(t);
                kneeAngle_current.PTUR(k,n) = kneeAngle_current_data.PTUR(t);
                
            end
            
            AFoffset_current.thigh.SVDLS(k) = AFoffset_output(k).thigh.SVDLS(t);
            AFoffset_current.thigh.PTUR(k) = AFoffset_output(k).thigh.PTUR(t);
            
            AFoffset_current.shank.SVDLS(k) = AFoffset_output(k).shank.SVDLS(t);
            AFoffset_current.shank.PTUR(k) = AFoffset_output(k).shank.PTUR(t);
            
        end

        for n=1:NumKneeAngles
            final_output_mean.kneeAngle.REF(t,n) = nanmean(kneeAngle_current.REF(:,n));
            final_output_mean.kneeAngle.SVDLS(t,n) = nanmean(kneeAngle_current.SVDLS(:,n));
            final_output_mean.kneeAngle.PTUR(t,n) = nanmean(kneeAngle_current.PTUR(:,n));
            final_output_mean.kneeAngle.PCTO(t,n) = nanmean(kneeAngle_current.PCTO(:,n));
            
            final_output_std.kneeAngle.REF(t,n) = std(kneeAngle_current.REF(:,n), 'omitnan');
            final_output_std.kneeAngle.SVDLS(t,n) = std(kneeAngle_current.SVDLS(:,n), 'omitnan');
            final_output_std.kneeAngle.PTUR(t,n) = std(kneeAngle_current.PTUR(:,n), 'omitnan');
            final_output_std.kneeAngle.PCTO(t,n) = std(kneeAngle_current.PCTO(:,n), 'omitnan');
            
            final_output_ub.kneeAngle.REF(t,n) = final_output_mean.kneeAngle.REF(t,n) + final_output_std.kneeAngle.REF(t,n);
            final_output_ub.kneeAngle.SVDLS(t,n) = final_output_mean.kneeAngle.SVDLS(t,n) + final_output_std.kneeAngle.SVDLS(t,n);
            final_output_ub.kneeAngle.PTUR(t,n) = final_output_mean.kneeAngle.PTUR(t,n) + final_output_std.kneeAngle.PTUR(t,n);
            final_output_ub.kneeAngle.PCTO(t,n) = final_output_mean.kneeAngle.PCTO(t,n) + final_output_std.kneeAngle.PCTO(t,n);
            
            final_output_lb.kneeAngle.REF(t,n) = final_output_mean.kneeAngle.REF(t,n) - final_output_std.kneeAngle.REF(t,n);
            final_output_lb.kneeAngle.SVDLS(t,n) = final_output_mean.kneeAngle.SVDLS(t,n) - final_output_std.kneeAngle.SVDLS(t,n);
            final_output_lb.kneeAngle.PTUR(t,n) = final_output_mean.kneeAngle.PTUR(t,n) - final_output_std.kneeAngle.PTUR(t,n);
            final_output_lb.kneeAngle.PCTO(t,n) = final_output_mean.kneeAngle.PCTO(t,n) - final_output_std.kneeAngle.PCTO(t,n);
            
        end
        
        final_output_mean.AFoffset.thigh.SVDLS(t) = nanmean(AFoffset_current.thigh.SVDLS);
        final_output_mean.AFoffset.thigh.PTUR(t) = nanmean(AFoffset_current.thigh.PTUR);
        final_output_mean.AFoffset.thigh.PCTO(t) = nanmean(AFoffset_current.thigh.PCTO);
        
        final_output_std.AFoffset.thigh.SVDLS(t) = std(AFoffset_current.thigh.SVDLS, 'omitnan');
        final_output_std.AFoffset.thigh.PTUR(t) = std(AFoffset_current.thigh.PTUR, 'omitnan');
        final_output_std.AFoffset.thigh.PCTO(t) = std(AFoffset_current.thigh.PCTO, 'omitnan');
        
        final_output_ub.AFoffset.thigh.SVDLS(t) = final_output_mean.AFoffset.thigh.SVDLS(t) + final_output_std.AFoffset.thigh.SVDLS(t);
        final_output_ub.AFoffset.thigh.PTUR(t) = final_output_mean.AFoffset.thigh.PTUR(t) + final_output_std.AFoffset.thigh.PTUR(t);
        final_output_ub.AFoffset.thigh.PCTO(t) = final_output_mean.AFoffset.thigh.PCTO(t) + final_output_std.AFoffset.thigh.PCTO(t);
        
        final_output_lb.AFoffset.thigh.SVDLS(t) = final_output_mean.AFoffset.thigh.SVDLS(t) - final_output_std.AFoffset.thigh.SVDLS(t);
        final_output_lb.AFoffset.thigh.PTUR(t) = final_output_mean.AFoffset.thigh.PTUR(t) - final_output_std.AFoffset.thigh.PTUR(t);
        final_output_lb.AFoffset.thigh.PCTO(t) = final_output_mean.AFoffset.thigh.PCTO(t) - final_output_std.AFoffset.thigh.PCTO(t);
        
        final_output_mean.AFoffset.shank.SVDLS(t) = nanmean(AFoffset_current.shank.SVDLS);
        final_output_mean.AFoffset.shank.PTUR(t) = nanmean(AFoffset_current.shank.PTUR);
        final_output_mean.AFoffset.shank.PCTO(t) = nanmean(AFoffset_current.shank.PCTO);
        
        final_output_std.AFoffset.shank.SVDLS(t) = std(AFoffset_current.shank.SVDLS, 'omitnan');
        final_output_std.AFoffset.shank.PTUR(t) = std(AFoffset_current.shank.PTUR, 'omitnan');
        final_output_std.AFoffset.shank.PCTO(t) = std(AFoffset_current.shank.PCTO, 'omitnan');
        
        final_output_ub.AFoffset.shank.SVDLS(t) = final_output_mean.AFoffset.shank.SVDLS(t) + final_output_std.AFoffset.shank.SVDLS(t);
        final_output_ub.AFoffset.shank.PTUR(t) = final_output_mean.AFoffset.shank.PTUR(t) + final_output_std.AFoffset.shank.PTUR(t);
        final_output_ub.AFoffset.shank.PCTO(t) = final_output_mean.AFoffset.shank.PCTO(t) + final_output_std.AFoffset.shank.PCTO(t);
        
        final_output_lb.AFoffset.shank.SVDLS(t) = final_output_mean.AFoffset.shank.SVDLS(t) - final_output_std.AFoffset.shank.SVDLS(t);
        final_output_lb.AFoffset.shank.PTUR(t) = final_output_mean.AFoffset.shank.PTUR(t) - final_output_std.AFoffset.shank.PTUR(t);
        final_output_lb.AFoffset.shank.PCTO(t) = final_output_mean.AFoffset.shank.PCTO(t) - final_output_std.AFoffset.shank.PCTO(t);
        
    end

    inBetween = struct;
    
    inBetween.kneeAngle = struct;
    
    inBetween.AFoffset = struct;
    
    inBetween.AFoffset.thigh = struct;
    inBetween.AFoffset.thigh.SVDLS = vertcat(final_output_ub.AFoffset.thigh.SVDLS, flipud(final_output_lb.AFoffset.thigh.SVDLS));
    inBetween.AFoffset.thigh.PTUR = vertcat(final_output_ub.AFoffset.thigh.PTUR, flipud(final_output_lb.AFoffset.thigh.PTUR));
    inBetween.AFoffset.thigh.PCTO = vertcat(final_output_ub.AFoffset.thigh.PCTO, flipud(final_output_lb.AFoffset.thigh.PCTO));
    
    inBetween.AFoffset.shank = struct;
    inBetween.AFoffset.shank.SVDLS = vertcat(final_output_ub.AFoffset.shank.SVDLS, flipud(final_output_lb.AFoffset.shank.SVDLS));
    inBetween.AFoffset.shank.PTUR = vertcat(final_output_ub.AFoffset.shank.PTUR, flipud(final_output_lb.AFoffset.shank.PTUR));
    inBetween.AFoffset.shank.PCTO = vertcat(final_output_ub.AFoffset.shank.PCTO, flipud(final_output_lb.AFoffset.shank.PCTO));
    
    N = input_size.N;
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N);
    path_avg = path_img + "\Average_plots\";
    
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\PCT paper drafts\Images\";
    path_avg_final = path_img_final;
    
    yl = [0 -6 -20];
    yu = [60 4 0];
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(3,1, 'TileSpacing', 'compact');
    for n=1:NumKneeAngles
        inBetween.kneeAngle.SVDLS = vertcat(final_output_ub.kneeAngle.SVDLS(:,n), flipud(final_output_lb.kneeAngle.SVDLS(:,n)));
        inBetween.kneeAngle.PTUR = vertcat(final_output_ub.kneeAngle.PTUR(:,n), flipud(final_output_lb.kneeAngle.PTUR(:,n)));
        inBetween.kneeAngle.PCTO = vertcat(final_output_ub.kneeAngle.PCTO(:,n), flipud(final_output_lb.kneeAngle.PCTO(:,n)));
        nexttile;
        xline(ToeOff, 'k', 'LineWidth', 1.5, 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_mean.kneeAngle.REF(:,n), 'Color', color_info.reference, 'LineWidth', 1, 'DisplayName', 'REF');
        plot(fraction_completed_gait_cycle, final_output_mean.kneeAngle.SVDLS(:,n), 'Color', color_info.SVDLS,  'LineWidth', 1, 'DisplayName', 'SVD-LS');
        plot(fraction_completed_gait_cycle, final_output_mean.kneeAngle.PTUR(:,n), 'Color', color_info.PTUR,  'LineWidth', 1, 'DisplayName', 'PCT-PT');
        plot(fraction_completed_gait_cycle, final_output_mean.kneeAngle.PCTO(:,n), 'Color', color_info.PCTO,  'LineWidth', 1, 'DisplayName', 'PCT');
        fill(x2, inBetween.kneeAngle.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
        fill(x2, inBetween.kneeAngle.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
        fill(x2(2:end-1), inBetween.kneeAngle.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
        fill([0 0.6 0.6 0], [yl(n) yl(n) yu(n) yu(n)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        grid;
        if n==1
            if display_flag.legend_bool
                legend('Location', 'best');
            end
        end
        if n==3
            xlabel("Fraction completed gait cycle");
        end
        ylabel('Degrees');
        title(kneeAngle_str_list(n));
        ax = gca;
        ax.FontSize = 20;
    end
    if display_flag.title_bool
        suptitle("Average Knee Angle Plots , PCT configurations = " + string(input_size.NumPCT), 'FontSize', 20);
    end
    filename_KneeAngleAverage = path_avg + "KneeAngleAverage_plot.png";
    exportgraphics(gcf, filename_KneeAngleAverage);
    filename_KneeAngleAverage_final = path_avg_final + "KneeAngle.png";
    exportgraphics(gcf, filename_KneeAngleAverage_final);
    close(gcf);
    
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,1, 'TileSpacing', 'compact');
    
    nexttile;
    xline(ToeOff, 'k', 'LineWidth', 1.5, 'DisplayName', 'Toe Off');
    hold on;
    plot(fraction_completed_gait_cycle, final_output_mean.AFoffset.thigh.SVDLS, 'Color', color_info.SVDLS,  'LineWidth', 1, 'DisplayName', 'SVD-LS');
    plot(fraction_completed_gait_cycle, final_output_mean.AFoffset.thigh.PTUR, 'Color', color_info.PTUR,  'LineWidth', 1, 'DisplayName', 'PCT-PT');
    plot(fraction_completed_gait_cycle, final_output_mean.AFoffset.thigh.PCTO, 'Color', color_info.PCTO,  'LineWidth', 1, 'DisplayName', 'PCT');
    fill(x2, inBetween.AFoffset.thigh.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
    fill(x2, inBetween.AFoffset.thigh.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
    fill(x2(2:end-1), inBetween.AFoffset.thigh.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
    fill([0 0.6 0.6 0], [0 0 4 4], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    grid;
    if display_flag.legend_bool
        legend('Location', 'best');
    end
    ylabel("Error (cm)");
    title("Thigh");
    ax = gca;
    ax.FontSize = 20;
    ylim([0, 4]);
    
    nexttile;
    xline(ToeOff, 'k', 'LineWidth', 1.5, 'DisplayName', 'Toe Off');
    hold on;
    plot(fraction_completed_gait_cycle, final_output_mean.AFoffset.shank.SVDLS, 'Color', color_info.SVDLS,  'LineWidth', 1, 'DisplayName', 'SVD-LS');
    plot(fraction_completed_gait_cycle, final_output_mean.AFoffset.shank.PTUR, 'Color', color_info.PTUR,  'LineWidth', 1, 'DisplayName', 'PTUR');
    plot(fraction_completed_gait_cycle, final_output_mean.AFoffset.shank.PCTO, 'Color', color_info.PCTO,  'LineWidth', 1, 'DisplayName', 'PCTO');
    fill(x2, inBetween.AFoffset.shank.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
    fill(x2, inBetween.AFoffset.shank.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
    fill(x2(2:end-1), inBetween.AFoffset.shank.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
    fill([0 0.6 0.6 0], [0 0 4 4], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    grid;
    xlabel("Fraction completed gait cycle");
    ylabel("Error (cm)");
    title("Shank");
    ax = gca;
    ax.FontSize = 20;
    ylim([0, 4]);
    
    if display_flag.title_bool
        suptitle("Average Anatomical Frame offset Plots , PCT configurations = " + string(input_size.NumPCT), 'FontSize', 20);        
    end
    
    filename_AFOAverage = path_avg + "AFOAverage_plot.png";
    exportgraphics(gcf, filename_AFOAverage);
    filename_AFOAverage_final = path_avg_final + "AFoffset.png";
    exportgraphics(gcf, filename_AFOAverage_final);
    close(gcf);
    
end

function final_output_thigh_mean = generate_ouput_stats_thigh(final_output, input_size, display_flag)

    final_output_thigh_mean = struct;
    
    input_size.Num_AL_Thigh = 4;
    
    color_info = display_flag.color_info;
    
    final_output_thigh_mean.T_estimated = struct;
    final_output_thigh_mean.T_estimated.SVDLS = nan(input_size.NumFrames, 3, 1);
    final_output_thigh_mean.T_estimated.PTUR = nan(input_size.NumFrames, 3, 1);
    final_output_thigh_mean.T_estimated.PCTO = nan(input_size.NumFrames, 3, 1);

    final_output_thigh_mean.CM_offset = struct;
    final_output_thigh_mean.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_thigh_mean.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_thigh_mean.CM_offset.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_std.CM_offset = struct;
    final_output_std.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_std.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_std.CM_offset.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_ub.CM_offset = struct;
    final_output_ub.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_ub.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_ub.CM_offset.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_lb.CM_offset = struct;
    final_output_lb.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_lb.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_lb.CM_offset.PCTO = nan(input_size.NumFrames, 1);

    final_output_thigh_mean.cluster_offset = struct;
    final_output_thigh_mean.cluster_offset.SVDLS = nan(input_size.NumFrames, input_size.N);
    final_output_thigh_mean.cluster_offset.PTUR = nan(input_size.NumFrames, input_size.N);
    final_output_thigh_mean.cluster_offset.PCTO = nan(input_size.NumFrames, input_size.N);

    final_output_thigh_mean.AL_offset = struct;
    final_output_thigh_mean.AL_offset.SVDLS = nan(input_size.NumFrames, input_size.Num_AL_Thigh);
    final_output_thigh_mean.AL_offset.PTUR = nan(input_size.NumFrames, input_size.Num_AL_Thigh);
    final_output_thigh_mean.AL_offset.PCTO = nan(input_size.NumFrames, input_size.Num_AL_Thigh);

    NumEigenvectors = 3;
    final_output_thigh_mean.eigenvector_offset = struct;
    final_output_thigh_mean.eigenvector_offset.SVDLS = nan(input_size.NumFrames, NumEigenvectors);
    final_output_thigh_mean.eigenvector_offset.PTUR = nan(input_size.NumFrames, NumEigenvectors);
    final_output_thigh_mean.eigenvector_offset.PCTO = nan(input_size.NumFrames, NumEigenvectors);
    
    NumConfigurations = min([length(final_output)-1, input_size.NumConfigurations]);
    
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    x2 = vertcat(fraction_completed_gait_cycle', flipud(fraction_completed_gait_cycle'));
    x2_short = vertcat(fraction_completed_gait_cycle(2:end)', flipud(fraction_completed_gait_cycle(2:end)'));
    ToeOff = 0.6;

    for t=1:input_size.NumFrames
        disp(t);
        
        T_estimated_current = struct;
        T_estimated_current.SVDLS = nan(NumConfigurations, 3, 1);
        T_estimated_current.PTUR = nan(NumConfigurations, 3, 1);
        T_estimated_current.PCTO = nan(NumConfigurations, 3, 1);
        
        CM_offset_current = struct;
        CM_offset_current.SVDLS = nan(NumConfigurations, 1);
        CM_offset_current.PTUR = nan(NumConfigurations, 1);
        CM_offset_current.PCTO = nan(NumConfigurations, 1);

        cluster_offset_current = struct;
        cluster_offset_current.SVDLS = nan(NumConfigurations, input_size.N);
        cluster_offset_current.PTUR = nan(NumConfigurations, input_size.N);
        cluster_offset_current.PCTO = nan(NumConfigurations, input_size.N);

        AL_offset_current = struct;
        AL_offset_current.SVDLS = nan(NumConfigurations, input_size.Num_AL_Thigh);
        AL_offset_current.PTUR = nan(NumConfigurations, input_size.Num_AL_Thigh);
        AL_offset_current.PCTO = nan(NumConfigurations, input_size.Num_AL_Thigh);

        eigenvector_offset_current = struct;
        eigenvector_offset_current.SVDLS = nan(NumConfigurations, NumEigenvectors);
        eigenvector_offset_current.PTUR = nan(NumConfigurations, NumEigenvectors);
        eigenvector_offset_current.PCTO = nan(NumConfigurations, NumEigenvectors);

        for k=1:NumConfigurations
            disp(k);
            T_estimated_current.SVDLS(k,:,1) = final_output(k).output_containers.T_estimated.SVDLS(t);
            T_estimated_current.PTUR(k,:,1) = final_output(k).output_containers.T_estimated.PTUR(t);
            T_estimated_current.PCTO(k,:,1) = final_output(k).output_containers.T_estimated.PCTO(t);
            CM_offset_current.SVDLS(k) = final_output(k).output_containers.CM_offset.SVDLS(t);
            CM_offset_current.PTUR(k) = final_output(k).output_containers.CM_offset.PTUR(t);
            
            if ~final_output(k).discontinuity_flag_PCT
                
                CM_offset_current.PCTO(k) = final_output(k).output_containers.CM_offset.PCTO(t);
                
                for n=1:input_size.N
                    cluster_offset_current.PCTO(k,n) = final_output(k).output_containers.cluster_offset.PCTO(t,n);
                end
                
                for n=1:input_size.Num_AL_Thigh
                    AL_offset_current.PCTO(k,n) = final_output(k).output_containers.AL_offset.PCTO(t,n);
                end
                
                for n=1:NumEigenvectors
                    eigenvector_offset_current.PCTO(k,n) = final_output(k).output_containers.eigenvector_offset.PCTO(t,n);
                end
                
            end
            
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
        final_output_thigh_mean.CM_offset.SVDLS(t) = nanmean(CM_offset_current.SVDLS);
        final_output_thigh_mean.CM_offset.PTUR(t) = nanmean(CM_offset_current.PTUR);
        final_output_thigh_mean.CM_offset.PCTO(t) = nanmean(CM_offset_current.PCTO);
        
        final_output_std.CM_offset.SVDLS(t) = std(CM_offset_current.SVDLS, 'omitnan');
        final_output_std.CM_offset.PTUR(t) = std(CM_offset_current.PTUR, 'omitnan');
        final_output_std.CM_offset.PCTO(t) = std(CM_offset_current.PCTO, 'omitnan');
        
        final_output_ub.CM_offset.SVDLS(t) = final_output_thigh_mean.CM_offset.SVDLS(t) + final_output_std.CM_offset.SVDLS(t);
        final_output_ub.CM_offset.PTUR(t) = final_output_thigh_mean.CM_offset.PTUR(t) + final_output_std.CM_offset.PTUR(t);
        final_output_ub.CM_offset.PCTO(t) = final_output_thigh_mean.CM_offset.PCTO(t) + final_output_std.CM_offset.PCTO(t);
        
        final_output_lb.CM_offset.SVDLS(t) = final_output_thigh_mean.CM_offset.SVDLS(t) - final_output_std.CM_offset.SVDLS(t);
        final_output_lb.CM_offset.PTUR(t) = final_output_thigh_mean.CM_offset.PTUR(t) - final_output_std.CM_offset.PTUR(t);
        final_output_lb.CM_offset.PCTO(t) = final_output_thigh_mean.CM_offset.PCTO(t) - final_output_std.CM_offset.PCTO(t);

        for n=1:input_size.N
            
            final_output_thigh_mean.cluster_offset.SVDLS(t,n) = nanmean(cluster_offset_current.SVDLS(:,n));
            final_output_thigh_mean.cluster_offset.PTUR(t,n) = nanmean(cluster_offset_current.PTUR(:,n));
            final_output_thigh_mean.cluster_offset.PCTO(t,n) = nanmean(cluster_offset_current.PCTO(:,n));
            
            final_output_std.cluster_offset.SVDLS(t,n) = std(cluster_offset_current.SVDLS(:,n), 'omitnan');
            final_output_std.cluster_offset.PTUR(t,n) = std(cluster_offset_current.PTUR(:,n), 'omitnan');
            final_output_std.cluster_offset.PCTO(t,n) = std(cluster_offset_current.PCTO(:,n), 'omitnan');
            
            final_output_ub.cluster_offset.SVDLS(t,n) = final_output_thigh_mean.cluster_offset.SVDLS(t,n) + final_output_std.cluster_offset.SVDLS(t,n);
            final_output_ub.cluster_offset.PTUR(t,n) = final_output_thigh_mean.cluster_offset.PTUR(t,n) + final_output_std.cluster_offset.PTUR(t,n);
            final_output_ub.cluster_offset.PCTO(t,n) = final_output_thigh_mean.cluster_offset.PCTO(t,n) + final_output_std.cluster_offset.PCTO(t,n);
            
            final_output_lb.cluster_offset.SVDLS(t,n) = final_output_thigh_mean.cluster_offset.SVDLS(t,n) - final_output_std.cluster_offset.SVDLS(t,n);
            final_output_lb.cluster_offset.PTUR(t,n) = final_output_thigh_mean.cluster_offset.PTUR(t,n) - final_output_std.cluster_offset.PTUR(t,n);
            final_output_lb.cluster_offset.PCTO(t,n) = final_output_thigh_mean.cluster_offset.PCTO(t,n) - final_output_std.cluster_offset.PCTO(t,n);
            
        end

        for n=1:input_size.Num_AL_Thigh
            
            final_output_thigh_mean.AL_offset.SVDLS(t,n) = nanmean(AL_offset_current.SVDLS(:,n));
            final_output_thigh_mean.AL_offset.PTUR(t,n) = nanmean(AL_offset_current.PTUR(:,n));
            final_output_thigh_mean.AL_offset.PCTO(t,n) = nanmean(AL_offset_current.PCTO(:,n));
            
            final_output_std.AL_offset.SVDLS(t,n) = std(AL_offset_current.SVDLS(:,n), 'omitnan');
            final_output_std.AL_offset.PTUR(t,n) = std(AL_offset_current.PTUR(:,n), 'omitnan');
            final_output_std.AL_offset.PCTO(t,n) = std(AL_offset_current.PCTO(:,n), 'omitnan');
            
            final_output_ub.AL_offset.SVDLS(t,n) = final_output_thigh_mean.cluster_offset.SVDLS(t,n) + final_output_std.AL_offset.SVDLS(t,n);
            final_output_ub.AL_offset.PTUR(t,n) = final_output_thigh_mean.AL_offset.PTUR(t,n) + final_output_std.AL_offset.PTUR(t,n);
            final_output_ub.AL_offset.PCTO(t,n) = final_output_thigh_mean.AL_offset.PCTO(t,n) + final_output_std.AL_offset.PCTO(t,n);
            
            final_output_lb.AL_offset.SVDLS(t,n) = final_output_thigh_mean.AL_offset.SVDLS(t,n) - final_output_std.AL_offset.SVDLS(t,n);
            final_output_lb.AL_offset.PTUR(t,n) = final_output_thigh_mean.AL_offset.PTUR(t,n) - final_output_std.AL_offset.PTUR(t,n);
            final_output_lb.AL_offset.PCTO(t,n) = final_output_thigh_mean.AL_offset.PCTO(t,n) - final_output_std.AL_offset.PCTO(t,n);
            
        end

        for n=1:NumEigenvectors
            
            final_output_thigh_mean.eigenvector_offset.SVDLS(t,n) = nanmean(eigenvector_offset_current.SVDLS(:,n));
            final_output_thigh_mean.eigenvector_offset.PTUR(t,n) = nanmean(eigenvector_offset_current.PTUR(:,n));
            final_output_thigh_mean.eigenvector_offset.PCTO(t,n) = nanmean(eigenvector_offset_current.PCTO(:,n));
            
            final_output_std.eigenvector_offset.SVDLS(t,n) = std(eigenvector_offset_current.SVDLS(:,n), 'omitnan');
            final_output_std.eigenvector_offset.PTUR(t,n) = std(eigenvector_offset_current.PTUR(:,n), 'omitnan');
            final_output_std.eigenvector_offset.PCTO(t,n) = std(eigenvector_offset_current.PCTO(:,n), 'omitnan');
            
            final_output_ub.eigenvector_offset.SVDLS(t,n) = final_output_thigh_mean.eigenvector_offset.SVDLS(t,n) + final_output_std.eigenvector_offset.SVDLS(t,n);
            final_output_ub.eigenvector_offset.PTUR(t,n) = final_output_thigh_mean.eigenvector_offset.PTUR(t,n) + final_output_std.eigenvector_offset.PTUR(t,n);
            final_output_ub.eigenvector_offset.PCTO(t,n) = final_output_thigh_mean.eigenvector_offset.PCTO(t,n) + final_output_std.eigenvector_offset.PCTO(t,n);
            
            final_output_lb.eigenvector_offset.SVDLS(t,n) = final_output_thigh_mean.eigenvector_offset.SVDLS(t,n) - final_output_std.eigenvector_offset.SVDLS(t,n);
            final_output_lb.eigenvector_offset.PTUR(t,n) = final_output_thigh_mean.eigenvector_offset.PTUR(t,n) - final_output_std.eigenvector_offset.PTUR(t,n);
            final_output_lb.eigenvector_offset.PCTO(t,n) = final_output_thigh_mean.eigenvector_offset.PCTO(t,n) - final_output_std.eigenvector_offset.PCTO(t,n);
            
        end
    end
    
    PCT_configurations = 0;
    for k=1:NumConfigurations
        if ~final_output(k).discontinuity_flag_PCT
            PCT_configurations = PCT_configurations + 1;
        end
    end
    
    inBetween = struct;
    
    inBetween.CM_offset = struct;
    inBetween.CM_offset.SVDLS = vertcat(final_output_ub.CM_offset.SVDLS, flipud(final_output_lb.CM_offset.SVDLS));
    inBetween.CM_offset.PTUR = vertcat(final_output_ub.CM_offset.PTUR, flipud(final_output_lb.CM_offset.PTUR));
    inBetween.CM_offset.PCTO = vertcat(final_output_ub.CM_offset.PCTO, flipud(final_output_lb.CM_offset.PCTO));
    
    inBetween.cluster_offset = struct;
    inBetween.AL_offset = struct;
    inBetween.eigenvector_offset = struct;
    
    N = input_size.N;
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N) + "\Thigh\";
    path_avg = path_img + "\Average_plots\";
    
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\PCT paper drafts\Images\";
    path_avg_final = path_img_final;

    figure('units','normalized','outerposition',[0 0 1 1]);
    xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    plot(fraction_completed_gait_cycle, final_output_thigh_mean.CM_offset.SVDLS, 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
    plot(fraction_completed_gait_cycle, final_output_thigh_mean.CM_offset.PTUR, 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
    plot(fraction_completed_gait_cycle, final_output_thigh_mean.CM_offset.PCTO, 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
    fill(x2, inBetween.CM_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
    fill(x2, inBetween.CM_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
    fill(x2(2:end-1), inBetween.CM_offset.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
    fill([0 0.6 0.6 0], [0 0 3.5 3.5], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim([0, 3.5]);
    grid;
    if display_flag.legend_bool
        legend('Location', 'best');
    end
    if display_flag.title_bool
        title('Average CM offset during gait cycle for Thigh' + " , PCT configurations = " + string(PCT_configurations));
    end
    ax = gca;
    ax.FontSize = 20;
    xlabel('Fraction completed gait cycle');
    ylabel('Error (cm)');
    filename_CMAverage = path_avg + "CMAverage_plot.png";
    exportgraphics(gcf, filename_CMAverage);
    filename_CMAverage_final = path_avg_final + "CM_offset_thigh.png";
    exportgraphics(gcf, filename_CMAverage_final);
    close(gcf);

    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for n=1:input_size.N
        inBetween.cluster_offset.SVDLS = vertcat(final_output_ub.cluster_offset.SVDLS(:,n), flipud(final_output_lb.cluster_offset.SVDLS(:,n)));
        inBetween.cluster_offset.PTUR = vertcat(final_output_ub.cluster_offset.PTUR(:,n), flipud(final_output_lb.cluster_offset.PTUR(:,n)));
        inBetween.cluster_offset.PCTO = vertcat(final_output_ub.cluster_offset.PCTO(:,n), flipud(final_output_lb.cluster_offset.PCTO(:,n)));
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.cluster_offset.SVDLS(:,n), 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.cluster_offset.PTUR(:,n), 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.cluster_offset.PCTO(:,n), 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
        fill(x2, inBetween.cluster_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
        fill(x2, inBetween.cluster_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
        fill(x2(2:end-1), inBetween.cluster_offset.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
        fill([0 0.6 0.6 0], [0 0 4 4], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim([0, 4]);
        grid;
        if n==2
            if display_flag.legend_bool
                legend('Location', 'best');
            end
        end
        if n>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(n,2) == 0)
            ylabel("Error (cm)");
        end
        title("Point IDX = " + string(n));
        ax = gca;
        ax.FontSize = 20;
    end
    if display_flag.title_bool
        suptitle("Average Cluster Point Error Plots for Thigh, PCT configurations = " + string(PCT_configurations), 'FontSize', 20);
    end
    filename_CPAverage = path_avg + "CPAverage_plot.png";
    exportgraphics(gcf, filename_CPAverage);
    filename_CPAverage_final = path_avg_final + "ClusterPt_offset_thigh.png";
    exportgraphics(gcf, filename_CPAverage_final);
    close(gcf);

    AL_name_list = ["Lateral Epicondyle";
                     "Medial Epicondyle";
                     "Greater Trochanter";
                     "Femoral Head"];

    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for n=1:input_size.Num_AL_Thigh
        final_output_ub.AL_offset.SVDLS = final_output_thigh_mean.AL_offset.SVDLS(2:end,n) + final_output_std.AL_offset.SVDLS(2:end,n);
        final_output_ub.AL_offset.PTUR = final_output_thigh_mean.AL_offset.PTUR(2:end,n) + final_output_std.AL_offset.PTUR(2:end,n);
        final_output_ub.AL_offset.PCTO = final_output_thigh_mean.AL_offset.PCTO(2:end,n) + final_output_std.AL_offset.PCTO(2:end,n);
        
        final_output_lb.AL_offset.SVDLS = final_output_thigh_mean.AL_offset.SVDLS(2:end,n) - final_output_std.AL_offset.SVDLS(2:end,n);
        final_output_lb.AL_offset.PTUR = final_output_thigh_mean.AL_offset.PTUR(2:end,n) - final_output_std.AL_offset.PTUR(2:end,n);
        final_output_lb.AL_offset.PCTO = final_output_thigh_mean.AL_offset.PCTO(2:end,n) - final_output_std.AL_offset.PCTO(2:end,n);
        
        inBetween.AL_offset.SVDLS = vertcat(final_output_ub.AL_offset.SVDLS, flipud(final_output_lb.AL_offset.SVDLS));
        inBetween.AL_offset.PTUR = vertcat(final_output_ub.AL_offset.PTUR, flipud(final_output_lb.AL_offset.PTUR));
        inBetween.AL_offset.PCTO = vertcat(final_output_ub.AL_offset.PCTO, flipud(final_output_lb.AL_offset.PCTO));
        
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.AL_offset.SVDLS(:,n), 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.AL_offset.PTUR(:,n), 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.AL_offset.PCTO(:,n), 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
        fill(x2_short, inBetween.AL_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
        fill(x2_short, inBetween.AL_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
        fill(x2_short, inBetween.AL_offset.PCTO, color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
        fill([0 0.6 0.6 0], [0 0 4 4], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim([0, 4]);
        grid;
        if n==2
            if display_flag.legend_bool
                legend('Location', 'best');
            end
        end
        if n>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(n,2) == 0)
            ylabel("Error (cm)");
        end
        title(AL_name_list(n));
        ax = gca;
        ax.FontSize = 20;
    end
    if display_flag.title_bool
        suptitle("Average Reconstruction Error Plots for Anatomical Landmarks for Thigh, PCT configurations = " + string(PCT_configurations), 'FontSize', 20);
    end
    filename_ALAverage = path_avg + "ALAverage_plot.png";
    exportgraphics(gcf, filename_ALAverage);
    filename_ALAverage_final = path_avg_final + "AL_offset_thigh.png";
    exportgraphics(gcf, filename_ALAverage_final);
    close(gcf);

    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(3,1, 'TileSpacing', 'compact');
    for n=1:NumEigenvectors
        final_output_ub.eigenvector_offset.SVDLS = final_output_thigh_mean.eigenvector_offset.SVDLS(2:end,n) + final_output_std.eigenvector_offset.SVDLS(2:end,n);
        final_output_ub.eigenvector_offset.PTUR = final_output_thigh_mean.eigenvector_offset.PTUR(2:end,n) + final_output_std.eigenvector_offset.PTUR(2:end,n);
        final_output_ub.eigenvector_offset.PCTO = final_output_thigh_mean.eigenvector_offset.PCTO(2:end,n) + final_output_std.eigenvector_offset.PCTO(2:end,n);
        
        final_output_lb.eigenvector_offset.SVDLS = final_output_thigh_mean.eigenvector_offset.SVDLS(2:end,n) - final_output_std.eigenvector_offset.SVDLS(2:end,n);
        final_output_lb.eigenvector_offset.PTUR = final_output_thigh_mean.eigenvector_offset.PTUR(2:end,n) - final_output_std.eigenvector_offset.PTUR(2:end,n);
        final_output_lb.eigenvector_offset.PCTO = final_output_thigh_mean.eigenvector_offset.PCTO(2:end,n) - final_output_std.eigenvector_offset.PCTO(2:end,n);
        
        inBetween.eigenvector_offset.SVDLS = vertcat(final_output_ub.eigenvector_offset.SVDLS, flipud(final_output_lb.eigenvector_offset.SVDLS));
        inBetween.eigenvector_offset.PTUR = vertcat(final_output_ub.eigenvector_offset.PTUR, flipud(final_output_lb.eigenvector_offset.PTUR));
        inBetween.eigenvector_offset.PCTO = vertcat(final_output_ub.eigenvector_offset.PCTO, flipud(final_output_lb.eigenvector_offset.PCTO));
        
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle(2:end), final_output_thigh_mean.eigenvector_offset.SVDLS(2:end,n), 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.eigenvector_offset.PTUR(:,n), 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
        plot(fraction_completed_gait_cycle, final_output_thigh_mean.eigenvector_offset.PCTO(:,n), 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
        fill(x2_short, inBetween.eigenvector_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
        fill(x2_short, inBetween.eigenvector_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
        fill(x2_short, inBetween.eigenvector_offset.PCTO, color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
        fill([0 0.6 0.6 0], [0 0 10 10], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim([0, 10]);
        grid;
        if n==1
            if display_flag.legend_bool
                legend('Location', 'best');
            end
        end
        if n==3
            xlabel("Fraction completed gait cycle");
        end
        ylabel('Error (degress)');
        title("Eigenvector IDX = " + string(n));
        ax = gca;
        ax.FontSize = 20;
    end
    if display_flag.title_bool
        suptitle("Average Reconstruction Error Plots for Eigenvectors for Thigh, PCT configurations = " + string(PCT_configurations), 'FontSize', 20);
    end
    filename_EVAverage = path_avg + "EVAverage_plot.png";
    exportgraphics(gcf, filename_EVAverage);
    filename_EVAverage_final = path_avg_final + "Eigenvector_offset_thigh.png";
    exportgraphics(gcf, filename_EVAverage_final);
    close(gcf);

end

function final_output_shank_mean = generate_ouput_stats_shank(final_output_shank, input_size, display_flag)

    final_output_shank_mean = struct;
    
    color_info = display_flag.color_info;
    
    final_output_shank_mean.T_estimated = struct;
    final_output_shank_mean.T_estimated.SVDLS = nan(input_size.NumFrames, 3, 1);
    final_output_shank_mean.T_estimated.PTUR = nan(input_size.NumFrames, 3, 1);
    final_output_shank_mean.T_estimated.PCTO = nan(input_size.NumFrames, 3, 1);

    final_output_shank_mean.CM_offset = struct;
    final_output_shank_mean.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_shank_mean.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_shank_mean.CM_offset.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_std.CM_offset = struct;
    final_output_std.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_std.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_std.CM_offset.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_ub.CM_offset = struct;
    final_output_ub.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_ub.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_ub.CM_offset.PCTO = nan(input_size.NumFrames, 1);
    
    final_output_lb.CM_offset = struct;
    final_output_lb.CM_offset.SVDLS = nan(input_size.NumFrames, 1);
    final_output_lb.CM_offset.PTUR = nan(input_size.NumFrames, 1);
    final_output_lb.CM_offset.PCTO = nan(input_size.NumFrames, 1);

    final_output_shank_mean.cluster_offset = struct;
    final_output_shank_mean.cluster_offset.SVDLS = nan(input_size.NumFrames, input_size.N);
    final_output_shank_mean.cluster_offset.PTUR = nan(input_size.NumFrames, input_size.N);
    final_output_shank_mean.cluster_offset.PCTO = nan(input_size.NumFrames, input_size.N);

    final_output_shank_mean.AL_offset = struct;
    final_output_shank_mean.AL_offset.SVDLS = nan(input_size.NumFrames, input_size.Num_AL_Shank);
    final_output_shank_mean.AL_offset.PTUR = nan(input_size.NumFrames, input_size.Num_AL_Shank);
    final_output_shank_mean.AL_offset.PCTO = nan(input_size.NumFrames, input_size.Num_AL_Shank);

    NumEigenvectors = 3;
    final_output_shank_mean.eigenvector_offset = struct;
    final_output_shank_mean.eigenvector_offset.SVDLS = nan(input_size.NumFrames, NumEigenvectors);
    final_output_shank_mean.eigenvector_offset.PTUR = nan(input_size.NumFrames, NumEigenvectors);
    final_output_shank_mean.eigenvector_offset.PCTO = nan(input_size.NumFrames, NumEigenvectors);
    
    NumConfigurations = min([length(final_output_shank)-1, input_size.NumConfigurations]);
    
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    x2 = vertcat(fraction_completed_gait_cycle', flipud(fraction_completed_gait_cycle'));
    x2_short = vertcat(fraction_completed_gait_cycle(2:end)', flipud(fraction_completed_gait_cycle(2:end)'));
    ToeOff = 0.6;

    for t=1:input_size.NumFrames
        disp(t);
        
        T_estimated_current = struct;
        T_estimated_current.SVDLS = nan(NumConfigurations, 3, 1);
        T_estimated_current.PTUR = nan(NumConfigurations, 3, 1);
        T_estimated_current.PCTO = nan(NumConfigurations, 3, 1);
        
        CM_offset_current = struct;
        CM_offset_current.SVDLS = nan(NumConfigurations, 1);
        CM_offset_current.PTUR = nan(NumConfigurations, 1);
        CM_offset_current.PCTO = nan(NumConfigurations, 1);

        cluster_offset_current = struct;
        cluster_offset_current.SVDLS = nan(NumConfigurations, input_size.N);
        cluster_offset_current.PTUR = nan(NumConfigurations, input_size.N);
        cluster_offset_current.PCTO = nan(NumConfigurations, input_size.N);

        AL_offset_current = struct;
        AL_offset_current.SVDLS = nan(NumConfigurations, input_size.Num_AL_Shank);
        AL_offset_current.PTUR = nan(NumConfigurations, input_size.Num_AL_Shank);
        AL_offset_current.PCTO = nan(NumConfigurations, input_size.Num_AL_Shank);

        eigenvector_offset_current = struct;
        eigenvector_offset_current.SVDLS = nan(NumConfigurations, NumEigenvectors);
        eigenvector_offset_current.PTUR = nan(NumConfigurations, NumEigenvectors);
        eigenvector_offset_current.PCTO = nan(NumConfigurations, NumEigenvectors);

        for k=1:NumConfigurations
            disp(k);
            T_estimated_current.SVDLS(k,:,1) = final_output_shank(k).output_containers.T_estimated.SVDLS(t);
            T_estimated_current.PTUR(k,:,1) = final_output_shank(k).output_containers.T_estimated.PTUR(t);
            CM_offset_current.SVDLS(k) = final_output_shank(k).output_containers.CM_offset.SVDLS(t);
            CM_offset_current.PTUR(k) = final_output_shank(k).output_containers.CM_offset.PTUR(t);
            
            if ~final_output_shank(k).discontinuity_flag_PCT
                
                CM_offset_current.PCTO(k) = final_output_shank(k).output_containers.CM_offset.PCTO(t);
                
                for n=1:input_size.N
                    cluster_offset_current.PCTO(k,n) = final_output_shank(k).output_containers.cluster_offset.PCTO(t,n);
                end
                
                for n=1:input_size.Num_AL_Shank
                    AL_offset_current.PCTO(k,n) = final_output_shank(k).output_containers.AL_offset.PCTO(t,n);
                end
                
                for n=1:NumEigenvectors
                    eigenvector_offset_current.PCTO(k,n) = final_output_shank(k).output_containers.eigenvector_offset.PCTO(t,n);
                end
                
            end
            
            for n=1:input_size.N
                cluster_offset_current.SVDLS(k,n) = final_output_shank(k).output_containers.cluster_offset.SVDLS(t,n);
                cluster_offset_current.PTUR(k,n) = final_output_shank(k).output_containers.cluster_offset.PTUR(t,n);
            end

            for n=1:input_size.Num_AL_Shank
                AL_offset_current.SVDLS(k,n) = final_output_shank(k).output_containers.AL_offset.SVDLS(t,n);
                AL_offset_current.PTUR(k,n) = final_output_shank(k).output_containers.AL_offset.PTUR(t,n);
            end

            for n=1:NumEigenvectors
                eigenvector_offset_current.SVDLS(k,n) = final_output_shank(k).output_containers.eigenvector_offset.SVDLS(t,n);
                eigenvector_offset_current.PTUR(k,n) = final_output_shank(k).output_containers.eigenvector_offset.PTUR(t,n);
            end

        end
        final_output_shank_mean.CM_offset.SVDLS(t) = nanmean(CM_offset_current.SVDLS);
        final_output_shank_mean.CM_offset.PTUR(t) = nanmean(CM_offset_current.PTUR);
        final_output_shank_mean.CM_offset.PCTO(t) = nanmean(CM_offset_current.PCTO);
        
        final_output_std.CM_offset.SVDLS(t) = std(CM_offset_current.SVDLS, 'omitnan');
        final_output_std.CM_offset.PTUR(t) = std(CM_offset_current.PTUR, 'omitnan');
        final_output_std.CM_offset.PCTO(t) = std(CM_offset_current.PCTO, 'omitnan');
        
        final_output_ub.CM_offset.SVDLS(t) = final_output_shank_mean.CM_offset.SVDLS(t) + final_output_std.CM_offset.SVDLS(t);
        final_output_ub.CM_offset.PTUR(t) = final_output_shank_mean.CM_offset.PTUR(t) + final_output_std.CM_offset.PTUR(t);
        final_output_ub.CM_offset.PCTO(t) = final_output_shank_mean.CM_offset.PCTO(t) + final_output_std.CM_offset.PCTO(t);
        
        final_output_lb.CM_offset.SVDLS(t) = final_output_shank_mean.CM_offset.SVDLS(t) - final_output_std.CM_offset.SVDLS(t);
        final_output_lb.CM_offset.PTUR(t) = final_output_shank_mean.CM_offset.PTUR(t) - final_output_std.CM_offset.PTUR(t);
        final_output_lb.CM_offset.PCTO(t) = final_output_shank_mean.CM_offset.PCTO(t) - final_output_std.CM_offset.PCTO(t);

        for n=1:input_size.N
            
            final_output_shank_mean.cluster_offset.SVDLS(t,n) = nanmean(cluster_offset_current.SVDLS(:,n));
            final_output_shank_mean.cluster_offset.PTUR(t,n) = nanmean(cluster_offset_current.PTUR(:,n));
            final_output_shank_mean.cluster_offset.PCTO(t,n) = nanmean(cluster_offset_current.PCTO(:,n));
            
            final_output_std.cluster_offset.SVDLS(t,n) = std(cluster_offset_current.SVDLS(:,n), 'omitnan');
            final_output_std.cluster_offset.PTUR(t,n) = std(cluster_offset_current.PTUR(:,n), 'omitnan');
            final_output_std.cluster_offset.PCTO(t,n) = std(cluster_offset_current.PCTO(:,n), 'omitnan');
            
            final_output_ub.cluster_offset.SVDLS(t,n) = final_output_shank_mean.cluster_offset.SVDLS(t,n) + final_output_std.cluster_offset.SVDLS(t,n);
            final_output_ub.cluster_offset.PTUR(t,n) = final_output_shank_mean.cluster_offset.PTUR(t,n) + final_output_std.cluster_offset.PTUR(t,n);
            final_output_ub.cluster_offset.PCTO(t,n) = final_output_shank_mean.cluster_offset.PCTO(t,n) + final_output_std.cluster_offset.PCTO(t,n);
            
            final_output_lb.cluster_offset.SVDLS(t,n) = final_output_shank_mean.cluster_offset.SVDLS(t,n) - final_output_std.cluster_offset.SVDLS(t,n);
            final_output_lb.cluster_offset.PTUR(t,n) = final_output_shank_mean.cluster_offset.PTUR(t,n) - final_output_std.cluster_offset.PTUR(t,n);
            final_output_lb.cluster_offset.PCTO(t,n) = final_output_shank_mean.cluster_offset.PCTO(t,n) - final_output_std.cluster_offset.PCTO(t,n);
            
        end

        for n=1:input_size.Num_AL_Shank
            
            final_output_shank_mean.AL_offset.SVDLS(t,n) = nanmean(AL_offset_current.SVDLS(:,n));
            final_output_shank_mean.AL_offset.PTUR(t,n) = nanmean(AL_offset_current.PTUR(:,n));
            final_output_shank_mean.AL_offset.PCTO(t,n) = nanmean(AL_offset_current.PCTO(:,n));
            
            final_output_std.AL_offset.SVDLS(t,n) = std(AL_offset_current.SVDLS(:,n), 'omitnan');
            final_output_std.AL_offset.PTUR(t,n) = std(AL_offset_current.PTUR(:,n), 'omitnan');
            final_output_std.AL_offset.PCTO(t,n) = std(AL_offset_current.PCTO(:,n), 'omitnan');
            
            final_output_ub.AL_offset.SVDLS(t,n) = final_output_shank_mean.AL_offset.SVDLS(t,n) + final_output_std.AL_offset.SVDLS(t,n);
            final_output_ub.AL_offset.PTUR(t,n) = final_output_shank_mean.AL_offset.PTUR(t,n) + final_output_std.AL_offset.PTUR(t,n);
            final_output_ub.AL_offset.PCTO(t,n) = final_output_shank_mean.AL_offset.PCTO(t,n) + final_output_std.AL_offset.PCTO(t,n);
            
            final_output_lb.AL_offset.SVDLS(t,n) = final_output_shank_mean.AL_offset.SVDLS(t,n) - final_output_std.AL_offset.SVDLS(t,n);
            final_output_lb.AL_offset.PTUR(t,n) = final_output_shank_mean.AL_offset.PTUR(t,n) - final_output_std.AL_offset.PTUR(t,n);
            final_output_lb.AL_offset.PCTO(t,n) = final_output_shank_mean.AL_offset.PCTO(t,n) - final_output_std.AL_offset.PCTO(t,n);
            
        end

        for n=1:NumEigenvectors
            
            final_output_shank_mean.eigenvector_offset.SVDLS(t,n) = nanmean(eigenvector_offset_current.SVDLS(:,n));
            final_output_shank_mean.eigenvector_offset.PTUR(t,n) = nanmean(eigenvector_offset_current.PTUR(:,n));
            final_output_shank_mean.eigenvector_offset.PCTO(t,n) = nanmean(eigenvector_offset_current.PCTO(:,n));
            
            final_output_std.eigenvector_offset.SVDLS(t,n) = std(eigenvector_offset_current.SVDLS(:,n), 'omitnan');
            final_output_std.eigenvector_offset.PTUR(t,n) = std(eigenvector_offset_current.PTUR(:,n), 'omitnan');
            final_output_std.eigenvector_offset.PCTO(t,n) = std(eigenvector_offset_current.PCTO(:,n), 'omitnan');
            
            final_output_ub.eigenvector_offset.SVDLS(t,n) = final_output_shank_mean.eigenvector_offset.SVDLS(t,n) + final_output_std.eigenvector_offset.SVDLS(t,n);
            final_output_ub.eigenvector_offset.PTUR(t,n) = final_output_shank_mean.eigenvector_offset.PTUR(t,n) + final_output_std.eigenvector_offset.PTUR(t,n);
            final_output_ub.eigenvector_offset.PCTO(t,n) = final_output_shank_mean.eigenvector_offset.PCTO(t,n) + final_output_std.eigenvector_offset.PCTO(t,n);
            
            final_output_lb.eigenvector_offset.SVDLS(t,n) = final_output_shank_mean.eigenvector_offset.SVDLS(t,n) - final_output_std.eigenvector_offset.SVDLS(t,n);
            final_output_lb.eigenvector_offset.PTUR(t,n) = final_output_shank_mean.eigenvector_offset.PTUR(t,n) - final_output_std.eigenvector_offset.PTUR(t,n);
            final_output_lb.eigenvector_offset.PCTO(t,n) = final_output_shank_mean.eigenvector_offset.PCTO(t,n) - final_output_std.eigenvector_offset.PCTO(t,n);
            
        end
    end
    
    PCT_configurations = 0;
    for k=1:NumConfigurations
        if ~final_output_shank(k).discontinuity_flag_PCT
            PCT_configurations = PCT_configurations + 1;
        end
    end
    
    inBetween = struct;
    
    inBetween.CM_offset = struct;
    inBetween.CM_offset.SVDLS = vertcat(final_output_ub.CM_offset.SVDLS, flipud(final_output_lb.CM_offset.SVDLS));
    inBetween.CM_offset.PTUR = vertcat(final_output_ub.CM_offset.PTUR, flipud(final_output_lb.CM_offset.PTUR));
    inBetween.CM_offset.PCTO = vertcat(final_output_ub.CM_offset.PCTO, flipud(final_output_lb.CM_offset.PCTO));
    
    inBetween.cluster_offset = struct;
    inBetween.AL_offset = struct;
    inBetween.eigenvector_offset = struct;
    
    N = input_size.N;
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N) + "\Shank\";
    path_avg = path_img + "\Average_plots\";
    
    path_img_final = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\PCT paper drafts\Images\";
    path_avg_final = path_img_final;

    figure('units','normalized','outerposition',[0 0 1 1]);
    xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    plot(fraction_completed_gait_cycle, final_output_shank_mean.CM_offset.SVDLS, 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
    plot(fraction_completed_gait_cycle, final_output_shank_mean.CM_offset.PTUR, 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
    plot(fraction_completed_gait_cycle, final_output_shank_mean.CM_offset.PCTO, 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
    fill(x2, inBetween.CM_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
    fill(x2, inBetween.CM_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
    fill(x2(2:end-1), inBetween.CM_offset.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
    fill([0 0.6 0.6 0], [0 0 3.5 3.5], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim([0, 3.5]);
    grid;
    if display_flag.legend_bool
        legend('Location', 'best');
    end
    if display_flag.title_bool
        title('Average CM offset during gait cycle for Shank' + " , PCT configurations = " + string(PCT_configurations));
    end
    xlabel('Fraction completed gait cycle');
    ylabel('Error (cm)');
    ax = gca;
    ax.FontSize = 20;
    filename_CMAverage = path_avg + "CMAverage_plot.png";
    exportgraphics(gcf, filename_CMAverage);
    filename_CMAverage_final = path_avg_final + "CM_offset_shank.png";
    exportgraphics(gcf, filename_CMAverage_final);
    close(gcf);

    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for n=1:input_size.N
        inBetween.cluster_offset.SVDLS = vertcat(final_output_ub.cluster_offset.SVDLS(:,n), flipud(final_output_lb.cluster_offset.SVDLS(:,n)));
        inBetween.cluster_offset.PTUR = vertcat(final_output_ub.cluster_offset.PTUR(:,n), flipud(final_output_lb.cluster_offset.PTUR(:,n)));
        inBetween.cluster_offset.PCTO = vertcat(final_output_ub.cluster_offset.PCTO(:,n), flipud(final_output_lb.cluster_offset.PCTO(:,n)));
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_shank_mean.cluster_offset.SVDLS(:,n), 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
        plot(fraction_completed_gait_cycle, final_output_shank_mean.cluster_offset.PTUR(:,n), 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
        plot(fraction_completed_gait_cycle, final_output_shank_mean.cluster_offset.PCTO(:,n), 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
        fill(x2, inBetween.cluster_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
        fill(x2, inBetween.cluster_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
        fill(x2(2:end-1), inBetween.cluster_offset.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
        fill([0 0.6 0.6 0], [0 0 4 4], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim([0, 4]);
        grid;
        if n==2
            if display_flag.legend_bool
                legend('Location', 'best');
            end
        end
        if n>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(n,2) == 0)
            ylabel("Error (cm)");
        end
        title("Point IDX = " + string(n));
        ax = gca;
        ax.FontSize = 20;
    end
    if display_flag.title_bool
        suptitle("Average Cluster Point Error Plots for Shank, PCT configurations = " + string(PCT_configurations), 'FontSize', 20);
    end
    filename_CPAverage = path_avg + "CPAverage_plot.png";
    exportgraphics(gcf, filename_CPAverage);
    filename_CPAverage_final = path_avg_final + "ClusterPt_offset_shank.png";
    exportgraphics(gcf, filename_CPAverage_final);
    close(gcf);

    AL_name_list = ["Lateral Malleolus";
                 "Medial Malleolus";
                 "Fibular Head";
                 "Tibial Tuberosity"];

    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for n=1:input_size.Num_AL_Shank
        inBetween.AL_offset.SVDLS = vertcat(final_output_ub.AL_offset.SVDLS(:,n), flipud(final_output_lb.AL_offset.SVDLS(:,n)));
        inBetween.AL_offset.PTUR = vertcat(final_output_ub.AL_offset.PTUR(:,n), flipud(final_output_lb.AL_offset.PTUR(:,n)));
        inBetween.AL_offset.PCTO = vertcat(final_output_ub.AL_offset.PCTO(:,n), flipud(final_output_lb.AL_offset.PCTO(:,n)));
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_shank_mean.AL_offset.SVDLS(:,n), 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
        plot(fraction_completed_gait_cycle, final_output_shank_mean.AL_offset.PTUR(:,n), 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
        plot(fraction_completed_gait_cycle, final_output_shank_mean.AL_offset.PCTO(:,n), 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
        fill(x2, inBetween.AL_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
        fill(x2, inBetween.AL_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
        fill(x2(2:end-1), inBetween.AL_offset.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
        fill([0 0.6 0.6 0], [0 0 4 4], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim([0, 4]);
        grid;
        if n==2
            if display_flag.legend_bool
                legend('Location', 'best');
            end
        end
        if n>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(n,2) == 0)
            ylabel("Error (cm)");
        end
        title(AL_name_list(n));
        ax = gca;
        ax.FontSize = 20;
    end
    if display_flag.title_bool
        suptitle("Average Reconstruction Error Plots for Anatomical Landmarks for Shank, PCT configurations = " + string(PCT_configurations), 'FontSize', 20);
    end
    filename_ALAverage = path_avg + "ALAverage_plot.png";
    exportgraphics(gcf, filename_ALAverage);
    filename_ALAverage_final = path_avg_final + "AL_offset_shank.png";
    exportgraphics(gcf, filename_ALAverage_final);
    close(gcf);

    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(3,1, 'TileSpacing', 'compact');
    for n=1:NumEigenvectors
        inBetween.eigenvector_offset.SVDLS = vertcat(final_output_ub.eigenvector_offset.SVDLS(:,n), flipud(final_output_lb.eigenvector_offset.SVDLS(:,n)));
        inBetween.eigenvector_offset.PTUR = vertcat(final_output_ub.eigenvector_offset.PTUR(:,n), flipud(final_output_lb.eigenvector_offset.PTUR(:,n)));
        inBetween.eigenvector_offset.PCTO = vertcat(final_output_ub.eigenvector_offset.PCTO(:,n), flipud(final_output_lb.eigenvector_offset.PCTO(:,n)));
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        plot(fraction_completed_gait_cycle, final_output_shank_mean.eigenvector_offset.SVDLS(:,n), 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
        plot(fraction_completed_gait_cycle, final_output_shank_mean.eigenvector_offset.PTUR(:,n), 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
        plot(fraction_completed_gait_cycle, final_output_shank_mean.eigenvector_offset.PCTO(:,n), 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
        fill(x2, inBetween.eigenvector_offset.SVDLS, color_info.SVDLS, 'FaceAlpha', 0.2, 'EdgeColor', color_info.SVDLS, 'HandleVisibility', 'off');
        fill(x2, inBetween.eigenvector_offset.PTUR, color_info.PTUR, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PTUR, 'HandleVisibility', 'off');
        fill(x2(2:end-1), inBetween.eigenvector_offset.PCTO(2:end-1), color_info.PCTO, 'FaceAlpha', 0.2, 'EdgeColor', color_info.PCTO, 'HandleVisibility', 'off');
        fill([0 0.6 0.6 0], [0 0 10 10], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim([0, 10]);
        grid;
        if n==1
            if display_flag.legend_bool
                legend('Location', 'best');
            end
        end
        if n==3
            xlabel("Fraction completed gait cycle");
        end
        ylabel('Error (degress)');
        title("Eigenvector IDX = " + string(n));
        ax = gca;
        ax.FontSize = 20;
    end
    if display_flag.title_bool
        suptitle("Average Reconstruction Error Plots for Eigenvectors for Shank, PCT configurations = " + string(PCT_configurations), 'FontSize', 20);
    end
    filename_EVAverage = path_avg + "EVAverage_plot.png";
    exportgraphics(gcf, filename_EVAverage);
    filename_EVAverage_final = path_avg_final + "Eigenvector_offset_shank.png";
    exportgraphics(gcf, filename_EVAverage_final);
    close(gcf);

end

function generate_output_plots_thigh(input_struct, output_containers_thigh)
    input_size = input_struct.input_size;
    input_data = input_struct.input_data_thigh;

    configuration_idx = input_data.configuration_idx;
    N = input_size.N;
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    ToeOff = 0.6;
    
    color_info = input_struct.color_info;
    alg_list = ["SVDLS", "PTUR", "PCTO"];
    alg_legend_list = ["SVD-LS", "PCT-PT", "PCT"];

    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N) + "\Thigh\";
    
    path_CM = path_img + "CM_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    for n=1:length(alg_list)
        current_alg = alg_list(n);
        current_alg_legend_name = alg_legend_list(n);
        current_data = output_containers_thigh.CM_offset.(current_alg);
        if current_alg == "PTNR"
            scatter(fraction_completed_gait_cycle, current_data, 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
        else
            plot(fraction_completed_gait_cycle, current_data, 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
        end
    end
    yL =  ylim;
    fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim(yL);
    grid;
    legend('Location', 'best');
    xlabel("Fraction completed gait cycle");
    ylabel("Error (cm)");
    %title("Center of Mass offset for Thigh, IDX = " + string(configuration_idx));
    ax = gca;
    ax.FontSize = 20;
    filename_CM = path_CM + "CM_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_CM);
    close(gcf);

    path_Eigenvector = path_img + "Eigenvector_plots\";
    figure('units','normalized','outerposition',[0 0 1 1])
    tiledlayout(3,1, 'TileSpacing', 'compact');
    for j=1:3
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_alg_legend_name = alg_legend_list(n);
            current_data = output_containers_thigh.eigenvector_offset.(current_alg);
            plot(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
        end
        yL =  ylim;
        fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim(yL);
        grid;
        if j==1
            legend('Location', 'best');
        end
        if j==3
            xlabel("Fraction completed gait cycle");
        end
        ylabel("Error (degrees)");
        title("Eigenvector " + string(j));
        ax = gca;
        ax.FontSize = 20;
    end
    %suptitle("Eigenvector Offset for Thigh, IDX = " + string(configuration_idx), 'FontSize', 20);
    filename_Eigenvector = path_Eigenvector + "Eigenvector_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_Eigenvector);
    close(gcf);
    
    AL_name_list = ["Lateral Epicondyle";
                 "Medial Epicondyle";
                 "Greater Trochanter";
                 "Femoral Head"];
             
    path_AL = path_img + "AL_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for j=1:4
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_alg_legend_name = alg_legend_list(n);
            current_data = output_containers_thigh.AL_offset.(current_alg);
            if current_alg == "PTNR"
                scatter(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            else
                plot(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            end
        end
        yL =  ylim;
        fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim(yL);
        grid;
        if j==2
            legend('Location', 'best');
        end
        if j>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(j,2) == 0)
            ylabel("Error (cm)");
        end
        title(AL_name_list(j));
        ax = gca;
        ax.FontSize = 20;
    end
    %suptitle("Anatomical Landmark Offset for Thigh, IDX = " + string(configuration_idx), 'FontSize', 20);
    filename_AL = path_AL + "AL_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_AL);
    close(gcf);
    
    path_ClusterPt = path_img + "ClusterPt_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for j=1:N
        %subplot(2, N/2, j);
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_alg_legend_name = alg_legend_list(n);
            current_data = output_containers_thigh.cluster_offset.(current_alg);
            if current_alg == "PTNR"
                scatter(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            else
                plot(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            end
        end
        yL =  ylim;
        fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim(yL);
        grid;
        if j==3
            legend('Location', 'best');
        end
        if j>3
            xlabel("Fraction completed gait cycle");
        end
        if or(j==1, j==4)
            ylabel("Error (cm)");
        end
        title("Cluster Point " + string(j));
        ax = gca;
        ax.FontSize = 20;
    end
    %suptitle("Cluster Point Global Position Reconstruction Error for Thigh, IDX = " + string(configuration_idx), 'FontSize', 20);
    filename_ClusterPt = path_ClusterPt + "ClusterPt_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_ClusterPt);
    close(gcf);

end

function generate_output_plots_shank(input_struct, output_containers_shank)
    input_size = input_struct.input_size;
    input_data = input_struct.input_data_shank;

    configuration_idx = input_data.configuration_idx;
    N = input_size.N;
    fraction_completed_gait_cycle = input_size.fraction_completed_gait_cycle;
    ToeOff = 0.6;
    
    color_info = input_struct.color_info;
    alg_list = ["SVDLS", "PTUR", "PCTO"];
    alg_legend_list = ["SVD-LS", "PCT-PT", "PCT"];

    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N) + "\Shank\";
    
    path_CM = path_img + "CM_plots\";
    figure('units','normalized','outerposition',[0 0 1 1])
    xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
    hold on;
    for n=1:length(alg_list)
        current_alg = alg_list(n);
        current_alg_legend_name = alg_legend_list(n);
        current_data = output_containers_shank.CM_offset.(current_alg);
        if current_alg == "PTNR"
            scatter(fraction_completed_gait_cycle, current_data, 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
        else
            plot(fraction_completed_gait_cycle, current_data, 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
        end
    end
    yL =  ylim;
    fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim(yL);
    grid;
    legend('Location', 'best');
    xlabel("Fraction completed gait cycle");
    ylabel("Error (cm)");
    %title("Center of Mass offset for Shank, IDX = " + string(configuration_idx));
    ax = gca;
    ax.FontSize = 20;
    filename_CM = path_CM + "CM_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_CM);
    close(gcf);

    path_Eigenvector = path_img + "Eigenvector_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(3,1, 'TileSpacing', 'compact');
    for j=1:3
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_alg_legend_name = alg_legend_list(n);
            current_data = output_containers_shank.eigenvector_offset.(current_alg);
            plot(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
        end
        yL =  ylim;
        fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim(yL);
        grid;
        if j==1
            legend('Location', 'best');
        end
        if j==3
            xlabel("Fraction completed gait cycle");
        end
        ylabel("Error (degrees)");
        title("Eigenvector " + string(j));
        ax = gca;
        ax.FontSize = 20;
    end
    %suptitle("Eigenvector Offset for Shank, IDX = " + string(configuration_idx), 'FontSize', 20);
    filename_Eigenvector = path_Eigenvector + "Eigenvector_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_Eigenvector);
    close(gcf);

    AL_name_list = ["Lateral Malleolus";
                 "Medial Malleolus";
                 "Fibular Head";
                 "Tibial Tuberosity"];
             
    path_AL = path_img + "AL_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for j=1:4
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_alg_legend_name = alg_legend_list(n);
            current_data = output_containers_shank.AL_offset.(current_alg);
            if current_alg == "PTNR"
                scatter(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            else
                plot(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            end
        end
        yL =  ylim;
        fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim(yL);
        grid;
        if j==2
            legend('Location', 'best');
        end
        if j>2
            xlabel("Fraction completed gait cycle");
        end
        if ~(mod(j,2) == 0)
            ylabel("Error (cm)");
        end
        title(AL_name_list(j));
        ax = gca;
        ax.FontSize = 20;
    end
    %suptitle("Anatomical Landmark Offset for Shank, IDX = " + string(configuration_idx), 'FontSize', 20);
    filename_AL = path_AL + "AL_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_AL);
    close(gcf);
    
    path_ClusterPt = path_img + "ClusterPt_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,2, 'TileSpacing', 'compact');
    for j=1:N
        nexttile;
        xline(ToeOff, 'k--', 'DisplayName', 'Toe Off');
        hold on;
        for n=1:length(alg_list)
            current_alg = alg_list(n);
            current_alg_legend_name = alg_legend_list(n);
            current_data = output_containers_shank.cluster_offset.(current_alg);
            if current_alg == "PTNR"
                scatter(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            else
                plot(fraction_completed_gait_cycle, current_data(:,j), 'Color', color_info.(current_alg), 'LineWidth', 1, 'DisplayName', current_alg_legend_name);
            end
        end
        yL =  ylim;
        fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
        ylim(yL);
        grid;
        if j==3
            legend('Location', 'best');
        end
        if j>3
            xlabel("Fraction completed gait cycle");
        end
        if or(j==1, j==4)
            ylabel("Error (cm)");
        end
        title("Cluster Point " + string(j));
        ax = gca;
        ax.FontSize = 20;
    end
    %suptitle("Cluster Point Global Position Reconstruction Error for Shank, IDX = " + string(configuration_idx), 'FontSize', 20);
    filename_ClusterPt = path_ClusterPt + "ClusterPt_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_ClusterPt);
    close(gcf);

end

function AFoffset = generate_AL_origin_offset_plot(output_container_current_PCT, input_struct)
    input_size = input_struct.input_size;
    N = input_size.N;
    configuration_idx = input_struct.input_data_thigh.configuration_idx;
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N) + "\";

    thigh_data = struct;
    thigh_data.REF = output_container_current_PCT.thigh.AF_Translation_reference.SVDLS;
    thigh_data.SVDLS = output_container_current_PCT.thigh.AF_Translation_estimated.SVDLS - thigh_data.REF;
    thigh_data.PTUR = output_container_current_PCT.thigh.AF_Translation_estimated.PTUR - thigh_data.REF;
    thigh_data.PCTO = output_container_current_PCT.thigh.AF_Translation_estimated.PCTO - thigh_data.REF;
    
    thigh_offset = struct;
    thigh_offset.SVDLS = vecnorm(thigh_data.SVDLS,2,2);
    thigh_offset.PTUR = vecnorm(thigh_data.PTUR,2,2);
    thigh_offset.PCTO = vecnorm(thigh_data.PCTO,2,2);
    
    shank_data = struct;
    shank_data.REF = output_container_current_PCT.shank.AF_Translation_reference.SVDLS;
    shank_data.SVDLS = output_container_current_PCT.shank.AF_Translation_estimated.SVDLS - shank_data.REF;
    shank_data.PTUR = output_container_current_PCT.shank.AF_Translation_estimated.PTUR - shank_data.REF;
    shank_data.PCTO = output_container_current_PCT.shank.AF_Translation_estimated.PCTO - shank_data.REF;
    
    shank_offset = struct;
    shank_offset.SVDLS = vecnorm(shank_data.SVDLS,2,2);
    shank_offset.PTUR = vecnorm(shank_data.PTUR,2,2);
    shank_offset.PCTO = vecnorm(shank_data.PCTO,2,2);
    
    AFoffset = struct;
    AFoffset.thigh = thigh_offset;
    AFoffset.shank = shank_offset;
    
    color_info = input_struct.color_info;

    path_AFoffset = path_img + "AFoffset_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(2,1, 'TileSpacing', 'compact');
    
    nexttile;
    xline(0.6, 'k--', 'DisplayName', 'ToeOff');
    hold on;
    plot(input_size.fraction_completed_gait_cycle, thigh_offset.SVDLS, 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
    plot(input_size.fraction_completed_gait_cycle, thigh_offset.PTUR, 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
    plot(input_size.fraction_completed_gait_cycle, thigh_offset.PCTO, 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
    yL =  ylim;
    fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim(yL);
    grid;
    legend('Location', 'best');
    ylabel("Error (cm)")
    title("Thigh");
    ax = gca;
    ax.FontSize = 20;
    
    nexttile;
    xline(0.6, 'k--', 'DisplayName', 'ToeOff');
    hold on;
    plot(input_size.fraction_completed_gait_cycle, shank_offset.SVDLS, 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
    plot(input_size.fraction_completed_gait_cycle, shank_offset.PTUR, 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
    plot(input_size.fraction_completed_gait_cycle, shank_offset.PCTO, 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
    yL =  ylim;
    fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim(yL);
    grid;
    legend('Location', 'best');
    ylabel("Error (cm)")
    xlabel("Fraction completed gait cycle");
    title("Shank", 'FontSize', 20);
    ax = gca;
    ax.FontSize = 20;
    
    %suptitle("Anatomical Frame offset for IDX = " + string(configuration_idx), 'FontSize', 20);
    
    filename_AFoffset = path_AFoffset + "AFoffset_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_AFoffset);
    close(gcf);
end

function generate_knee_angle_plot(knee_angle_current, input_struct)
    input_size = input_struct.input_size;
    N = input_size.N;
    configuration_idx = input_struct.input_data_thigh.configuration_idx;
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(N) + "\";
    
    path_KneeAngle = path_img + "KneeAngle_plots\";
    figure('units','normalized','outerposition',[0 0 1 1]);
    tiledlayout(3,1, 'TileSpacing', 'compact');
    
    color_info = input_struct.color_info;
    
    nexttile;
    xline(0.6, 'k--', 'DisplayName', 'ToeOff');
    hold on;
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.reference.flexion, color_info.reference, 'LineWidth', 1, 'DisplayName', 'Reference');
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.SVDLS.flexion, 'Color', color_info.SVDLS, 'LineWidth', 1, 'DisplayName', 'SVD-LS');
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.PTUR.flexion, 'Color', color_info.PTUR, 'LineWidth', 1, 'DisplayName', 'PCT-PT');
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.PCTO.flexion, 'Color', color_info.PCTO, 'LineWidth', 1, 'DisplayName', 'PCT');
    yL =  ylim;
    fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim(yL);
    
    ylabel("degrees");
    grid;
    legend('Location', 'best');
    title("Flexion Extension");
    ax = gca;
    ax.FontSize = 20;
    
    nexttile;
    xline(0.6, 'k--');
    hold on;
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.reference.adduction, color_info.reference, 'LineWidth', 1);
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.SVDLS.adduction, 'Color', color_info.SVDLS, 'LineWidth', 1);
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.PTUR.adduction, 'Color', color_info.PTUR, 'LineWidth', 1);
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.PCTO.adduction, 'Color', color_info.PCTO, 'LineWidth', 1);
    yL =  ylim;
    fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim(yL);
    ylabel("degrees");
    grid;
    title("Adduction");
    ax = gca;
    ax.FontSize = 20;
    
    nexttile;
    xline(0.6, 'k--');
    hold on;
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.reference.externalRotation, color_info.reference, 'LineWidth', 1);
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.SVDLS.externalRotation, 'Color', color_info.SVDLS, 'LineWidth', 1);
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.PTUR.externalRotation, 'Color', color_info.PTUR, 'LineWidth', 1);
    plot(input_size.fraction_completed_gait_cycle, knee_angle_current.estimated.PCTO.externalRotation, 'Color', color_info.PCTO, 'LineWidth', 1);
    yL =  ylim;
    fill([0 0.6 0.6 0], [yL(1) yL(1) yL(2) yL(2)], 'k', 'FaceAlpha', 0.05, 'HandleVisibility', 'off');
    ylim(yL);
    ylabel("degrees");
    xlabel("fraction completed gait cycle");
    grid;
    title("External Rotation");
    ax = gca;
    ax.FontSize = 20;
    
    %suptitle("Knee Angles for IDX = " + string(configuration_idx), 'FontSize', 20);
    
    ax = gca;
    ax.FontSize = 20;
    
    filename_KneeAngle = path_KneeAngle + "KneeAngle_plot_" + string(configuration_idx) + ".png";
    exportgraphics(gcf, filename_KneeAngle);
    close(gcf);
    
end

function knee_angle = estimate_knee_angle_gs(output_container_current, input_size, measurement_type, algorithm_str)
    knee_angle = struct;
    current_side = "left";
    
    NumMeasurements = input_size.NumFrames;
    
    flexion_list = nan(NumMeasurements,1);
    adduction_list = nan(NumMeasurements,1);
    externalRotation_list = nan(NumMeasurements,1);
    
    for n=1:NumMeasurements
        disp(n);
        
        measurement_str = "AF_Rotation_" + measurement_type;
        output_data = struct;
        output_data.thigh = output_container_current.thigh.(measurement_str).(algorithm_str);
        output_data.shank = output_container_current.shank.(measurement_str).(algorithm_str);
        
        unit_vector_hat_val = struct;
        unit_vector_hat_val.thigh = generate_unit_vectors_local(output_data.thigh(n,:,:));
        unit_vector_hat_val.shank = generate_unit_vectors_local(output_data.shank(n,:,:));

        i = unit_vector_hat_val.shank.x;
        j = unit_vector_hat_val.shank.y;
        k = unit_vector_hat_val.shank.z;

        I = unit_vector_hat_val.thigh.x;
        J = unit_vector_hat_val.thigh.y;
        K = unit_vector_hat_val.thigh.z;

        I_gs = K;
        J_gs = I;
        K_gs = J;

        i_gs = k;
        j_gs = i;
        k_gs = j;

        e1 = I_gs;
        e3 = k_gs;
        e2 = cross(e3,e1);

        e1_hat = e1/norm(e1);
        e2_hat = e2/norm(e2);
        e3_hat = e3/norm(e3);

        alpha = rad2deg(asin(-dot(e2_hat,K_gs)));
        beta = rad2deg(acos(dot(I_gs,k_gs)));
        gamma = rad2deg(asin(dot(e2_hat,i_gs)));

        current_flexion = alpha;
        if current_side == "right"
            current_adduction = beta - 90;
            current_externalRotation = -gamma;
        else
            current_adduction = 90 - beta;
            current_externalRotation = gamma;
        end

        flexion_list(n) = current_flexion;
        adduction_list(n) = current_adduction;
        externalRotation_list(n) = current_externalRotation;
    end
    knee_angle.flexion = flexion_list;
    knee_angle.adduction = adduction_list;
    knee_angle.externalRotation = externalRotation_list;
    
end

function unit_vectors_local = generate_unit_vectors_local(R_input)
    R = reshape(R_input, [3,3]);
    
    unit_vectors_local = struct;
    unit_vectors_local.x = R(:,1);
    unit_vectors_local.y = R(:,2);
    unit_vectors_local.z = R(:,3);
end

function input_size = generate_input_size_struct(N, NumConfigurations, box_markers_gait)
    input_size = struct;
    
    fraction_completed_gait_cycle = box_markers_gait(1).fraction_completed_gait_cycle;
    NumFrames = length(fraction_completed_gait_cycle);
    Num_AL_Shank = 4;
    
    input_size.N = N;
    input_size.fraction_completed_gait_cycle = fraction_completed_gait_cycle;
    input_size.NumFrames = NumFrames;
    input_size.Num_AL_Shank = Num_AL_Shank;
    input_size.NumAngles = 3;
    input_size.NumConfigurations = NumConfigurations;
end

function create_animation_RF(input_struct, output_container_current, reference_frame_str)
    
    path_img = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\PCT_plots\N=" + string(input_struct.input_size.N) + "\";
    
    movingView_flag = input_struct.visualization_settings.movingView_flag_val;
    animation_flag = input_struct.visualization_settings.animation_flag;
    fig_scale = input_struct.visualization_settings.fig_scale;
    
    substruct_str = struct;
    
    if reference_frame_str == "TF"
        
        substruct_str.R_estimated =  "R_estimated";
        substruct_str.T_estimated =  "T_estimated";
        
        substruct_str.R_reference =  "R_reference";
        substruct_str.T_reference =  "T_reference";
        
        folder_name_start = "Movement_Videos_PCT_TF\movement_subNo_";
        if movingView_flag
            folder_name_start = "Movement_MovingPtView_Videos_PCT_TF\movement_subNo_";
        end
        folder_name_start = path_img + folder_name_start;
        
    elseif reference_frame_str == "AF"
        
        substruct_str.R_estimated = "AF_Rotation_estimated";
        substruct_str.T_estimated = "AF_Translation_estimated";
        
        substruct_str.R_reference = "AF_Rotation_reference";
        substruct_str.T_reference = "AF_Translation_reference";
        
        folder_name_start = "Movement_Videos_PCT_AF\movement_subNo_";
        if movingView_flag
            folder_name_start = "Movement_MovingPtView_Videos_PCT_AF\movement_subNo_";
        end
        folder_name_start = path_img + folder_name_start;
        
    end
    
    if animation_flag
        
        input_data_thigh = input_struct.input_data_thigh;
        pose_gait_thigh = input_data_thigh.pose_gait_thigh;
        cylinder_data_thigh = input_data_thigh.cylinder_data_thigh;
        mass_distribution_reference_standard_pose_thigh = input_data_thigh.mass_distribution_reference_thigh;
        box_markers_gait_thigh = input_data_thigh.box_markers_gait_thigh;
        
        input_data_shank = input_struct.input_data_shank;
        pose_gait_shank = input_data_shank.pose_gait_shank;
        cylinder_data_shank = input_data_shank.cylinder_data_shank;
        mass_distribution_reference_standard_pose_shank = input_data_shank.mass_distribution_reference_shank;
        box_markers_gait_shank = input_data_shank.box_markers_gait_shank;
        
        configuration_idx = input_data_shank.configuration_idx;
        
        input_size = input_struct.input_size;
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
            
            unit_vector_local.thigh = struct;
            unit_vector_local.thigh.REF = generate_unit_vectors_local(output_container_current.thigh.(substruct_str.R_reference).SVDLS(k, :, :));
            unit_vector_local.thigh.SVDLS = generate_unit_vectors_local(output_container_current.thigh.(substruct_str.R_estimated).SVDLS(k, :, :));
            unit_vector_local.thigh.PTUR = generate_unit_vectors_local(output_container_current.thigh.(substruct_str.R_estimated).PTUR(k, :, :));
            
            unit_vector_local.shank = struct;
            unit_vector_local.shank.REF = generate_unit_vectors_local(output_container_current.shank.(substruct_str.R_reference).SVDLS(k, :, :));
            unit_vector_local.shank.SVDLS = generate_unit_vectors_local(output_container_current.shank.(substruct_str.R_estimated).SVDLS(k, :, :));
            unit_vector_local.shank.PTUR = generate_unit_vectors_local(output_container_current.shank.(substruct_str.R_estimated).PTUR(k, :, :));
            
            current_pose_thigh = pose_gait_thigh(k);
            cylinder_data_pose_thigh = generate_cylinder_data_pose(cylinder_data_thigh, current_pose_thigh);
            X_thigh = cylinder_data_pose_thigh.X;
            Y_thigh = cylinder_data_pose_thigh.Y;
            Z_thigh = cylinder_data_pose_thigh.Z;
            
            current_pose_shank = pose_gait_shank(k);
            cylinder_data_pose_shank = generate_cylinder_data_pose(cylinder_data_shank, current_pose_shank);
            X_shank = cylinder_data_pose_shank.X;
            Y_shank = cylinder_data_pose_shank.Y;
            Z_shank = cylinder_data_pose_shank.Z;
            
            mass_distribution_noisy_standard_pose_thigh = generate_noisy_distribution_STA(mass_distribution_reference_standard_pose_thigh, k, box_markers_gait_thigh, N);
            mass_distribution_reference_thigh = apply_current_pose(mass_distribution_reference_standard_pose_thigh, current_pose_thigh, N);
            mass_distribution_noisy_thigh = apply_current_pose(mass_distribution_noisy_standard_pose_thigh, current_pose_thigh, N);
            
            mass_distribution_noisy_standard_pose_shank = generate_noisy_distribution_STA(mass_distribution_reference_standard_pose_shank, k, box_markers_gait_shank, N);
            mass_distribution_reference_shank = apply_current_pose(mass_distribution_reference_standard_pose_shank, current_pose_shank, N);
            mass_distribution_noisy_shank = apply_current_pose(mass_distribution_noisy_standard_pose_shank, current_pose_shank, N);
            
            clf
            global_frame_plotter;
            surf(X_thigh, Y_thigh, Z_thigh, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', '0.1');
            hold on;
            surf(X_shank, Y_shank, Z_shank, 'EdgeColor', 'none', 'FaceColor', 'b', 'FaceAlpha', '0.1');
            
            technical_frame_plotter(output_container_current.thigh.(substruct_str.T_reference).SVDLS(k, :, 1), unit_vector_local.thigh.REF, "REF");
            technical_frame_plotter(output_container_current.thigh.(substruct_str.T_estimated).SVDLS(k, :, 1), unit_vector_local.thigh.SVDLS, "SVDLS");
            technical_frame_plotter(output_container_current.thigh.(substruct_str.T_estimated).PTUR(k, :, 1), unit_vector_local.thigh.PTUR, "PTUR");
            
            technical_frame_plotter(output_container_current.shank.(substruct_str.T_reference).SVDLS(k, :, 1), unit_vector_local.shank.REF, "REF");
            technical_frame_plotter(output_container_current.shank.(substruct_str.T_estimated).SVDLS(k, :, 1), unit_vector_local.shank.SVDLS, "SVDLS");
            technical_frame_plotter(output_container_current.shank.(substruct_str.T_estimated).PTUR(k, :, 1), unit_vector_local.shank.PTUR, "PTUR");
            
            for j=1:N
                scatter3(mass_distribution_reference_thigh(j).X, mass_distribution_reference_thigh(j).Y, mass_distribution_reference_thigh(j).Z, 'k', 'filled');
                scatter3(mass_distribution_noisy_thigh(j).X, mass_distribution_noisy_thigh(j).Y, mass_distribution_noisy_thigh(j).Z, 'm');
                
                scatter3(mass_distribution_reference_shank(j).X, mass_distribution_reference_shank(j).Y, mass_distribution_reference_shank(j).Z, 'k', 'filled');
                scatter3(mass_distribution_noisy_shank(j).X, mass_distribution_noisy_shank(j).Y, mass_distribution_noisy_shank(j).Z, 'm');
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

function cylinder_data_thigh = generate_cylinder_data_thigh
    t = 0:0.1:50;
    n = length(t);
    np = 100;
    r = 5.5 + 0.08*t;
    [X, Y, Z] = cylinder(r, np-1);
    Zscaled = 50*Z;

    cylinder_data_thigh = struct;
    cylinder_data_thigh.X = X;
    cylinder_data_thigh.Y = Y;
    cylinder_data_thigh.Z = Zscaled;
    cylinder_data_thigh.np = np;
    cylinder_data_thigh.n = n;
end

function cylinder_data_shank = generate_cylinder_data_shank
    shank_length = 35;
    t = 0:0.1:shank_length;
    n = length(t);
    np = 100;
    r = 5.5 - 0.0429*t;
    [X, Y, Z] = cylinder(r, np-1);
    Zscaled = -shank_length*Z;

    cylinder_data_shank = struct;
    cylinder_data_shank.X = X;
    cylinder_data_shank.Y = Y;
    cylinder_data_shank.Z = Zscaled;
    cylinder_data_shank.np = np;
    cylinder_data_shank.n = n;
end

function fig_scale_info = generate_fig_scale_info_PCT(fig_scale_flag, input_struct, output_container_current)
    if fig_scale_flag
        plot_edge_width = 50;
        view_vector = [-37.5 30];

        input_size = input_struct.input_size;
        num_measurements = input_size.NumFrames;
        
        cylinder_x_mean_list_thigh = nan(num_measurements,1);
        cylinder_y_mean_list_thigh = nan(num_measurements,1);
        cylinder_z_mean_list_thigh = nan(num_measurements,1);
        
        cylinder_x_mean_list_shank = nan(num_measurements,1);
        cylinder_y_mean_list_shank = nan(num_measurements,1);
        cylinder_z_mean_list_shank = nan(num_measurements,1);
        
        for n=1:num_measurements
            
            cylinder_x_mean_list_thigh(n) = output_container_current.thigh.T_estimated.SVDLS(n,1);
            cylinder_y_mean_list_thigh(n) = output_container_current.thigh.T_estimated.SVDLS(n,2);
            cylinder_z_mean_list_thigh(n) = output_container_current.thigh.T_estimated.SVDLS(n,3);
            
            cylinder_x_mean_list_shank(n) = output_container_current.shank.T_estimated.SVDLS(n,1);
            cylinder_y_mean_list_shank(n) = output_container_current.shank.T_estimated.SVDLS(n,2);
            cylinder_z_mean_list_shank(n) = output_container_current.shank.T_estimated.SVDLS(n,3);

        end

        x_max_limit = max([max(cylinder_x_mean_list_thigh) , max(cylinder_x_mean_list_shank)]) + plot_edge_width;
        y_max_limit = max([max(cylinder_y_mean_list_thigh) , max(cylinder_y_mean_list_shank)]) + plot_edge_width;
        z_max_limit = max([max(cylinder_z_mean_list_thigh) , max(cylinder_z_mean_list_shank)]) + plot_edge_width;

        x_min_limit = min([min(cylinder_x_mean_list_thigh) , min(cylinder_x_mean_list_shank)]) - plot_edge_width;
        y_min_limit = min([min(cylinder_y_mean_list_thigh) , min(cylinder_y_mean_list_shank)]) - plot_edge_width;
        z_min_limit = min([min(cylinder_z_mean_list_thigh) , min(cylinder_z_mean_list_shank)]) - plot_edge_width;

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