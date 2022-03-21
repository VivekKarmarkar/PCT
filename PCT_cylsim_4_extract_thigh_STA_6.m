current_struct_name_idx = 2;
STA_plot_bool = true;
plot_flag = false;

N = 6;
folder_name = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\STA_data\Matlab\";
str_list = ["UNIMELB", "LMAM"];
current_str_idx = current_struct_name_idx;
current_str = str_list(current_str_idx);
struct_name = "dataSample_" + current_str + "_treadmillWalking";
file_name = struct_name + ".mat";
path_name = folder_name + file_name;
load(path_name);

NumTrials = length(dataSample_LMAM_treadmillWalking.subj01.trial01.L_thigh.mrk.M01);
frames = 1:NumTrials;
completionPercent = frames/NumTrials;

HeelStrike = dataSample_LMAM_treadmillWalking.subj01.trial01.evt.HeelStrike;
NumGaitCycle = length(HeelStrike) - 1;
ToeOff = HeelStrike(1:end-1) + floor(0.6*diff(HeelStrike));
Switch = HeelStrike(1:end-1) + floor(0.7*diff(HeelStrike));

lower_limit_cm = 1;
upper_limit_cm = 3;

X = nan(2*NumGaitCycle,1);
X(1:2:end) = HeelStrike(1:end-1);
X(2:2:end) = HeelStrike(1:end-1) + floor(0.5*diff(HeelStrike));
X = [X;HeelStrike(end)];
Y = repmat([lower_limit_cm; upper_limit_cm], 3, 1);
Y = [Y;lower_limit_cm];

NumMarkers = length(fields(dataSample_LMAM_treadmillWalking.subj01.trial01.L_thigh.mrk));
marker_idx_list = 1:NumMarkers;
h = nan(NumMarkers, 1);
theta = nan(NumMarkers, 1);

box_markers = struct;

box_markers(1).name = "Anterior Low";
box_markers(1).idx_bool = false(NumMarkers, 1);

box_markers(2).name = "Anterior Mid";
box_markers(2).idx_bool = false(NumMarkers, 1);

box_markers(3).name = "Anterior High";
box_markers(3).idx_bool = false(NumMarkers, 1);

box_markers(4).name = "Lateral Low";
box_markers(4).idx_bool = false(NumMarkers, 1);

box_markers(5).name = "Lateral Mid";
box_markers(5).idx_bool = false(NumMarkers, 1);

box_markers(6).name = "Lateral High";
box_markers(6).idx_bool = false(NumMarkers, 1);

for j=1:NumMarkers
    marker_str = "M" + num2str(j, '%.2d');
    marker_data = dataSample_LMAM_treadmillWalking.subj01.trial01.L_thigh.mrk.(marker_str);
    x1 = marker_data(1,1);
    y1 = marker_data(1,2);
    z1 = marker_data(1,3);
    h(j) = y1/10;
    theta(j) = eval_theta(x1,z1);
    box_val = eval_marker(theta(j), h(j));
    idx_bool_new = box_markers(box_val).idx_bool;
    idx_bool_new(j) = true;
    box_markers(box_val).idx_bool = idx_bool_new;
end
figure;
xline(120, '--r', 'LineWidth', 1.5);
hold on;
scatter(theta, h);
xlabel("Theta (degrees)");
ylabel("Height (cm)");
title("Marker distribution");
yline(15, '--k');
yline(30, '--k');
xline(60, '--k')
grid;

input_data.marker_idx_list = marker_idx_list;
input_data.box_markers = box_markers;
input_data.dataSample_LMAM_treadmillWalking = dataSample_LMAM_treadmillWalking;
input_data.NumTrials = NumTrials;
input_data.completionPercent = completionPercent;
input_data.plot_flag = plot_flag;

x_offset_total = nan(NumTrials, 6);
y_offset_total = nan(NumTrials, 6);
z_offset_total = nan(NumTrials, 6);
for k=1:N
    box_idx_current = k;
    input_data.box_idx_current = box_idx_current;

    vector_normalized_mean = generate_box_marker_data(input_data);
    
    box_markers(k).vector = vector_normalized_mean;
    box_markers(k).x = vector_normalized_mean(:,1);
    box_markers(k).y = vector_normalized_mean(:,2);
    box_markers(k).z = vector_normalized_mean(:,3);
    
    x_offset_total(:, box_idx_current) = box_markers(k).x;
    y_offset_total(:, box_idx_current) = box_markers(k).y;
    z_offset_total(:, box_idx_current) = box_markers(k).z;
    
end
x_offset_cm = mean(x_offset_total, 2);
y_offset_cm = mean(y_offset_total, 2);
z_offset_cm = mean(z_offset_total, 2);
norm_offset_cm_6 = sqrt(x_offset_cm.^2 + y_offset_cm.^2 + z_offset_cm.^2);

figure
yline(0,'-k', 'DisplayName', 'x axis');
hold on;
plot(x_offset_cm/10, 'r', 'LineWidth', 0.5, 'DisplayName', 'x CM offset');
plot(y_offset_cm/10, 'g', 'LineWidth', 1.5, 'DisplayName', 'y CM offset');
plot(z_offset_cm/10, 'b', 'LineWidth', 0.5, 'DisplayName', 'z CM offset');
for j=1:4
    xline(HeelStrike(j), 'k--', 'LineWidth', 1, 'DisplayName', 'Heel Strike')
end
for j=1:3
    xline(ToeOff(j), 'm--', 'DisplayName', 'Toe Off')
    xline(Switch(j), 'cyan--', 'LineWidth', 1, 'DisplayName', 'Switch')
    patch([ToeOff(j) HeelStrike(j+1) HeelStrike(j+1) ToeOff(j)], [-1 -1 2 2], 'black', 'FaceAlpha', 0.1, 'DisplayName', 'delta CM y negative')
end
grid;
ylim([-1,2]);
xlim([0, HeelStrike(end)]);
xlabel("Frame number");
ylabel("x (red),   y (green),   z (blue)");
title("Comparison of CM offset")

figure
plot(norm_offset_cm_6/10, 'r', 'LineWidth', 1.5, 'DisplayName', 'CM offset norm');
hold on;
for j=1:4
    xline(HeelStrike(j), 'k--', 'DisplayName', 'Heel Strike')
end
grid;
line(X,Y);
xlabel("Frame number");
ylabel("CM norm offset");
title("CM norm offset and Bouunding Line");
ylim([0, upper_limit_cm + 0.5]);
legend;


box_markers_gait = generate_box_markers_gait_cycle_data(box_markers, HeelStrike);
NumMarkers = 6;
NumTrials_gait = length(box_markers_gait(1).x);

ToeOff_gait = 0.6 * NumTrials_gait;
Switch_gait = 0.7 * NumTrials_gait;

upper_bound_adaptive_CM = 1.8;
lower_bound_adaptive_CM = 0.7;
max_pt_adaptive_CM = 0.5;

X_gait = [0; max_pt_adaptive_CM * NumTrials_gait; NumTrials_gait];
Y_gait = [lower_bound_adaptive_CM; upper_bound_adaptive_CM; lower_bound_adaptive_CM];

x_offset_total_gait = nan(NumTrials_gait, NumMarkers);
y_offset_total_gait = nan(NumTrials_gait, NumMarkers);
z_offset_total_gait = nan(NumTrials_gait, NumMarkers);
for k=1:NumMarkers
    box_idx_current = k;
    input_data.box_idx_current = box_idx_current;
    
    x_offset_total_gait(:, box_idx_current) = box_markers_gait(k).x;
    y_offset_total_gait(:, box_idx_current) = box_markers_gait(k).y;
    z_offset_total_gait(:, box_idx_current) = box_markers_gait(k).z;
    
end
x_offset_cm_gait = mean(x_offset_total_gait, 2);
y_offset_cm_gait = mean(y_offset_total_gait, 2);
z_offset_cm_gait = mean(z_offset_total_gait, 2);
norm_offset_cm_6_gait = sqrt(x_offset_cm_gait.^2 + y_offset_cm_gait.^2 + z_offset_cm_gait.^2);

figure
yline(0,'-k', 'DisplayName', 'x axis');
hold on;
plot(x_offset_cm_gait/10, 'r', 'LineWidth', 0.5, 'DisplayName', 'x CM offset');
plot(y_offset_cm_gait/10, 'g', 'LineWidth', 1.5, 'DisplayName', 'y CM offset');
plot(z_offset_cm_gait/10, 'b', 'LineWidth', 0.5, 'DisplayName', 'z CM offset');

xline(ToeOff_gait, 'm--', 'DisplayName', 'Toe Off');
xline(Switch_gait, 'black--', 'LineWidth', 1, 'DisplayName', 'Switch');
patch([ToeOff_gait NumTrials_gait NumTrials_gait ToeOff_gait], [-1 -1 2 2], 'black', 'FaceAlpha', 0.1, 'DisplayName', 'delta CM y negative');

grid;
ylim([-1,2]);
xlabel("Frame number");
ylabel("x (red),   y (green),   z (blue)");
title("Comparison of CM offset (Average Gait Cycle)")

figure
plot(norm_offset_cm_6_gait/10, 'r', 'LineWidth', 1.5, 'DisplayName', 'CM offset norm');
hold on;
grid;
line(X_gait, Y_gait);
xline(ToeOff_gait, 'm--', 'DisplayName', 'Toe Off');
xline(Switch_gait, 'black--', 'LineWidth', 1, 'DisplayName', 'Switch');
xlabel("Frame number");
ylabel("CM norm offset");
title("CM norm offset and Bounding Line (Average Gait Cycle)");
ylim([0, upper_bound_adaptive_CM + 0.5]);
legend;

generate_STA_plots(box_markers_gait, STA_plot_bool);

save('box_markers_gait_6.mat', 'box_markers_gait');

function generate_STA_plots(box_markers_gait, STA_plot_bool)
    if STA_plot_bool
        frames_current_array = 1:length(box_markers_gait(1).x);
        frames_current_normalized = frames_current_array - frames_current_array(1);
        frames_current_length = length(frames_current_normalized) - 1;
        fraction_completed = frames_current_normalized/frames_current_length;
        
        ToeOff_gait = 0.6;
        NumTrials_gait = 1.0;
        
        NumMarkers = length(box_markers_gait);
        for j=1:NumMarkers
            figure
            
            yline(0, '--k', 'x axis');
            hold on;
            
            plot(fraction_completed, box_markers_gait(j).x/10, 'r', 'LineWidth', 1, 'DisplayName', 'x');
            plot(fraction_completed, box_markers_gait(j).y/10, 'g', 'LineWidth', 1.5, 'DisplayName', 'y');
            plot(fraction_completed, box_markers_gait(j).z/10, 'b', 'LineWidth', 1, 'DisplayName', 'z');
            
            xline(ToeOff_gait, 'm--', 'LineWidth', 1.5, 'DisplayName', 'Toe Off');
            
            xlabel("Fraction of Gait cycle completed");
            ylabel("Noise (cm)");
            title(box_markers_gait(j).name + " STA  (Average for Gait cycle)");
            legend;
            grid;
        end
    end
end

function box_markers_gait = generate_box_markers_gait_cycle_data(box_markers, HeelStrike)
    NumBoxes = length(box_markers);
    
    box_markers_gait = struct;
    
    if HeelStrike(1) == 0
        HeelStrike(1) = 1;
    end
    
    for k=1:NumBoxes
        box_marker_idx = k;
        [vector_gait, ~, fraction_completed_union] = generate_average_gait_cycle_data(box_markers, box_marker_idx, HeelStrike);
        
        box_markers_gait(box_marker_idx).name = box_markers(box_marker_idx).name;
        box_markers_gait(box_marker_idx).fraction_completed_gait_cycle = fraction_completed_union;
        box_markers_gait(box_marker_idx).vector = vector_gait;
        box_markers_gait(box_marker_idx).x = vector_gait(:,1);
        box_markers_gait(box_marker_idx).y = vector_gait(:,2);
        box_markers_gait(box_marker_idx).z = vector_gait(:,3);
    end
end

function [vector_gait, gait_cycle_data, fraction_completed_union] = generate_average_gait_cycle_data(box_markers, box_marker_idx, HeelStrike)
    gait_cycle_data = struct;
    fraction_completed_union = 0;
    
    if HeelStrike(1) == 0
        HeelStrike(1) = 1;
    end
    
    NumGaitCycles = length(HeelStrike) - 1;
    for j=1:NumGaitCycles
        frames_current_array = HeelStrike(j):HeelStrike(j+1);
        frames_current_normalized = frames_current_array - frames_current_array(1);
        frames_current_length = length(frames_current_normalized) - 1;
        
        gait_cycle_data(j).frames = frames_current_array;
        gait_cycle_data(j).fraction_completed = frames_current_normalized/frames_current_length;
        gait_cycle_data(j).offset_vector = box_markers(box_marker_idx).vector(frames_current_array, :);
        gait_cycle_data(j).vector = gait_cycle_data(j).offset_vector - gait_cycle_data(j).offset_vector(1,:);
        
        fraction_completed_union = union(fraction_completed_union, gait_cycle_data(j).fraction_completed);
        
    end
    
    NumDataPts = length(fraction_completed_union);
    x_all = nan(NumDataPts, NumGaitCycles);
    y_all = nan(NumDataPts, NumGaitCycles);
    z_all = nan(NumDataPts, NumGaitCycles);
    for j=1:NumGaitCycles
        fraction_completed_current = gait_cycle_data(j).fraction_completed;
        x_current = gait_cycle_data(j).vector(:,1);
        y_current = gait_cycle_data(j).vector(:,2);
        z_current = gait_cycle_data(j).vector(:,3);
        
        x_interpolated = interp1(fraction_completed_current, x_current, fraction_completed_union)';
        y_interpolated = interp1(fraction_completed_current, y_current, fraction_completed_union)';
        z_interpolated = interp1(fraction_completed_current, z_current, fraction_completed_union)';
        vector_interpolated = horzcat(x_interpolated, y_interpolated, z_interpolated);
        
        x_all(:,j) = x_interpolated;
        y_all(:,j) = y_interpolated;
        z_all(:,j) = z_interpolated;
        
        gait_cycle_data(j).vector_interpolated = vector_interpolated;
        gait_cycle_data(j).fraction_completed_union = fraction_completed_union;
    end
    
    x_gait = mean(x_all, 2);
    y_gait = mean(y_all, 2);
    z_gait = mean(z_all, 2);
    vector_gait = horzcat(x_gait, y_gait, z_gait);
    
end

function vector_normalized_mean = generate_box_marker_data(input_data)
    
    marker_idx_list = input_data.marker_idx_list;
    box_idx_current = input_data.box_idx_current;
    box_markers = input_data.box_markers;
    dataSample_LMAM_treadmillWalking = input_data.dataSample_LMAM_treadmillWalking;
   
    NumTrials = input_data.NumTrials;
    completionPercent = input_data.completionPercent;
    plot_flag = input_data.plot_flag;
    
    box_markers_current = marker_idx_list(box_markers(box_idx_current).idx_bool);
    NumBoxMarkers_current = length(box_markers_current);
    
    HeelStrike = dataSample_LMAM_treadmillWalking.subj01.trial01.evt.HeelStrike;
    NumGaitCycles = length(HeelStrike);
    
    HeelStrike_1 = HeelStrike(1);
    ToeOff_1 = floor(0.6*dataSample_LMAM_treadmillWalking.subj01.trial01.evt.HeelStrike(2));

    plot_vert_str = ["x (cm)", "y (cm)", "z (cm)"];
    if plot_flag
        figure
        for k=1:3
            h(k) = subplot(3,1,k);
            yline(0, 'k--');
            hold on;
            for j=1:length(HeelStrike)
                xline(HeelStrike(j), 'k--', 'LineWidth', 1);
            end
            grid;
        end
    end

    x = nan(NumTrials, NumBoxMarkers_current);
    y = nan(NumTrials, NumBoxMarkers_current);
    z = nan(NumTrials, NumBoxMarkers_current);
    x_normalized = nan(NumTrials, NumBoxMarkers_current);
    y_normalized = nan(NumTrials, NumBoxMarkers_current);
    z_normalized = nan(NumTrials, NumBoxMarkers_current);
    for j=1:NumBoxMarkers_current
        marker_idx = box_markers_current(j);
        marker_str = "M" + num2str(marker_idx, '%.2d');
        marker_data = dataSample_LMAM_treadmillWalking.subj01.trial01.L_thigh.mrk.(marker_str);
        x(:,j) = marker_data(:,1);
        y(:,j) = marker_data(:,2);
        z(:,j) = marker_data(:,3);
        x_normalized(:,j) = marker_data(:,1) - marker_data(1,1);
        y_normalized(:,j) = marker_data(:,2) - marker_data(1,2);
        z_normalized(:,j) = marker_data(:,3) - marker_data(1,3);
        vector_normalized = horzcat(x_normalized(:,j), y_normalized(:,j), z_normalized(:,j));
        
        x_current = x_normalized(:,j);
        y_current = y_normalized(:,j);
        z_current = z_normalized(:,j);
        x_truncated = x_current(1:ToeOff_1);
        y_truncated = y_current(1:ToeOff_1);
        z_truncated = z_current(1:ToeOff_1);
        vector_truncated = horzcat(x_truncated, y_truncated, z_truncated);
        
        if plot_flag
            for k=1:3
                h(k) = subplot(3,1,k);
                hLine(k) = plot(vector_normalized(:,k)/10, 'b');
                %hLine(k) = plot(completionPercent, vector_normalized(:,k), 'b');
                %hLine(k) = plot(vector_truncated(:,k), 'b');
                ylabel(plot_vert_str(k));
            end
        end

    end

    x_normalized_mean = mean(x_normalized, 2);
    y_normalized_mean = mean(y_normalized, 2);
    z_normalized_mean = mean(z_normalized, 2);
    vector_normalized_mean = horzcat(x_normalized_mean, y_normalized_mean, z_normalized_mean);
    
    x_truncated_mean = mean(x_truncated, 2);
    y_truncated_mean = mean(y_truncated, 2);
    z_truncated_mean = mean(z_truncated, 2);
    vector_truncated_mean = horzcat(x_truncated_mean, y_truncated_mean, z_truncated_mean);
    
    if plot_flag
        for k=1:3
            h(k) = subplot(3,1,k);
            hLine(k) = plot(vector_normalized_mean(:,k)/10, 'r', 'LineWidth', 1.5);
            %hLine(k) = plot(completionPercent, vector_normalized_mean(:,k), 'r', 'LineWidth', 1.5);
            %hLine(k) = plot(vector_truncated_mean(:,k), 'r', 'LineWidth', 1.5);
            if k==3
                xlabel("Frame number");
                %xlabel("fraction of task completed");
            end
        end
        suptitle(box_markers(box_idx_current).name + " STA  [Red = Average , Blue = Markers]");
    end

end

function theta_val = eval_theta(x,z)
    z_reflected = -z;

    if x>0 && z_reflected>0
        theta_val = rad2deg(atan(z_reflected/x));
    elseif x<0 && z_reflected>0
        theta_val = 180 + rad2deg(atan(z_reflected/x));
    elseif x<0 && z_reflected<0
        theta_val = 180 + rad2deg(atan(z_reflected/x));
    else
        theta_val = 360 + rad2deg(atan(z_reflected/x));
    end
end

function box_val = eval_marker(theta, h)
    if theta>0 && theta<60
        if h<15
            box_val = 1;
        elseif h>15 && h<30
            box_val = 2;
        else
            box_val = 3;
        end
    else
        if h<15
            box_val = 4;
        elseif h>15 && h<30
            box_val = 5;
        else
            box_val = 6;
        end
    end
end