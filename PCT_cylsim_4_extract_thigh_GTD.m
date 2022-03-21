folder_name_GTD = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\STA_data\Matlab\";
file_name_GTD = "dataSample_LMAM_treadmillWalking.mat";
path_name_GTD = folder_name_GTD + file_name_GTD;
load(path_name_GTD);

[pose_gait, gait_cycle_data, fraction_completed_gait_cycle] = generate_average_gait_cycle_data(dataSample_LMAM_treadmillWalking);
save('ground_truth_data_thigh.mat', 'pose_gait');

function [pose_gait, gait_cycle_data, fraction_completed_union] = generate_average_gait_cycle_data(dataSample_LMAM_treadmillWalking)
    pose_gait = struct;
    gait_cycle_data = struct;
    fraction_completed_union = 0;
    
    HeelStrike = dataSample_LMAM_treadmillWalking.subj01.trial01.evt.HeelStrike;
    
    if HeelStrike(1) == 0
        HeelStrike(1) = 1;
    end
    
    T = dataSample_LMAM_treadmillWalking.subj01.trial01.L_thigh.gta;
    size_T = size(T);
    NumColsT = size_T(2);
    
    R = dataSample_LMAM_treadmillWalking.subj01.trial01.L_thigh.gRa;
    size_R = size(R);
    NumColsR = size_R(2);
    
    NumGaitCycles = length(HeelStrike) - 1;
    for j=1:NumGaitCycles
        frames_current_array = HeelStrike(j):HeelStrike(j+1);
        frames_current_normalized = frames_current_array - frames_current_array(1);
        frames_current_length = length(frames_current_normalized) - 1;
        
        gait_cycle_data(j).frames = frames_current_array;
        gait_cycle_data(j).fraction_completed = frames_current_normalized/frames_current_length;
        gait_cycle_data(j).vector_T = T(frames_current_array, :);
        gait_cycle_data(j).vector_R = R(frames_current_array, :);
        
        fraction_completed_union = union(fraction_completed_union, gait_cycle_data(j).fraction_completed);
        
    end
    
    fraction_completed_union = transpose(fraction_completed_union);
    NumDataPts = length(fraction_completed_union);
    
    T_gait = nan(NumDataPts, NumColsT);
    T_all = nan(NumDataPts, NumGaitCycles, NumColsT);
    
    R_gait = nan(NumDataPts, NumColsR);
    R_all = nan(NumDataPts, NumGaitCycles, NumColsR);
   
    for j=1:NumGaitCycles
        fraction_completed_current = gait_cycle_data(j).fraction_completed;
        
        vector_interpolated_T = nan(NumDataPts, NumColsT);
        for n=1:NumColsT
            col_current_T = gait_cycle_data(j).vector_T(:,n);
            col_interpolated_T = interp1(fraction_completed_current, col_current_T, fraction_completed_union)';
            vector_interpolated_T(:,n) = col_interpolated_T;
            T_all(:,j,n) = vector_interpolated_T(:, n);
        end
        
        vector_interpolated_R = nan(NumDataPts, NumColsR);
        for n=1:NumColsR
            col_current_R = gait_cycle_data(j).vector_R(:,n);
            col_interpolated_R = interp1(fraction_completed_current, col_current_R, fraction_completed_union)';
            vector_interpolated_R(:,n) = col_interpolated_R;
            R_all(:,j,n) = vector_interpolated_R(:, n);
        end
        
        gait_cycle_data(j).vector_interpolated_T = vector_interpolated_T;
        gait_cycle_data(j).vector_interpolated_R = vector_interpolated_R;
        gait_cycle_data(j).fraction_completed_union = fraction_completed_union;
    end
    
    for n=1:NumColsT
        T_current = T_all(:,:,n);
        T_gait(:,n) = mean(T_current, 2);
    end
    
    for n=1:NumColsR
        R_current = R_all(:,:,n);
        R_gait(:,n) = mean(R_current, 2);
    end
    
    for j=1:NumDataPts
        pose_gait(j).T = T_gait(j,:)';
        pose_gait(j).R = [R_gait(j,1:3)', R_gait(j,4:6)', R_gait(j,7:9)'];
    end
    
end