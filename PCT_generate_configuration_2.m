N = 4;
NumTrials = 300000;

folderName_settings = "C:\Users\vkarmarkar\OneDrive - University of Iowa\Desktop\Research\Code\Test Code\MarkerConfigurationConstrainedSettings_" + string(N);

cylinder_data = generate_cylinder_data;
X = cylinder_data.X;
Y = cylinder_data.Y;
Z = cylinder_data.Z;

np_constrained = 34;
X_constrained = X(:,1:np_constrained);
Y_constrained = Y(:,1:np_constrained);
Z_constrained = Z(:,1:np_constrained);

theta = rad2deg(atan(Y./X));
first_jump_idx = find(Y(1,:)==max(Y(1,:)));
second_jump_idx = find(Y(1,:)==min(Y(1,:)));
theta(:,first_jump_idx:second_jump_idx) = theta(:,first_jump_idx:second_jump_idx) + 180;
theta(:,second_jump_idx:end) = theta(:,second_jump_idx:end) + 360;

configuration_count = 0;

for j=1:NumTrials
    disp(j);
    mass_distribution_reference = generate_random_distribution_constrained(cylinder_data, N);
    intermarker_distance_bool = check_intermarker_distance(mass_distribution_reference, N);
    if intermarker_distance_bool
          configuration_count = configuration_count + 1;
          fileName_settings = folderName_settings + "\FixedMarkerDistribution_" + string(configuration_count);
          save(fileName_settings, 'mass_distribution_reference', '-v7.3');
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

function mass_distribution = generate_random_distribution_constrained(cylinder_data, N)
    X = cylinder_data.X;
    Y = cylinder_data.Y;
    Z = cylinder_data.Z;
    np = cylinder_data.np;
    n = cylinder_data.n;
    
    theta = rad2deg(atan(Y./X));
    first_jump_idx = find(Y(1,:)==max(Y(1,:)));
    second_jump_idx = find(Y(1,:)==min(Y(1,:)));
    theta(:,first_jump_idx:second_jump_idx) = theta(:,first_jump_idx:second_jump_idx) + 180;
    theta(:,second_jump_idx:end) = theta(:,second_jump_idx:end) + 360;
    
    np_constrained = floor(np/3) + 1;
    height_division_size = floor(n/(N/2));
    
    dZ = 0.1;
    deltaZ_max = 6;
    deltaZ_random = 0.5;
    Nmax = deltaZ_max/dZ;
    Nrandom = deltaZ_random/dZ;
    Nsubtract = Nmax+Nrandom+1;

    plane_idx_anterior = randi([1, floor(np_constrained/2)], 1, N/2);
    plane_idx_lateral = randi([floor(np_constrained/2) + 1, np_constrained], 1, N/2);
    
    height_lower_idx = Nsubtract;
    height_upper_idx = n-Nsubtract;
    height_idx = nan(N/2+1, 1);
    height_idx(1) = height_lower_idx;
    for k=1:(N/2-1)
        height_idx(k+1) = k*height_division_size;
    end
    height_idx(end) = height_upper_idx;

    mass_distribution = struct;
    for k=1:N
        if k <= N/2
            current_plane_idx = plane_idx_anterior(k);
            current_height_idx = randi([height_idx(k), height_idx(k+1)], 1, 1);
        else
            k_shifted = k - N/2;
            current_plane_idx = plane_idx_lateral(k_shifted);
            current_height_idx = randi([height_idx(k_shifted), height_idx(k_shifted + 1)], 1, 1);
        end
        mass_distribution(k).mass = 1;
        mass_distribution(k).plane_idx = current_plane_idx;
        mass_distribution(k).height_idx = current_height_idx;
        mass_distribution(k).X = X(current_height_idx, current_plane_idx);
        mass_distribution(k).Y = Y(current_height_idx, current_plane_idx);
        mass_distribution(k).Z = Z(current_height_idx, current_plane_idx);
        mass_distribution(k).r = sqrt((mass_distribution(k).X)^2 + (mass_distribution(k).Y)^2);
        mass_distribution(k).theta = theta(current_height_idx, current_plane_idx);
        mass_distribution(k).vector = [mass_distribution(k).X; mass_distribution(k).Y; mass_distribution(k).Z];
    end

end

function intermarker_distance_bool = check_intermarker_distance(mass_distribution, N)
    intermarker_distance_bool = false;
    height_diff_threshold = 15;
    theta_diff_threshold = 60;
    height_list_anterior = nan(N/2,1);
    height_list_lateral = nan(N/2,1);
    theta_diff_list = nan(N/2,1);
    for k=1:N/2
        k_shifted = k + N/2;
        height_list_anterior(k) = mass_distribution(k).Z;
        height_list_lateral(k) = mass_distribution(k_shifted).Z;
        theta_diff_list(k) = mass_distribution(k_shifted).theta - mass_distribution(k).theta;
    end
    height_list_anterior_diff_min = min(diff(height_list_anterior));
    height_list_lateral_diff_min = min(diff(height_list_lateral));
    theta_diff_list_min = min(theta_diff_list);
    
    if height_list_anterior_diff_min > height_diff_threshold && height_list_lateral_diff_min > height_diff_threshold
        if theta_diff_list_min > theta_diff_threshold
            intermarker_distance_bool = true;
        end
    end
end