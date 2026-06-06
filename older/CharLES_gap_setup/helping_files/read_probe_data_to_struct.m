function [y_plus_data_struct, time_probe_data_struct, space_probe_data_struct] = read_probe_data_to_struct(main_dir, num_of_perpendicular_probes, parameters_being_measured_by_space_probes, parameters_being_measured_by_time_probes)
    % Reading space probes
    space_probes_dir = [main_dir, 'space_probes/'];
    probe_names = get_probe_names(space_probes_dir);
    
    for probe_idx = 1:length(probe_names)
        probe_name = probe_names(probe_idx);
        probe_name = probe_name{1};
        [space_probe_data_struct.(probe_name), cascade_parameter_ind] = read_space_probe_xyz([space_probes_dir, probe_name, '/measurement.pbin'], num_of_perpendicular_probes(probe_idx));

       	file_names = get_space_file_names([space_probes_dir, probe_name]);
        for file_name_idx = 1:length(file_names)
            file_name = file_names(file_name_idx);
            file_name = file_name{1};
            file_location = [space_probes_dir, probe_name, '/', file_name];
            space_probe_data_struct.(probe_name).DATA(file_name_idx) = read_space_probe_data(file_location, num_of_perpendicular_probes(probe_idx), parameters_being_measured_by_space_probes, cascade_parameter_ind);
        end
    end

    % Reading time probes
    time_probes_dir = [main_dir, 'time_probes/'];
    probe_names = get_probe_names(time_probes_dir);
    
    for probe_idx = 1:length(probe_names)
        probe_name = probe_names(probe_idx);
        probe_name = probe_name{1};

        time_probe_data_struct.(probe_name) = read_time_probe_xyz([time_probes_dir, probe_name, '/transient', '/measurement.README']);

        file_names = get_time_file_names([time_probes_dir, probe_name, '/transient']);
        for file_name_idx = 1:length(file_names)
            file_name = file_names(file_name_idx);
            file_name = file_name{1};
            time_probe_data_struct.(probe_name).(parameters_being_measured_by_time_probes(file_name_idx)) = read_time_probe_data(time_probes_dir, probe_name, file_name);
            steady_idx = time_probe_data_struct.(probe_name).(parameters_being_measured_by_time_probes(file_name_idx)).steady_idx;
            time_probe_data_struct.(probe_name).(parameters_being_measured_by_time_probes(file_name_idx)).perturbations = time_probe_data_struct.(probe_name).(parameters_being_measured_by_time_probes(file_name_idx)).DATA - mean(time_probe_data_struct.(probe_name).(parameters_being_measured_by_time_probes(file_name_idx)).DATA(steady_idx:end,:), 1);
        end
    end

    % Reading y_plus data
    y_plus_file_names = dir([main_dir, 'y_plus/']);
    last_y_plus_filepath = [main_dir, 'y_plus/', y_plus_file_names(end).name];
    y_plus_data = readtable(last_y_plus_filepath, NumHeaderLines=6, FileType="text");
    y_plus_data_struct.y_plus.x = y_plus_data.Var3;
    y_plus_data_struct.y_plus.DATA = y_plus_data.Var5;
end

% Helper functions
function time_probes_xyz_struct = read_time_probe_xyz(readme_file_path)
    file_info = readtable(readme_file_path, 'NumHeaderLines', 2, FileType='text');
    axes = ['x', 'y', 'z'];
    for axes_idx = 1:length(axes)
        time_probes_xyz_struct.(axes(axes_idx)) = file_info.(axes_idx + 1);
    end
end

function time_probes_data_struct = read_time_probe_data(time_probes_dir, probe_name, file_name)
    % Transient data
    file_location = [time_probes_dir, probe_name, '/transient/', file_name];
    file_info = readmatrix(file_location, 'NumHeaderLines', 3, FileType="text");
    time_probes_data_struct.step = file_info(:,1);
    time_probes_data_struct.time = file_info(:,2);
    time_probes_data_struct.DATA = file_info(:,4:end); % Each column represents a different probe
    time_probes_data_struct.steady_idx = length(time_probes_data_struct.step) + 1;

    % Steady data
    file_location = [time_probes_dir, probe_name, '/steady/', file_name];
    file_info = readmatrix(file_location, 'NumHeaderLines', 3, FileType="text");
    time_probes_data_struct.step = [time_probes_data_struct.step; file_info(:,1)];
    time_probes_data_struct.time = [time_probes_data_struct.time; file_info(:,2)];
    time_probes_data_struct.DATA = [time_probes_data_struct.DATA; file_info(:,4:end)]; % Each column represents a different probe
end

function time_file_names = get_time_file_names(file_dir)
    file_info = dir(file_dir);
    time_file_names = {};
    for i = 3:length(file_info)
        file_name = file_info(i).name;
        if length(file_name) > 2 && file_name(end-5:end) ~= "README"
            time_file_names{end+1} = file_name; 
        end
    end
end

function probe_names = get_probe_names(main_dir)
    probe_results_contents = dir(main_dir);
    probe_names = cell(length(probe_results_contents)-2, 1);
    for i = 3:length(probe_results_contents)
        probe_names{i-2} = probe_results_contents(i).name;
    end
end

function space_file_names = get_space_file_names(file_dir)
    file_info = dir(file_dir);
    space_file_names = cell(length(file_info)-4, 1);
    for i = 3:length(file_info)
        file_name = file_info(i).name;
        if length(file_name) > 2 && file_name(end-3:end) == ".pcd"
            space_file_names{i-2} = file_name; 
        end
    end
end

function [probe_xyz_struct, ind] = read_space_probe_xyz(pbin, num_of_perpendicular_probes)
    fh = fopen(pbin,'r');
    endian = 'l';
    magic_number = fread(fh,1,'int64',endian); % magic number
    if magic_number ~= 1235813 
        endian = 'b';
    end
    fread(fh,1,'int64',endian); % skip version
    np = fread(fh,1,'int64',endian); % number of points
    assert(2 == fread(fh,1,'int64',endian)); % 0:no data, 1:delta, 2:index
    buf = fread(fh,3*np,'double',endian); % always double for consistency w/ other pbins
    ind = fread(fh,np,'int64',endian); % global index (should equal ordering in ascii file)
    fclose(fh);
    
    probe_xyz = zeros(np,3);
    for i=1:np
	for j=1:3
            probe_xyz(1+ind(i),j) = buf(3*(i-1)+j);
        end
    end

    axes = ['x', 'y', 'z'];
    for axes_idx = 1:length(axes)
        probe_xyz_struct.(axes(axes_idx)) = reshape(probe_xyz(:,axes_idx), num_of_perpendicular_probes, []);
    end
end

function probe_data_struct = read_space_probe_data(pcd, num_of_perpendicular_probes, variable_names, ind)
    fh = fopen(pcd,'r');
    endian = 'l';
    magic_number = fread(fh,1,'int64',endian); % magic number
    if magic_number ~= 1235813
        endian = 'b';
    end
    fread(fh,1,'int64',endian); % skip version
    np = fread(fh,1,'int64',endian); % number of points
    nv = fread(fh,1,'int64',endian); % number of variables
    prec = fread(fh,1,'int64',endian); % 0 float/ 1 double
    
    if (prec == 1)
        buf = fread(fh,nv*np,'double',endian);
    else
        buf = fread(fh,nv*np,'float',endian);
    end
    fclose(fh);

    probe_data = zeros(nv,np);
    for j = 1:nv
        for i = 1:np
            probe_data(j,1+ind(i)) = buf(i+(j-1)*np);
        end
    end

    probe_data = probe_data';
    
    for variable_name_idx = 1:length(variable_names)
        probe_data_struct.(variable_names(variable_name_idx)) = reshape(probe_data(:,variable_name_idx), num_of_perpendicular_probes, []);
    end
end

