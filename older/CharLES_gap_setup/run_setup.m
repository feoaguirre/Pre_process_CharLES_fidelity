%% Setup
clear; clc; close all;

addpath('./helping_files');

% Getting input
input_file_name = "input_parameters.txt";
[parameters, input_file_data, input_file_title_row_indices] = ...
    setup_functions.read_parameters_from_txt(input_file_name);

folder_name = parameters.Run_name;
folder_path = "../" + folder_name;

if isfolder(folder_path)
    error("Simulation wasn't setup because the folder %s already exists", folder_name);
end

% Creating analysis folder
analysis_folder_path = folder_path + "_analysis/";
mkdir(analysis_folder_path);

analysis_files = ["read_probe_data_to_struct.m", "apply_fft.m", "analyze_results.m"];
for i = 1:length(analysis_files)
    copyfile("./helping_files/" + analysis_files(i), analysis_folder_path);
end

mkdir(analysis_folder_path + "analysis_images/");

% Creating system of folders for running the simulation
subfolders = ["/logs", "/helping_files/space_probe_coordinates", ...
              "/results/images", "/results/mesh_images"];

for i = 1:length(subfolders)
    mkdir(folder_path + subfolders(i));
end

% Calculating additional parameters
if parameters.HALF_W == 0
    parameters.HALF_W = 0.5 * parameters.HCP_DELTA_VAL / (2^parameters.MAX_REFINEMENT_LEVEL);
else
    parameters.HALF_W = parameters.HALF_W / 2; % B/c user inputs full value of W
end

if parameters.WALL_BC == "ISOTHERMAL"
    adiabatic_flag = 0;
    parameters.WALL_BC = parameters.WALL_BC + " T_WALL " + ...
                     string(parameters.T_WALL * parameters.T_INF);
else
    adiabatic_flag = 1;
end
if parameters.WALL_MODELING == "Yes"
    parameters.WALL_BC = "WM_ALG_" + parameters.WALL_BC;
else
    parameters.WALL_BC = "WALL_" + parameters.WALL_BC;
end

parameters.TOTAL_DOMAIN_LENGTH = parameters.L + parameters.L1 + ...
                                 parameters.L2 + parameters.L3;

parameters.TRANSIENT_SIMTIME_VAL = parameters.TRANSIENT_SIMTIME_VAL * ...
    parameters.TOTAL_DOMAIN_LENGTH / parameters.U_INF;

parameters.STEADY_SIMTIME_VAL = parameters.STEADY_SIMTIME_VAL * ...
    parameters.TOTAL_DOMAIN_LENGTH / parameters.U_INF + ...
    parameters.TRANSIENT_SIMTIME_VAL;

% Creating space probe files
for i = 1:length(parameters.space_probe_list)
    probe = parameters.space_probe_list{i};

    write_path = folder_path + "/helping_files/space_probe_coordinates/" + probe.file_name;

    [x1, x2, y1, y2, number_of_probes_along_x, number_of_probes_along_y] = setup_functions.get_space_probe_corners(probe, parameters);

    generate_space_probes(x1, x2, y1, y2, number_of_probes_along_x, number_of_probes_along_y, write_path);
end

% Creating text file with list of parameters
for i = 1:length(input_file_data)
    input_file_data{i} = {string(strjoin(input_file_data{i}, ", ")) + newline};
end

input_file_data = [ ...
    input_file_data(1:input_file_title_row_indices(2)-1); ...
    {{sprintf('Domain width (W), %f, double\n', 2*parameters.HALF_W)}}; ...
    input_file_data(input_file_title_row_indices(2):end) ...
];

fid = fopen("../" + parameters.Run_name + "_parameters.txt", 'w');
input_file_data = string(input_file_data); 
fprintf(fid, '%s', input_file_data{:});
fclose(fid);

% Creating input and run files
input_and_run_file_names = [ ...
    "surfer.in", "stitch.in", "transient_charles_ig.in", ...
    "steady_charles_ig.in", "run-surfer.sh", "run-stitch.sh", ...
    "run-transient-charles.sh", "run-steady-charles.sh", "move_logs.sh" ];

param_names = fieldnames(parameters);

for i = 1:length(input_and_run_file_names)

    file_name = input_and_run_file_names(i);
    template_file_path = "./helping_files/template-" + file_name;
    file_write_path = folder_path + "/" + file_name;
    fprintf("Creating file %s\n", file_write_path);

    file_data = readlines(template_file_path);

    % Handle .in files
    if endsWith(file_name, ".in")
        for j = 1:length(file_data)
            line = file_data(j);
            tokens = regexp(line, 'DEFINE (.+) = ([\.\d]+)', 'tokens');
            if ~isempty(tokens)
                param_name = string(tokens{1}{1});
                if any(strcmp(param_name, param_names))
                    value = string(parameters.(param_name));
                    file_data(j) = regexprep(line, tokens{1}{2}, value, 'once');
                end
            end
        end

    % Handle .sh files
    elseif endsWith(file_name, ".sh")
        for j = 1:length(file_data)
            line = file_data(j);
            tokens = regexp(line, '#PBS -N (.+)', 'tokens');
            if ~isempty(tokens)
                new_name = parameters.Run_name + "_" + tokens{1}{1};
                file_data(j) = regexprep(line, tokens{1}{1}, new_name, 'once');
                break;
            end
        end
    end

    % Special handling for steady_charles_ig.in and transient_charles_ig.in - Adding probes
    if contains(file_name, "charles_ig.in")
        if file_name == "transient_charles_ig.in"
            run_step = "transient";
            
        elseif file_name == "steady_charles_ig.in"
            run_step = "steady";

            % Space probes
            for k = 1:length(parameters.space_probe_list)
                probe = parameters.space_probe_list{k};
                probe_name = erase(probe.file_name, "." + extractAfter(probe.file_name, "."));
                new_line = sprintf(['POINTCLOUD_PROBE NAME=./results/space_probes/%s/measurement ' ...
                    'INTERVAL $(SPACE_PROBES_WRITE_INTERVAL) ' ...
                    'GEOM=FILE ./helping_files/space_probe_coordinates/%s ' ...
                    'FORMAT=BINARY VARS=comp(avg(u),0) comp(avg(u),1) avg(p) avg(rho) avg(T) ' ...
                    'comp(curl(avg(u)),2) avg(div(u))'], ...
                    probe_name, probe.file_name);
                file_data(end+1) = string(new_line);
            end
        end

        % Time probes
        for k = 1:length(parameters.time_probe_list)
            probe = parameters.time_probe_list{k};
            [x1, x2, y1, y2, num_of_probes] = setup_functions.get_time_probe_corners(probe, parameters);
            new_line = sprintf(['PROBE NAME=./results/time_probes/%s/%s/measurement ' ...
                'INTERVAL $(TIME_PROBES_WRITE_INTERVAL) ' ...
                'GEOM=LINE %0.8f %0.8f 0 %0.8f %0.8f 0 %d ' ...
                'VARS=comp(u,0) comp(u,1) p rho T comp(vorticity(),2) div(u)'], ...
                probe.name, run_step, x1, y1, x2, y2, num_of_probes);
            file_data(end+1) = string(new_line);
        end
    end

    % Write file
    writelines(file_data, file_write_path);
end

disp("Finished Simulation Setup");