classdef setup_functions

   methods(Static)

        % Structure based equivalents of the classes I built in python version
        function obj = SpaceProbeObj(file_name, num_of_parallel_probes, num_of_perpendicular_probes, distance_from_wall)
            obj = struct( ...
                'file_name', string(file_name), ...
                'num_of_parallel_probes', num_of_parallel_probes, ...
                'num_of_perpendicular_probes', num_of_perpendicular_probes, ...
                'distance_from_wall', distance_from_wall ...
            );
        end
        
        function obj = TimeProbeObj(name, num_of_probes)
            obj = struct( ...
                'name', string(name), ...
                'num_of_probes', num_of_probes ...
            );
        end
        
        % Parameter name translation
        function out = translate_parameter_name(name)
            switch string(name)
                case "Run name"
                    out = "Run_name";
                case "Pregap slip length (L1/delta_star)"
                    out = "L1";
                case "Pregap noslip length (L2/delta_star)"
                    out = "L2";
                case "Gap length (L/delta_star)"
                    out = "L";
                case "Gap depth (D/delta_star)"
                    out = "D";
                case "Gap width (W/delta_star - 0 for 2D)"
                    out = "HALF_W";
                case "Postgap length (L3/delta_star)"
                    out = "L3";
                case "Domain height (H/delta_star)"
                    out = "H";
                case "Mesh size at far field in units of delta_star (HCP_DELTA)"
                    out = "HCP_DELTA_VAL";
                case "Max refinement level (hat(gamma))"
                    out = "MAX_REFINEMENT_LEVEL";
                case "Height above solids with same refinement (D)"
                    out = "REFINEMENT_HEIGHT_ABOVE_PLATE";
                case "Number of cells in every following layer of refinement (NLAYERS)"
                    out = "NLAYERS_VAL";
                case "Number of smoothing iterations (NSMOOTH)"
                    out = "NSMOOTH_VAL";
                case "Reynolds delta_star"
                    out = "RE_DELTA_STAR";
                case "Mach"
                    out = "MACH";
                case "Inflow speed (U_inf)"
                    out = "U_INF";
                case "Inflow density (rho_inf)"
                    out = "RHO_INF";
                case "Inflow temperature (T_inf)"
                    out = "T_INF";
                case "gamma"
                    out = "G";
                case "Prandtl number"
                    out = "PR_LAM_VAL";
                case "Mu power law exponent"
                    out = "MU_POWER_LAW_VAL";
                case "Wall boundary condition (ADIABATIC/ISOTHERMAL)"
                    out = "WALL_BC";
                case "Wall temperature (T_wall/T_inf - ignored in adiabatic case)"
                    out = "T_WALL";
                case "Simulation time to run transient part (in Flow Through Time)"
                    out = "TRANSIENT_SIMTIME_VAL";
                case "Simulation time to run steady part (in Flow Through Time)"
                    out = "STEADY_SIMTIME_VAL";
                case "CFL"
                    out = "CFL_VAL";
                case "Check interval"
                    out = "CHECK_INTERVAL_VAL";
                case "Image write interval"
                    out = "IMAGE_WRITE_INTERVAL";
                case "Space probes write interval"
                    out = "SPACE_PROBES_WRITE_INTERVAL";
                case "Time probes write interval"
                    out = "TIME_PROBES_WRITE_INTERVAL";
                case "Sub Grid Scale Model (VREMAN/NONE)"
                    out = "SGS_MODEL_VAL";
                case "Use wall modeling? (Yes/No)"
                    out = "WALL_MODELING";
                otherwise
                    error("Unknown parameter name: %s", name);
            end
        end
        
        % Read input file
        function [params, file_data, title_row_indices] = read_parameters_from_txt(filename)
        
            fid = fopen(filename, 'r');
            raw_lines = textscan(fid, '%s', 'Delimiter', '\n');
            fclose(fid);
        
            raw_lines = raw_lines{1};
        
            % Split lines like Python
            file_data = cell(length(raw_lines),1);
            for i = 1:length(raw_lines)
                file_data{i} = strsplit(strtrim(raw_lines{i}), ", ");
            end
        
            params = struct();
            params.space_probe_list = {};
            params.time_probe_list = {};
        
            title_row_indices = [];
        
            for idx = 1:length(file_data)
                row = file_data{idx};
                if length(row) == 1 % Reading title rows
                    title_row_indices(end+1) = idx;
                    continue
        
                elseif length(row) == 2 % Reading space probe data
                    file_name = replace(row{1}, " ", "_") + ".csv";
                    vals = strsplit(row{2}, "-");
                    obj = setup_functions.SpaceProbeObj(file_name, ...
                        str2double(vals{1}), ...
                        str2double(vals{2}), ...
                        str2double(vals{3}));
                    params.space_probe_list{end+1} = obj;
                    continue
        
                elseif length(row) == 4 % Reading time probe data
                    name = replace(row{1}, " ", "_");
                    obj = setup_functions.TimeProbeObj(name, str2double(row{2}));
                    params.time_probe_list{end+1} = obj;
                    continue
                end
        
                % Reading other parameters
                key = row{1};
                value = row{2};
                number_type = row{3};
        
                switch number_type
                    case "double"
                        value = str2double(value);
                    case "int"
                        value = str2double(value);
                    case "array"
                        value = replace(value, "-", ",");
                end
        
                field_name = setup_functions.translate_parameter_name(key);
                params.(field_name) = value;
            end
        end
        
        
        % Space probe geometry
        function [x1, x2, y1, y2, nx, ny] = get_space_probe_corners(probe, parameters)
        
            eps_val = 1e-8;
        
            switch probe.file_name
        
                case "Pregap_noslip.csv"
                    x1 = -parameters.L2;
                    x2 = 0;
                    y1 = eps_val;
                    y2 = eps_val + probe.distance_from_wall;
                    nx = probe.num_of_parallel_probes;
                    ny = probe.num_of_perpendicular_probes;
        
                case "Gap_front_wall.csv"
                    x1 = eps_val;
                    x2 = eps_val + probe.distance_from_wall;
                    y1 = -parameters.D + eps_val;
                    y2 = 0;
                    nx = probe.num_of_perpendicular_probes;
                    ny = probe.num_of_parallel_probes;
        
                case "Gap.csv"
                    x1 = eps_val;
                    x2 = parameters.L - eps_val;
                    y1 = -parameters.D + eps_val;
                    y2 = -parameters.D + eps_val + probe.distance_from_wall;
                    nx = probe.num_of_parallel_probes;
                    ny = probe.num_of_perpendicular_probes;
        
                case "Gap_back_wall.csv"
                    x1 = parameters.L - eps_val - probe.distance_from_wall;
                    x2 = parameters.L - eps_val;
                    y1 = -parameters.D + eps_val;
                    y2 = 0;
                    nx = probe.num_of_perpendicular_probes;
                    ny = probe.num_of_parallel_probes;
        
                case "Postgap.csv"
                    x1 = parameters.L;
                    x2 = parameters.L + parameters.L3;
                    y1 = eps_val;
                    y2 = eps_val + probe.distance_from_wall;
                    nx = probe.num_of_parallel_probes;
                    ny = probe.num_of_perpendicular_probes;
        
                otherwise
                    error("Unknown space probe: %s", probe.file_name);
            end
        end
        
        
        % Time probe geometry
        function [x1, x2, y1, y2, n] = get_time_probe_corners(probe, parameters)
        
            eps_val = 1e-8;
        
            switch probe.name
        
                case "Shear_layer"
                    x1 = 0;
                    x2 = parameters.L;
                    y1 = eps_val;
                    y2 = eps_val;
                    n = probe.num_of_probes;
        
                case "Gap_front_wall"
                    x1 = eps_val;
                    x2 = eps_val;
                    y1 = -parameters.D + eps_val;
                    y2 = 0;
                    n = probe.num_of_probes;
        
                case "Gap_floor"
                    x1 = eps_val;
                    x2 = parameters.L - eps_val;
                    y1 = -parameters.D + eps_val;
                    y2 = -parameters.D + eps_val;
                    n = probe.num_of_probes;
        
                case "Gap_back_wall"
                    x1 = parameters.L - eps_val;
                    x2 = parameters.L - eps_val;
                    y1 = -parameters.D + eps_val;
                    y2 = 0;
                    n = probe.num_of_probes;
        
                case "Gap_mid_height"
                    x1 = eps_val;
                    x2 = parameters.L - eps_val;
                    y1 = -0.5 * parameters.D;
                    y2 = -0.5 * parameters.D;
                    n = probe.num_of_probes;
        
                otherwise
                    error("Unknown time probe: %s", probe.name);
            end
        end
   end
end