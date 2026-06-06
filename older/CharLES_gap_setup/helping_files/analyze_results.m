%{
1. num_of_perpendicular_probes can be taken from the input_parameters file 
and needs to be in Alphabetical order: 
Back wall, gap, front wall, postgap, pregap, shear layer

2. Check to make sure in which order parameters_being_measured appear.

3. All the data is read from the files and organized into the probe_data
structures. Each parameter is organizes such that each column of the matrix
represents a different value of the axis perpendicular to the wall the probes
are placed on. For example, for Pregap_noslip these are x values, and
for Gap_back_wall these are y values. All Data follows the same order,
meaning that for a given set of probes, rho is organized just like x, y, and z are.
%}

clear; close all; clc;

[~, run_name, extension] = fileparts(pwd);
run_name = append(run_name, extension);
clear extension;
run_name = string(run_name(1:end-9));
main_dir = 'results/';
num_of_perpendicular_probes = [201, 101, 101];
parameters_being_measured_by_space_probes = ["p", "rho", "T", "u", "v", "vorticity"]; % Organize according to the README generated for these probes
parameters_being_measured_by_time_probes = ["u", "v", "vorticity", "divu", "p", "rho", "T"]; % Organize alphabetically (by file order in the folder)

[y_plus_data, time_probe_data, space_probe_data] = read_probe_data_to_struct(main_dir, num_of_perpendicular_probes, parameters_being_measured_by_space_probes, parameters_being_measured_by_time_probes);
Gap_length = 25;

%% Checking integral values (delta_star, theta, H) along the pregap plate
pregap_integral_values = delta_star_analysis('Pregap_noslip', space_probe_data, 1);

function integral_values = delta_star_analysis(probes_being_analyzed, space_probe_data, BL_end_vorticity_percent)

    x = space_probe_data.(probes_being_analyzed).x(1, :);
    
    % Simulation values
    time_being_analyzed_idx = length(space_probe_data.(probes_being_analyzed).DATA);
    variable_data = space_probe_data.(probes_being_analyzed).DATA(time_being_analyzed_idx);
    integral_values.simulation.delta_star = [];
    integral_values.simulation.theta = [];
    integral_values.simulation.H = [];
    for i = 3:1:length(x)
        y = space_probe_data.(probes_being_analyzed).y(:, i);
        vorticity = variable_data.vorticity(:, i);
        u = variable_data.u(:, i);
        rho = variable_data.rho(:, i);
        max_vorticity = vorticity(1);
        vorticity_tag = vorticity - max_vorticity * 0.01 * BL_end_vorticity_percent;
        BL_end_idx = find(vorticity_tag(1:end-1) .* vorticity_tag(2:end) <= 0) + 1;
        rho_u_e = rho(BL_end_idx) * u(BL_end_idx);
        integrand = 1 - rho.*u ./ rho_u_e;
        integral_values.simulation.delta_star(end + 1) = trapz(y(1:BL_end_idx), integrand(1:BL_end_idx));
        integrand = rho.*u ./ rho_u_e .* (1 - u / u(BL_end_idx));
        integral_values.simulation.theta(end + 1) = trapz(y(1:BL_end_idx), integrand(1:BL_end_idx));
        integral_values.simulation.H = integral_values.simulation.delta_star ./ integral_values.simulation.theta;
    end
    
    % Plotting delta_star
    figure;
    hold on;
    plot(x(3:end), integral_values.simulation.delta_star, 'b');
    legend(location="northwest");
    xlabel("X");
    ylabel("\delta^*");
    grid;
    % exportgraphics(gca, "./analysis_images/delta_star.eps");
    
    % Plotting theta
    figure;
    hold on;
    plot(x(3:end), integral_values.simulation.theta, 'b');
    legend(location="northwest");
    xlabel("X");
    ylabel("\theta");
    grid;
    % exportgraphics(gca, "./analysis_images/theta.eps");

    % Plotting H
    figure;
    hold on;
    plot(x(3:end), integral_values.simulation.H, 'b', DisplayName='Simulation');
    legend(location="south");
    xlabel("X");
    ylabel("H");
    ylim([-1 5]);
    grid;
    % exportgraphics(gca, "./analysis_images/H.eps");

    % Plotting additional profiles
    % for probe_idx = 10:10:80
    %     figure;
    %     plot(variable_data.u(:, probe_idx), y, 'b');
    %     xlabel("u");
    %     ylabel("y");
    %     grid;
    % end
end

%% Plotting y+ along all of SOLIDS
figure;
hold on;
plot(y_plus_data.y_plus.x, y_plus_data.y_plus.DATA);
xlabel('x');
ylabel('y+');

%% Calculating and comparing p' FFT results for all time probes along the shear layer and the other probe strips
all_shear_layer_fft_results = compare_FFT_along_probe_strip(time_probe_data, "Shear_layer", Gap_length, "x");
% all_mid_height_fft_results = compare_FFT_along_probe_strip(time_probe_data, "Gap_mid_height", Gap_length, "x");
% all_front_wall_fft_results = compare_FFT_along_probe_strip(time_probe_data, "Gap_front_wall", Gap_length, "y");
% all_back_wall_fft_results = compare_FFT_along_probe_strip(time_probe_data, "Gap_back_wall", Gap_length, "y");

function all_probe_strip_fft_results = compare_FFT_along_probe_strip(time_probe_data, probe_strip_name, Gap_length, x_or_y)
        % Calulating FFT for all probes
    number_of_probes = size(time_probe_data.(probe_strip_name).p.DATA);
    number_of_probes = number_of_probes(2);
    
    steady_idx = time_probe_data.(probe_strip_name).p.steady_idx;

    time = time_probe_data.(probe_strip_name).p.time(steady_idx:end);
    if mod(length(time), 2) == 1
        time = time(2:end);
    end
    
    for time_probe_being_analyzed = 1:number_of_probes
        signal = time_probe_data.(probe_strip_name).p.perturbations(steady_idx:end,time_probe_being_analyzed);
        if mod(length(signal), 2) == 1
            signal = signal(2:end);
        end
        fft_results.("probe" + time_probe_being_analyzed) = apply_fft(time, signal, 10^-5);
        % Converting from St_deltaStar to St_L
        fft_results.("probe" + time_probe_being_analyzed).frequencies = fft_results.("probe" + time_probe_being_analyzed).frequencies * Gap_length;
        fft_results.("probe" + time_probe_being_analyzed).peak_frequencies = fft_results.("probe" + time_probe_being_analyzed).peak_frequencies * Gap_length;
    end
    
        % Preparing data for a contour plot
    array_lengths = zeros(1, number_of_probes);
    for time_probe_being_analyzed = 1:number_of_probes
        array_lengths(time_probe_being_analyzed) = length(fft_results.("probe" + time_probe_being_analyzed).frequencies);
    end
    max_length = max(array_lengths);
    
    x_matrix = (time_probe_data.(probe_strip_name).x/Gap_length)' + zeros(max_length, number_of_probes);
    y_matrix = (time_probe_data.(probe_strip_name).y)' + zeros(max_length, number_of_probes);
    if x_or_y == "x"
        matrix_to_plot = x_matrix;
        xlabel_txt = "x/L";
    elseif x_or_y == "y"
        matrix_to_plot = y_matrix;
        xlabel_txt = "y/\delta^*";
    end

    frequencies_matrix = zeros(max_length, number_of_probes);
    amplitudes_matrix = zeros(max_length, number_of_probes);
    
    for time_probe_being_analyzed = 1:number_of_probes
        frequencies = fft_results.("probe" + time_probe_being_analyzed).frequencies;
        log_amplitudes = log(fft_results.("probe" + time_probe_being_analyzed).amplitudes);
        log_amplitudes(isinf(log_amplitudes)) = -20; % So that I don't get -inf for amplitude 0
        frequencies_matrix(1:length(frequencies), time_probe_being_analyzed) = frequencies;
        amplitudes_matrix(1:length(log_amplitudes), time_probe_being_analyzed) = log_amplitudes;
    end

    all_probe_strip_fft_results = fft_results;
    all_probe_strip_fft_results.number_of_probes = number_of_probes;

        % Plotting
    figure;
    contourf(matrix_to_plot, frequencies_matrix, amplitudes_matrix, 100,'Linestyle','none');
    colorbar_obj = colorbar();
    colorbar_obj.Label.String = "log(|FFT(p')|)";
    xlabel(xlabel_txt);
    ylabel("St_L");
    ylim([0, 2]);
    clim([-10,-5.5]);
end

%% Plotting p' signal and fft of a single time probe
time_analysis(time_probe_data, all_shear_layer_fft_results, 7, 5);

function time_analysis(time_probe_data, all_shear_layer_fft_results, time_probe_being_analyzed, num_of_harmonic_lines)
    time = time_probe_data.Shear_layer.p.time;
    signal = time_probe_data.Shear_layer.p.perturbations(:, time_probe_being_analyzed);

    % Plotting change in time
    figure;
    hold on;
    plot(time, signal);
    xlabel('Time since beginning of simulation');
    ylabel("p'=p-|p|");
    
    % Plotting FFT of results
    fft_results = all_shear_layer_fft_results.("probe" + time_probe_being_analyzed);
    
    main_frequency = all_shear_layer_fft_results.probe1.peak_frequencies(all_shear_layer_fft_results.probe1.peak_amplitudes == max(all_shear_layer_fft_results.probe1.peak_amplitudes));

    fft_results.harmonic_frequencies = zeros(1, num_of_harmonic_lines);
    for harmonic_idx = 1:num_of_harmonic_lines
        fft_results.harmonic_frequencies(harmonic_idx) = main_frequency * harmonic_idx;
    end
    
    max_amplitude = max(fft_results.amplitudes);
    min_amplitude = min(fft_results.amplitudes(fft_results.amplitudes > 0));

    figure;
    hold on;
    plot(fft_results.frequencies, fft_results.amplitudes);
    for harmonic_idx = 1:num_of_harmonic_lines
        freq = fft_results.harmonic_frequencies(harmonic_idx);
        plot([freq freq], [min_amplitude max_amplitude], '--k');
    end
    plot(fft_results.peak_frequencies, fft_results.peak_amplitudes, 'or');
    xlabel("St_L");
    ylabel("|FFT(p')|");
    xlim([0 2]);
    set(gca, 'yscale', 'log');
end

%% Plotting curl(u) and div(u) along shear layer to see wave propagation, with option to plot other variables
plot_wave_propegation(time_probe_data, all_shear_layer_fft_results, "Shear_layer", Gap_length, 4, "x");
% plot_wave_propegation(time_probe_data, all_mid_height_fft_results, "Gap_mid_height", Gap_length, 4, "x");
% plot_wave_propegation(time_probe_data, all_front_wall_fft_results, "Gap_front_wall", Gap_length, 4, "y");

function plot_wave_propegation(time_probe_data, all_probe_strip_fft_results, probe_strip_name, Gap_length, periods_to_display, x_or_y)
    number_of_probes = all_probe_strip_fft_results.number_of_probes;
    time_length = length(time_probe_data.(probe_strip_name).p.time);
    time_matrix = time_probe_data.(probe_strip_name).p.time + zeros(time_length, number_of_probes);
    
    % Making plot show last periods of the disturbances in the measured signal
    main_frequency = all_probe_strip_fft_results.probe1.peak_frequencies(all_probe_strip_fft_results.probe1.peak_amplitudes == max(all_probe_strip_fft_results.probe1.peak_amplitudes));
    period_time = 1 / (main_frequency / Gap_length);
    t0 = time_matrix(end, 1) - periods_to_display * period_time;
    temp_time_array = time_matrix(:, 1) - t0;
    t0_idx = find(temp_time_array(1:end-1) .* temp_time_array(2:end) <= 0);
    t0_idx = t0_idx(1);
    L_time_matrix = (time_matrix(t0_idx:end, :) - t0) / Gap_length;
    
    x_matrix = (time_probe_data.(probe_strip_name).x/Gap_length)' + zeros(time_length, number_of_probes);
    x_matrix = x_matrix(t0_idx:end, :);
    y_matrix = (time_probe_data.(probe_strip_name).y)' + zeros(time_length, number_of_probes);
    y_matrix = y_matrix(t0_idx:end, :);
    if x_or_y == "x"
        matrix_to_plot = x_matrix;
        xlabel_txt = "x/L";
        drawlines_bool = 1;
    elseif x_or_y == "y"
        matrix_to_plot = y_matrix;
        xlabel_txt = "y/\delta^*";
        drawlines_bool = 0;
    end
    
    % Variables to be plotted
    pressure_matrix = time_probe_data.(probe_strip_name).p.perturbations(t0_idx:end, :);
    speed_matrix = time_probe_data.(probe_strip_name).v.perturbations(t0_idx:end, :);
    density_matrix = time_probe_data.(probe_strip_name).rho.perturbations(t0_idx:end, :);
    Temperature_matrix = time_probe_data.(probe_strip_name).T.DATA(t0_idx:end, :);
    vorticity_matrix = time_probe_data.(probe_strip_name).vorticity.perturbations(t0_idx:end, :);
    dilatation_matrix = time_probe_data.(probe_strip_name).divu.perturbations(t0_idx:end, :);
    
    % Plotting vorticity and dilatation
    if x_or_y == "x"
        second_x_matrix = x_matrix + 0.98;

        fig = figure;
        ax1 = axes('Parent', fig);
        hold on;
        contourf(ax1, x_matrix, L_time_matrix, vorticity_matrix, 100,'Linestyle','none');
        colorbar_obj = colorbar(ax1, "westoutside");
        colorbar_obj.Label.String = "vorticity";
        clim(ax1, [-0.8 0.8]);
        xlim(ax1, [0 2]);
        xticks([0 0.25 0.5 0.75 1 1.25 1.5 1.75 2]);
        xticklabels({'0', '0.25', '0.5', '0.75', '1', '0.75', '0.5', '0.25', '0'});
        xline(1,'k--');   % center line

        xlabel(xlabel_txt);
        ylabel("$t_L=(t-t_0) / \frac{L}{\delta^*}$", Interpreter="latex");

        vortex_firsthalf = drawline(0, 0.47, 3, 5, 0, 0.47);
        vortex_secondhalf_upper = drawline(0.47, 0.94, 5.25, 6.2, 0.47, 0.98);
        vortex_secondhalf_lower = drawline(0.47, 0.94, 4.75, 4.37, 0.47, 0.98);


        ax2 = axes('Parent', fig);
        hold on;
        contourf(ax2, second_x_matrix, L_time_matrix, flip(dilatation_matrix, 2), 100,'Linestyle','none');
        colorbar_obj2 = colorbar(ax2);
        colorbar_obj2.Label.String = "div(u)";
        clim(ax2, [-0.2 0.2]);
        ax2.Color = 'none';
        ax2.Visible = 'off';
        ax2.XTick = [];
        ax2.YTick = [];
        xlim(ax2, [0 2]);
        ax2.Position = ax1.Position;
        
        acousticwave_backward = drawline(1.1, 1.9, 6.7, 8.5, 1, 2);
        acousticwave_forward = drawline(1.1, 2, 8.2, 5.6, 1, 2);
    end

    function p = drawline(x0, x1, y0, y1, xstart, xend)
        p = polyfit([x0 x1], [y0 y1], 1);
        x = [xstart xend];
        plot(x, polyval(p, x),'--r');
    end

    
    %     % Plotting v'
    % figure;
    % contourf(matrix_to_plot, L_time_matrix, speed_matrix, 100,'Linestyle','none');
    % colorbar_obj = colorbar();
    % colorbar_obj.Label.String = "v'/U_\infty";
    % xlabel(xlabel_txt);
    % ylabel("$t_L=(t-t_0) / \frac{L}{\delta^*}$", Interpreter="latex");

    % % Drawing wave_lines
    % hold on;
    % forward_line.x = [0.05, 0.47];
    % forward_line.y = [0.47, 1];
    % forward_line.derivative = (forward_line.y(2) - forward_line.y(1)) / (forward_line.x(2) - forward_line.x(1));
    % forward_line.y0 = forward_line.y(1) - forward_line.derivative * forward_line.x(1);
    % 
    % backward_line.x = [0.57, 0.95];
    % backward_line.y = [1, 0.77];
    % backward_line.derivative = (backward_line.y(2) - backward_line.y(1)) / (backward_line.x(2) - backward_line.x(1));
    % backward_line.y0 = backward_line.y(1) - backward_line.derivative * backward_line.x(1);
    % 
    % drawlines(gca, forward_line, backward_line, periods_to_display, drawlines_bool);
    % 
    % function drawlines(ax, forward_line, backward_line, periods_to_display, drawlines_bool)
    %     if drawlines_bool % For probes placed along L
    %         for line_info = {forward_line, backward_line}
    %             line_info = line_info{1};
    %             for y0 = line_info.y0:1:(line_info.y0 + periods_to_display - 1)
    %                 y1 = y0 + line_info.derivative;
    %                 plot(ax, [0 1], [y0 y1], '--r');
    %             end
    %         end
    %         ylim([0 periods_to_display]);
    %     else % For probes placed along D
    %         line_info = forward_line;
    %         for y0 = line_info.y0:1:(line_info.y0 + periods_to_display - 1)
    %             plot(ax, [0 -10], [y0 y0], '--r');
    %         end
    %     end
    % end
    % 
    %     % Plotting p'
    % figure;
    % hold on;
    % contourf(matrix_to_plot, L_time_matrix, pressure_matrix, 100,'Linestyle','none');
    % colorbar_obj = colorbar();
    % colorbar_obj.Label.String = "p'/p_\infty";
    % xlabel(xlabel_txt);
    % ylabel("$t_L=(t-t_0) / \frac{L}{\delta^*}$", Interpreter="latex");
    % 
    % drawlines(gca, forward_line, backward_line, periods_to_display, drawlines_bool);
    % 
    % 
    %     % Plotting u
    % figure;
    % hold on;
    % contourf(matrix_to_plot, L_time_matrix, time_probe_data.(probe_strip_name).u.DATA(t0_idx:end, :), 100,'Linestyle','none');
    % colorbar_obj = colorbar();
    % colorbar_obj.Label.String = "u/U_\infty";
    % xlabel(xlabel_txt);
    % ylabel("$t_L=(t-t_0) / \frac{L}{\delta^*}$", Interpreter="latex");
    % 
    %     % Plotting rho'
    % figure;
    % hold on;
    % contourf(matrix_to_plot, L_time_matrix, density_matrix, 100,'Linestyle','none');
    % colorbar_obj = colorbar();
    % colorbar_obj.Label.String = "\rho'/\rho_\infty";
    % xlabel(xlabel_txt);
    % ylabel("$t_L=(t-t_0) / \frac{L}{\delta^*}$", Interpreter="latex");
    % 
    % drawlines(gca, forward_line, backward_line, periods_to_display, drawlines_bool);

    %     % Plotting T
    % figure;
    % hold on;
    % contourf(matrix_to_plot, L_time_matrix, Temperature_matrix, 100,'Linestyle','none');
    % colorbar_obj = colorbar();
    % colorbar_obj.Label.String = "T/T_\infty";
    % xlabel(xlabel_txt);
    % ylabel("$t_L=(t-t_0) / \frac{L}{\delta^*}$", Interpreter="latex");
    %
    % drawlines(gca, forward_line, backward_line, periods_to_display, drawlines_bool);
end

%% Creating contour plot of speed inside gap
plot_speed_contour(space_probe_data);

function plot_speed_contour(space_probe_data)
    probe_being_analyzed = "Gap";
    x = space_probe_data.(probe_being_analyzed).x;
    y = space_probe_data.(probe_being_analyzed).y;
    u = space_probe_data.(probe_being_analyzed).DATA(end).u;
    v = space_probe_data.(probe_being_analyzed).DATA(end).v;
    rho = space_probe_data.(probe_being_analyzed).DATA(end).rho;
    
    streamline_seeds.x = linspace(min(x(:)), max(x(:)), 7);
    streamline_seeds.y = linspace(min(y(:)), max(y(:)), 7);
    [streamline_seeds.x, streamline_seeds.y] = meshgrid(streamline_seeds.x, streamline_seeds.y);
    
    speed = sqrt(u.^2 + v.^2);
    
    figure;
    hold on;
    contourf(x, y, speed, 30, 'LineColor', 'none');
    colorbar;
    colorbar_obj = colorbar();
    colorbar_obj.Label.String = '(u^2+v^2)/U_\infty';
    streamline(x, y, u, v, streamline_seeds.x, streamline_seeds.y);
    xlabel('x');
    ylabel('y');
    title('Streamlines Over Velocity Magnitude');

    % Seeing Speed profile inside gap
    figure;
    contourf(x, y, u, 30, 'LineColor', 'none');
    colorbar;
    colorbar_obj = colorbar();
    colorbar_obj.Label.String = 'u/U_\infty';
    xlabel('x');
    ylabel('y');
    title('Speed Profile Over Velocity Magnitude');
    
    u_L_05 = u(:, 50);
    y_L_05 = y(:, 50);
    rho_L_05 = rho(:, 50);
    
    hold on;
    plot(u_L_05 + 12.5, y_L_05, 'r'); % Plus 12.5 is for proper placement in gap (center of length)
    plot([12.5 12.5], [-10 2], 'r');
    % plot([0 25], [0.8 0.8], '--r');
    % text(15, 1.2, 'u/U_\infty=0.5', Color='r');
    % plot([0 25], [0 0], '--r');
    % text(15, -0.4, 'u/U_\infty=0.26', Color='r');

    % Finding inflection point
    du_dy = gradient(u_L_05, y_L_05);
    d2u_dy2 = gradient(du_dy, y_L_05);
    zero_cross_indices = find(diff(sign(d2u_dy2)));
    
    figure;
    hold on;
    plot(u_L_05, y_L_05, 'r');
    plot(u_L_05(zero_cross_indices(end)), y_L_05(zero_cross_indices(end)), 'or');
    text(u_L_05(zero_cross_indices(end)) - 0.2, y_L_05(zero_cross_indices(end)) - 1, append("u:", string(u_L_05(zero_cross_indices(end))), " y:", string(y_L_05(zero_cross_indices(end)))));
    xlabel('u/U_\infty');
    ylabel('y');

    % finding compressible inflection point
    rho_du_dy = rho_L_05 .* du_dy;
    drho_du_dydy = gradient(rho_du_dy, y_L_05);
    compressible_zero_cross_indices = find(diff(sign(drho_du_dydy)));

    figure;
    hold on;
    plot(u_L_05, y_L_05, 'r');
    plot(u_L_05(compressible_zero_cross_indices(end)), y_L_05(compressible_zero_cross_indices(end)), 'or');
    text(u_L_05(compressible_zero_cross_indices(end)) - 0.2, y_L_05(compressible_zero_cross_indices(end)) - 1, append("u:", string(u_L_05(compressible_zero_cross_indices(end))), " y:", string(y_L_05(compressible_zero_cross_indices(end)))));
    xlabel('u/U_\infty');
    ylabel('y');
    title('With Compressible Correction')
end

save("results.mat");