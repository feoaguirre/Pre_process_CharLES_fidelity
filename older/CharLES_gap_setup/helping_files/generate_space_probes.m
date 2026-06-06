function probe_coordinates = generate_space_probes(x1, x2, y1, y2, num_of_probes_along_x, num_of_probes_along_y, filename)
    x1 = double(x1); x2 = double(x2); y1 = double(y1); y2 = double(y2); 
    num_of_probes_along_x = double(num_of_probes_along_x); 
    num_of_probes_along_y = double(num_of_probes_along_y);
    
    x = (linspace(x1, x2, num_of_probes_along_x)' +  zeros(1, num_of_probes_along_y))';
    y = linspace(y1, y2, num_of_probes_along_y)' + zeros(1, num_of_probes_along_x);
    
    x = reshape(x, [], 1);
    y = reshape(y, [], 1);
    z = zeros(length(x), 1);

    probe_coordinates = table(x, y, z);
    writetable(probe_coordinates, filename);
end



