%**************************************************************************
% GET_TEMP_EVOL.m
% Last edited by: pjh4 Nov 2020
%
% This function returns the temperature evolution at a grid point
% specified by row (x) and col (y)
%
% sim_value_file - string, where simulation values are saved
% batches - number of batches that gif should display
% batch_size - number of iterations in each batch
% x - row value of the grid point
% y - col value of the grid point
% time_arr - array that stores time values to be used for plotting
% temp_arr - array that stores temperatures calues to be used for plotting
%**************************************************************************

function [time_arr, temp_arr] = GET_TEMP_EVOL(sim_value_file, ...
    batches, batch_size, x, y)

temp_arr = zeros(batches, 1);
time_arr = zeros(batches, 1);

for b = 1:batches

    [grid, ~, ~] = LOAD_BATCH(sim_value_file, b, batch_size);

%     saved_data = strcat(sim_value_file, num2str(b), '_',...
%         'size', num2str(batch_size));
% 
%     disp(strcat('Loading data from:', saved_data));
%     load(saved_data, 'grid');

    % find the temperature from the first cycle in the batch
    temp_arr(b) = grid.temp(x, y, 1);
    time_arr(b) = grid.dt * batch_size * b;

end

end
