%**************************************************************************
% LOAD_BATCH.m
% Last edited by: pjh4 Nov 2020
%
% This function runs the simulation on batches
%
% grid - struct, stores values such as h, dt, psi, temp, vort, u, v
% flow - struct, stores flow parameters
% shape - struct, stores data on object in the flow
% sim_value_file - string, where simulation values are saved
% batch - batch whose data should be loaded
% batch_size - number of iterations in each batch
%**************************************************************************

function [grid, flow, shape] = LOAD_BATCH(sim_value_file, batch, batch_size)

saved_data = strcat(sim_value_file, num2str(batch), '_',...
    'size', num2str(batch_size));

fprintf('\n');
disp(['Loading data from:',' ', saved_data]);
load(saved_data, 'grid', 'flow', 'shape');

end