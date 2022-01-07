%**************************************************************************
% BATCH_SIM.m
% Last edited by: pjh4 Nov 2020
%
% This function runs the simulation on batches
%
% grid - struct, stores values such as h, dt, psi, temp, vort, u, v
% flow - struct, stores flow parameters
% shape - struct, stores data on object in the flow
% cycles_to_run - total number of cycles to run simulation for
% file_prefix - string, prefix for how the simulation data per batch 
%               will be saved
% start_batch - typically 1, change value if want to run simulation longer
%               on previously saved data
% batch_size - number of iterations in each batch
% iter_arr - array, stores how many iterations it would take for poisson
%            eqn to converge for each cycle
%**************************************************************************

function [grid, iter_arr] = BATCH_SIM(grid, flow, shape, file_prefix, cycles_to_run, ...
    start_batch, batch_size)

iter_arr = zeros(cycles_to_run, 1);

for c=2:cycles_to_run+1
    fprintf('\n');
    disp(['Running Cycle:',' ', num2str(c-1)]);
    
    ii = mod(c-2, batch_size)+2;
    
    % perform bulk computation
    [grid, iter_arr(c-1)] = BULK_COMPUTE(grid,flow, ii);
    
    if ii == batch_size + 1
        
        file_name = strcat(file_prefix, num2str(start_batch), '_',...
            'size', num2str(batch_size));
        
        fprintf('\n');
        
        disp(['Finished Batch', ' ', num2str(start_batch)]);
        disp(['Saving results to:', ' ', file_name]);
        
        % save values and shift them over in 3D array
        grid = STORE_VALUES(grid, flow, shape, file_name, ii);
        
        start_batch = start_batch + 1;
    end
    
end

% save last batch if cycles is not perfectly divisble by batch_size
if mod(cycles_to_run,batch_size) ~= 0
    
    file_name = strcat(file_prefix, num2str(batch), '_',...
        'size', num2str(batch_size), '.mat');
    
    disp(strcat('Finished batch ', num2str(batch)));
    disp(['Saving results to:', ' ', file_name]);
    save('file_name');
end


end