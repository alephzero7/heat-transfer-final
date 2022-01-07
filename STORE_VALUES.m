%**************************************************************************
% STORE_VALUES.m
% Last edited by: pjh4 Nov 2020
%
% This function stores the values in the simulation in a .mat
%
% grid - struct, stores values such as h, dt, psi, temp, vort, u, v
% flow - struct, stores flow parameters
% shape - struct, stores data on object in the flow
% file_name - string, where simulation values are saved
% cycle - the cycle in the batch whose data will saved (typically last 
%         cycle in the batch)
%**************************************************************************

% save data values in simulation
function grid = STORE_VALUES(grid, flow, shape, file_name, cycle)

    % temporarily store 3D arrays
    psi = grid.psi;
    vort = grid.vort;
    temp = grid.temp;
    u = grid.u;
    v = grid.v;
    
    % set the variables to values from last iteration
    grid.psi = psi(:,:,cycle);
    grid.vort = vort(:,:,cycle);
    grid.temp = temp(:,:,cycle);
    grid.u = grid.u(:,:,cycle-1);
    grid.v = grid.v(:,:,cycle-1);
    
    % save variables
    save(file_name, 'grid', 'flow', 'shape');
    
    % change variables back to 3D arrays 
    grid.psi = psi;
    grid.vort = vort;
    grid.temp = temp;
    grid.u = u;
    grid.v = v;
    
    % move values from last iteration to the first 
    grid.psi(:,:,1) = grid.psi(:,:,cycle);
    grid.vort(:,:,1) = grid.vort(:,:,cycle);
    grid.temp(:,:,1) = grid.temp(:,:,cycle);
    grid.u(:,:,1) = grid.u(:,:,cycle);
    grid.v(:,:,1) = grid.v(:,:,cycle);

end