%**************************************************************************
% BULK_COMPUTE.m
% Last edited by: pjh4 Nov 2020
%
% This function performs the bulk computations for the numerical stream 
% function vorticity method.
%
% grid - the struct containing data on the simulation
% flow - the struct containing data on the flow parameters
% cycle - the iteration in the simulation
% iter - number of iterations needed for Poisson equation to converge
%**************************************************************************

% perform bulk computations
function [grid, iter] = BULK_COMPUTE(grid, flow, cycle)

% update boundary condition for vorticity
grid = GET_VORT_INIT_CONDITION(grid, cycle);
disp('Finished updating vorticity boundary conditions');

% calculate velocity from previous time step
grid = BULK_VEL(grid, cycle-1);
disp('Finished computing vel');

% calculate bulk vorticity
grid = BULK_VORT(grid, flow, cycle);
disp('Finished computing vorticity');

% relax psi (poisson)
[grid, iter] = RELAX_PSI(grid, 1, "false", cycle);
disp('Finished computing psi');

% update outflow of vorticity and psi
grid = UPDATE_OUTFLOW(grid, cycle);
disp('Finished updating outflow conditions for vorticity & psi');

% calculate bulk and update outflow temperature
grid = BULK_AND_OUTFLOW_TEMP(grid, flow, cycle);
disp('Finished calculating bulk and outflow temperature');

end

% determine vorticity boundary conditions
function grid = GET_VORT_INIT_CONDITION(grid, cycle)

psi = grid.psi(:,:,cycle-1);
vort = grid.vort(:,:,cycle-1);

for ii=1:grid.rows
    
    % vorticity remains 0 at top and bottom lids
    for j = 2:grid.cols-1
        
        if grid.key(ii,j) == 1
            
            % vorticity = -2/h^2 * Sum(psi_neighbors)
            sum_psi = psi(ii-1,j) + psi(ii+1,j) + psi(ii,j-1) + psi(ii,j+1);
            vort(ii,j) = -2*sum_psi/(grid.h)^2 ;
        end
        
    end
end

grid.vort(:,:,cycle) = vort;

end

% calculate bulk velocity
function grid = BULK_VEL(grid, cycle)

u = zeros(grid.rows, grid.cols);
v = zeros(grid.rows, grid.cols);

psiValues = grid.psi(:,:,cycle);

% bulk computation for flow velocities
for ii = 2:grid.rows-1
    
    for j = 2:grid.cols-1
        
        % if grid point is on/in the object, then velocity is 0
        % due to no slip condition
        if grid.key(ii,j) == 0
            u(ii,j) = (psiValues(ii+1,j) - psiValues(ii-1,j))/(2*grid.h);
            v(ii,j) = (psiValues(ii, j-1) - psiValues(ii, j+1))/(2*grid.h);
        end
    end
end

grid.u(:,:,cycle) = u;
grid.v(:,:,cycle) = v;

end

% calculate bulk vorticity
function grid = BULK_VORT(grid,flow, cycle)

vort = grid.vort(:,:,cycle);
newVort = vort;
u = grid.u(:,:,cycle-1);
v = grid.v(:,:,cycle-1);

% bulk computation for vorticity
for ii = 2:grid.rows-1
    for j = 2:grid.cols-1
        
        % only compute vorticity for grid points not on/in surface
        if grid.key(ii,j) == 0
            
            if u(ii,j) < 0
                deltaUVortN = u(ii,j+1)*vort(ii,j+1) - u(ii,j)*vort(ii,j);
            else
                deltaUVortN = u(ii,j)*vort(ii,j) - u(ii,j-1)*vort(ii,j-1);
            end
            
            if v(ii,j) < 0
                deltaVVortN = v(ii+1,j)*vort(ii+1,j) - v(ii,j)*vort(ii,j);
                
            else
                deltaVVortN = v(ii,j)*vort(ii,j) - v(ii-1,j)*vort(ii-1,j);
                
            end
            
            delSq = (vort(ii+1,j)+vort(ii-1,j)+vort(ii,j-1)+vort(ii,j+1) -...
                4*vort(ii,j))/(grid.h^2);
            
            newVort(ii,j) = vort(ii,j) + grid.dt*(-deltaUVortN / grid.h - ...
                deltaVVortN/grid.h + flow.mu/flow.rho*delSq);
            
        end
    end
end

grid.vort(:,:,cycle) = newVort;

end

% calculate bulk temperature, and update temp at outflow
function grid = BULK_AND_OUTFLOW_TEMP(grid, flow, cycle)

temp = grid.temp(:,:,cycle-1);
newTemp = temp;
u = grid.u(:,:,cycle-1);
v = grid.v(:,:,cycle-1);

% bulk computation for temperature
for ii = 2:grid.rows-1
    
    for j = 2:grid.cols-1
        
        % only compute temp for grid points not on/in surface
        if grid.key(ii,j) == 0
            
            if u(ii,j) < 0
                u_deltaT = u(ii,j)*(temp(ii,j+1) - temp(ii,j));
            else
                u_deltaT = u(ii,j)*(temp(ii,j) - temp(ii,j-1));
            end
            
            if v(ii,j) < 0
                v_deltaT = v(ii,j)*(temp(ii+1,j) - temp(ii,j));
                
            else
                v_deltaT = v(ii,j)*(temp(ii,j) - temp(ii-1,j));
                
            end
            
            % laplacian
            delSq = (temp(ii+1,j)+temp(ii-1,j)+temp(ii,j-1)+temp(ii,j+1) -...
                4*temp(ii,j))/(grid.h^2);
            
            % calculate new temp
            newTemp(ii,j) = temp(ii,j) + grid.dt*(-u_deltaT/ grid.h - ...
                v_deltaT/grid.h + flow.alpha*delSq);
            
        end
    end
end


% calculate outflow for temperature using constant heat flux assumption
for ii = 2:grid.rows-1
    newTemp(ii,grid.cols) = 2*newTemp(ii,grid.cols-1)-newTemp(ii,grid.cols-2);
end

% update grid variable
grid.temp(:,:,cycle) = newTemp;

end


% update outflow conditions
function grid = UPDATE_OUTFLOW (grid, cycle)

psi = grid.psi(:,:,cycle);
vort = grid.vort(:,:,cycle);

% outflow condition for stream function
for ii = 2:grid.rows-1
    psi(ii,grid.cols) = 2*psi(ii,grid.cols-1)-psi(ii,grid.cols-2);
end

% outflow condition for vorticity
for ii = 2:grid.rows-1
    vort(ii,grid.cols) = vort(ii,grid.cols-1);
end

grid.psi(:,:,cycle) = psi;
grid.vort(:,:,cycle) = vort;

end

