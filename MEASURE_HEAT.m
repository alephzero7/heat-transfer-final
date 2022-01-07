%**************************************************************************
% MEASURE_HEAT.m
% Last edited by: pjh4 Nov 2020
%
% This function calculates heat transfer from the object (cylinder) surface 
% and from the outflow
%
% sim_value_file - string, where simulation values are saved
% batches - number of batches that gif should display
% batch_size - number of iterations in each batch
% q_outflow - heat transfer per unit length, calculated from outflow
% q_surf - heat transfer per unit length, calculated from object surface
% time_arr - array of time values, used for plotting
%**************************************************************************

% measure the heat transfer
function [q_outflow, q_surf, time_arr] = MEASURE_HEAT(sim_value_file, ...
    batches, batch_size)

q_outflow = zeros(batches,1);
q_surf = zeros(batches,1);
time_arr = zeros(batches,1);

for b = 1:batches
    % load simulation data
    [grid, flow, shape] = LOAD_BATCH(sim_value_file, b, batch_size);
    grid = FIND_ONE_AWAY(grid);
    
    % calculate heat transfer from outflow, and flow mapping circle around
    % cylinder
    [grid, q_outflow(b)] = HEAT_OUTFLOW(grid, flow);
    [grid, q_surf(b)] = HEAT_SURF(grid, shape, flow);
    
    time_arr(b) = grid.dt * batch_size * b;
end
end

% find points that are 1 away from the object surface
% (used for heat transfer calculation)
function grid = FIND_ONE_AWAY(grid)

grid.keyPlus1 = zeros(grid.rows, grid.cols);
grid.keyBound = zeros(grid.rows, grid.cols);

for ii = 2:grid.rows-1
    for j = 2:grid.cols - 1
        
        sum = grid.key(ii-1,j) + grid.key(ii+1,j) + grid.key(ii,j-1) + ...
            grid.key(ii, j+1);
        
        if grid.key(ii,j) == 0

            if sum < 4 && sum > 0
                grid.keyPlus1(ii,j) = 1;
            end
            
        else
            if sum < 4
                grid.keyBound(ii,j) = 1;
            end
        end
    end
end

end

% calculate heat transfer in outflow
function [grid, q] = HEAT_OUTFLOW(grid, flow)

u = grid.u(:,:);
temp = grid.temp(:,:);

% calculate energy already in the fluid coming in
mass_in = 0;
for ii = 2:grid.rows-1
    m_dot = flow.rho * flow.Uinf * grid.h;
    mass_in = mass_in + m_dot * flow.cp * temp(ii, 1);
end

% calculate energy in the fluid in outflow
mass_out = 0;
q_fourier = 0;
for ii = 2:grid.rows-1
    
    m_dot = flow.rho * u(ii,grid.cols-1) * grid.h;
    mass_out = mass_out + m_dot * flow.cp * (temp(ii, grid.cols));
    
    % -k dT/dx term
    q_fourier = -flow.k * (temp(ii, grid.cols) - temp(ii, grid.cols-1));
end

q = mass_out - mass_in + q_fourier;
end

% calculate heat transfer around surface
function [grid, q] = HEAT_SURF(grid, shape, flow)

temp = grid.temp(:,:);
u = grid.u(:,:);
v = grid.v(:,:);

sum_tempBound = 0;
sum_tempOuter = 0;

q_mass = 0;

for ii = 2:grid.rows-1
    for j = 2:grid.cols-1
        
        if (grid.keyBound(ii,j) == 1)
            sum_tempBound = sum_tempBound + temp(ii,j);
        end
        
        if (grid.keyPlus1(ii,j) == 1)
            sum_tempOuter = sum_tempOuter + temp(ii,j);
            
            % calculations to find contribution to heat transfer from
            % mass flow
            vel = [u(ii,j) v(ii,j)]; % velocity vector
            dx = ii - shape.cx;
            dy = j - shape.cy;
            n = [dx dy];
            unit_n = n ./ norm(n); % vector normal to object surface
            
            % velocity into/out of the surface can be determined by dot
            % product
            dot_prod = dot(vel, unit_n);
            m_dot = flow.rho * dot_prod * grid.h;
            q_mass = q_mass + m_dot * flow.cp * temp(ii, grid.cols);
            
        end
    end
end

bound_len = sum(sum(grid.keyBound));
boundPlus1 = sum(sum(grid.keyPlus1));

% take averages over the object boundary and the mapping one grid away
avgT_bound = sum_tempBound/bound_len;
avgT_plus1 = sum_tempOuter/boundPlus1;

% q/A = -k (dT/dx) Fourier Conduction
% q/L = -k (dT) / grid.h * grid.h * (# grid points on object boundary)
q_fourier = -flow.k*(avgT_plus1-avgT_bound) * bound_len; %[W/m]

q = q_fourier + q_mass;

end
