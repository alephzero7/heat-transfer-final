%**************************************************************************
% MAKE_GRID.m
% Last edited by: pjh4 Nov 2020
% Modified from file provided by sfr Oct 2020
%
% This function initializes values in the grid based on given shape
%
% grid - the struct containing data on the simulation
% flow - the struct containing data on the flow parameters
% shape - struct, stores data on object in the flow
% cycles - total number of cycles to run simulation for
% batch_size - number of iterations in each batch
%**************************************************************************

function [grid, shape] = MAKE_GRID(grid, shape, flow, cycles, batch_size)

grid.row_mat = (-(shape.cy-1):1:grid.rows-shape.cy).*ones(grid.cols,grid.rows);
grid.row_mat = transpose(grid.row_mat);
grid.col_mat = (-(shape.cx-1):1:grid.cols-shape.cx).*ones(grid.rows,grid.cols);

switch shape.name
    case 'cylinder'
        [grid, shape] = cylinder(grid,shape);
    case 'cylinder_plate'
        [grid, shape] = cylinder_plate(grid, shape);
    otherwise
        disp('Invalid shape name')
end

% do computations in batches
if cycles > batch_size
    cycles = batch_size;
end

% set up psi and vorticity grid
grid.psi = zeros(grid.rows, grid.cols, cycles+1);
grid.vort = zeros(grid.rows, grid.cols, cycles+1);

% create velocity grid
grid.u = zeros(grid.rows, grid.cols, cycles+1);
grid.v = zeros(grid.rows, grid.cols, cycles+1);

% set shape to 0
psi_init = grid.row_mat*flow.Uinf*grid.h;
psi_init(grid.key == 1) = 0;

% move initial psi into 3D array
grid.psi(:,:,1) = psi_init;
[grid, shape] = INIT_TEMP(grid, shape, cycles);

end

% cylinder
function [grid, shape] = cylinder(grid,shape)

shape.inds = find(grid.row_mat.^2 + grid.col_mat.^2 < (shape.len/2)^2);

% indicate which grid points are on/in the object
grid.key = zeros(grid.rows, grid.cols);
grid.key(shape.inds) = 1;
end

% cylinder with splitting plate
function [grid, shape] = cylinder_plate(grid,shape)

shape.inds = find(grid.row_mat.^2 + grid.col_mat.^2 < (shape.len/2)^2);

% indicate which grid points are on/in the object
grid.key = zeros(grid.rows, grid.cols);
grid.key(shape.inds) = 1;

% add a plate behind the cylinder
rear = shape.cx + floor(shape.len/2);

for ii = -floor(shape.len/20) + shape.cy: floor(shape.len/20) + shape.cy
    for j = rear : rear + shape.lenFactor * shape.len
        
        grid.key(ii,j) = 1;
        
    end
end
end

% initialize temperature
function [grid, shape] = INIT_TEMP(grid, shape, iter)

% create temperature grid (initially set everything to 300K)
grid.temp = 300*ones(grid.rows, grid.cols, iter+1);

% at object surface, set temperature to 400K
temp_init = 300*ones(grid.rows, grid.cols);
temp_init(grid.key == 1) = 400;
grid.temp(:,:,1) = temp_init;

end