%**************************************************************************
% MAIN_P2.m
% Last edited by: pjh4 Nov 2020
%
% This is the main file to run a 2D convective flow simulation using
% the stream function - vorticity method on a cylinder with a splitter
% plate attached on the rear. 
%
%**************************************************************************
%% Cylinder with plate - run section to perform simulation
INIT
shape.name = 'cylinder_plate';
shape.lenFactor = 2; % change this to adjust length of splitter plate

cycles = 50;
batch_size = 10; % choose batch size less than cycles to run simulation for

file_prefix = strcat('plateL', num2str(shape.lenFactor), '_', ...
    num2str(cycles),'_vals_batch');

[grid, shape] = MAKE_GRID(grid, shape, flow, cycles, batch_size);

PLOT_VAL(grid.psi(:,:,1), 'Boundary Condition Setup for \psi Distribution');

[grid, k] = RELAX_PSI(grid, 0, "", 1);
PLOT_VAL(grid.psi(:,:,1), 'Initial \psi Distribution')

tic
start_batch = 1;

[grid, iter_arr] = BATCH_SIM(grid, flow, shape, file_prefix, cycles, ...
    start_batch, batch_size);
toc

%% Make movie of simulation
% uncomment this and change values if want to load different data
% cycles = 10;
% batch_size = 5;
% file_prefix = strcat('sim', num2str(cycles),'_vals_batch');

interval = 1;
batches = ceil(cycles/batch_size); % change this if needed

gif_file_name = strcat('plate_',num2str(cycles),'iter_', ...
    num2str(interval), 'f.gif');

makegif(gif_file_name,[-0.5,0.5]*1e2, file_prefix, ...
batches, batch_size,interval)