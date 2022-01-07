%**************************************************************************
% MAIN.m
% Last edited by: pjh4 Nov 2020
%
% This is the main file to run a 2D convective flow simulation using
% the stream function - vorticity method on a cylinder. The code is 
% broken up into different sections based on the intended goal.
%
%**************************************************************************
%% Cylinder - run section to initialize simulation
INIT
% TEST

cycles = 50;
batch_size = 10;
file_prefix = strcat('sim', num2str(cycles),'_vals_batch');

[grid, shape] = MAKE_GRID(grid, shape, flow, cycles, batch_size);
PLOT_VAL(grid.psi(:,:,1), 'Boundary Condition Setup for \psi Distribution');

[grid, k] = RELAX_PSI(grid, 0, "", 1);
PLOT_VAL(grid.psi(:,:,1), 'Initial \psi Distribution')

%% Run whole simulation
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

interval = 2; % change this to choose how many batches to jump by
batches = ceil(cycles/batch_size); % change this if needed
gif_file_name = strcat('sim_',num2str(cycles),'iter_', ...
    num2str(interval), 'f.gif');

makegif(gif_file_name,[-0.5,0.5]*1e2, file_prefix, ...
batches, batch_size, interval)


%% Plot velocity distribution on top/rear of cylinder

% uncomment this and change values if want to load different data
% cycles = 10;
% batch_size = 5;
% file_prefix = strcat('sim', num2str(cycles),'_vals_batch');

x_range = [115 140];
y_range = [80 120];

% 0.00 seconds
PLOT_VEL_TEMP(file_prefix, 1, batch_size, x_range, y_range);

% % 0.06 seconds
% PLOT_VEL_TEMP(file_prefix, 10, batch_size, x_range, y_range);
% 
% % 0.66 seconds
% PLOT_VEL_TEMP(file_prefix, 100, batch_size, x_range, y_range);

%% Plot temperature evolution at different grid points

% uncomment this and change values if want to load different data
% cycles = 10;
% batch_size = 5;
% file_prefix = strcat('sim', num2str(cycles),'_vals_batch');

batches = ceil(cycles/batch_size); % change this if needed

% (141, 101) is a point near the top of cylinder
[time_arr, top_temp_evol] = GET_TEMP_EVOL(file_prefix, ...
    batches, batch_size, 141, 101);

% (101, 141) is point near the rear of cylinder
[time_arr2, rear_temp_evol] = GET_TEMP_EVOL(file_prefix, ...
    batches, batch_size, 101, 141);

% (101, 300) is point further downstream of cylinder
[time_arr3, down_temp_evol] = GET_TEMP_EVOL(file_prefix, ...
    batches, batch_size, 101, 300);

figure();
hold on;

p1= plot(time_arr, top_temp_evol, 'LineWidth', 1.2);
p2= plot(time_arr, rear_temp_evol, 'LineWidth', 1.2);
p3= plot(time_arr, down_temp_evol, 'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('Temperature (K)')
title('Temperature Evolution in Flow');

legend([p1 p2 p3], 'Top of cyl. (141, 101)', 'Rear of cyl. (101, 141)', ...
    'Downstream of cyl. (101, 300)', 'Location', 'northwest');


%% Plot heat transfer over time

% uncomment this and change values if want to load different data
% cycles = 10;
% batch_size = 5;
% file_prefix = strcat('sim', num2str(cycles),'_vals_batch');

batches = ceil(cycles/batch_size); % change this if needed

[heat_out, heat_cyl, time_arr] = MEASURE_HEAT(file_prefix, ...
    batches, batch_size);

mean(heat_out)
mean(heat_cyl)
mean(heat_out(floor(length(heat_out)/2):length(heat_out)));
mean(heat_cyl(floor(length(heat_out)/2):length(heat_out)));

figure();
hold on;

p1 = plot(time_arr, heat_out, 'LineWidth', 1.1);
p2 = plot(time_arr, heat_cyl, '--', 'LineWidth', 1.5);

xlabel('Time(s)');
ylabel('q/L (W/m)')

title('Heat Transfer Per Unit Length');
legend([p1, p2], 'Outflow - Inflow', 'Cylinder surface', 'Location', 'northwest');

