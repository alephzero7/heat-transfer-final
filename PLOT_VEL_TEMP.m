%**************************************************************************
% PLOT_VEL_TEMP.m
% Last edited by: pjh4 Nov 2020
%
% This function plots the velocity and temperature distribution
% at the top and rear surfaces of the object (cylinder)
%
% a third party function is used which can be found at:
% % https://www.mathworks.com/matlabcentral/fileexchange/58527-quiver-magnitude-dependent-color-in-2d-and-3d/
%
% grid - the struct containing data on the simulation
% flow - the struct containing data on the flow parameters
% shape - the struct containing data on the shape used in the simulation
% saved_results 
% filename - string, what you want the gif to be called
% colour_limit - 1x2 array, min and max limits on colorbar
% err - current error value
%**************************************************************************

function PLOT_VEL_TEMP(sim_value_file, batch, batch_size, x_range, y_range)

[grid, ~, shape] = LOAD_BATCH(sim_value_file, batch, batch_size);

time = grid.dt * (batch_size * (batch));
time = round(time, 2);

% 
x_loc = grid.col_mat(x_range(1):x_range(2),y_range(1):y_range(2));
y_loc = grid.row_mat(x_range(1):x_range(2),y_range(1):y_range(2));

figure();
sgtitle(strcat('Top of cylinder at ',{' '}, num2str(time), 's'));

%% PLOT TOP OF CYLINDER
subplot(2,1,1);
hold on;
contour(x_loc, y_loc,...
grid.key(x_range(1):x_range(2), y_range(1):y_range(2)), 1, 'k');
quiverC2D(x_loc, y_loc, ...
    grid.u(x_range(1):x_range(2), y_range(1):y_range(2)), ... 
    grid.v(x_range(1):x_range(2), y_range(1):y_range(2)), 'LineWidth', 1)

c = colorbar;
c.Label.String = 'Velocity (m/s)';

daspect([1 1 1]) % change aspect ratio of the figure
ylim([x_range(1) - shape.cx, x_range(2) - shape.cx]);

axis off;

subplot(2,1,2);

imagesc(...
    grid.temp(x_range(1):x_range(2), y_range(1):y_range(2)))
axis('on', 'image');
axis xy
c = colorbar;
c.Label.String = 'Temperature (K)';
caxis([300, 400]);
axis off;

%% PLOT REAR OF CYLINDER
figure();
sgtitle(strcat('Rear of cylinder at',{' '}, num2str(time), 's'));
subplot(1,2,1);
hold on;

x_loc = grid.col_mat(y_range(1):y_range(2),x_range(1):x_range(2));
y_loc = grid.row_mat(y_range(1):y_range(2),x_range(1):x_range(2));
contour(x_loc, y_loc, ...
grid.key(y_range(1):y_range(2), x_range(1):x_range(2)), 1, 'k');
quiverC2D(x_loc, y_loc, ...
    grid.u(y_range(1):y_range(2), x_range(1):x_range(2)), ...
    grid.v(y_range(1):y_range(2), x_range(1):x_range(2)), 'LineWidth', 1)

c = colorbar;
c.Label.String = 'Velocity (m/s)';

xlim([x_range(1) - shape.cx, x_range(2) - shape.cx]);
daspect([1 1 1]) % change aspect ratio of the figure
axis off;

subplot(1,2,2);

imagesc(grid.temp(y_range(1):y_range(2), x_range(1):x_range(2)))
axis('on', 'image');
axis xy
c = colorbar;
c.Label.String = 'Temperature (K)';
caxis([300, 400]);
axis off;

end
