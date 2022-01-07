%**************************************************************************
% makegif.m
% Last edited by: pjh4 Nov 2020
% (modified from file provided by sfr Oct 2020
%
% This function takes stored simulation data and creates a gif in the
% same folder as the main script.
%
% filename - string, what you want the gif to be called
% colour_limit - 1x2 array, min and max limits on colorbar for vorticity
% sim_value_file - string, where simulation values are saved
% batches - number of batches that gif should display
% batch_size - number of iterations in each batch
% interval - positive integer, set to greater than 1 if want to skip
%            every few batches in making the gif
%**************************************************************************

function makegif(filename,vort_colour_limit, sim_value_file,  ...
    batches, batch_size,interval)

for b = 1:interval:batches
    
    % load simulation data
    [grid, ~, ~] = LOAD_BATCH(sim_value_file, b, batch_size);
    
    % set time step
    time_step = grid.dt * interval * batch_size;
    
    % set frame in gif
    set_frame(grid, filename, vort_colour_limit, time_step, ...
        b, batch_size);
end

disp(['Saved gif to:', ' ', filename]);

end


function set_frame(grid, filename, colour_limit, time_step, ...
    batch, batch_size)

h = figure(222);
axis tight manual % this ensures that getframe() returns a consistent size

iter = batch * batch_size;
sgtitle(strcat('Iteration: ',num2str(iter)))

% plot vorticity
subplot(3,1,1);
imagesc(grid.vort(:,:));
colormap 'jet'
daspect([1 1 1]) % change aspect ratio of the figure
c = colorbar;
c.Label.String = 'Vorticity (s^{-1})';
set(gca,'FontWeight','bold')

if colour_limit == [0,0]
else
    caxis(colour_limit)
end
title('Vorticity');

% plot temperature
subplot(3,1,2);
imagesc(grid.temp(:,:));
colormap 'jet'
daspect([1 1 1])
c = colorbar;
c.Label.String = 'Temperature (K)';
set(gca,'FontWeight','bold')
caxis([300, 400]);
title('Temperature');

% plot stream function
subplot(3,1,3);
contourf(grid.psi(:,:), 30);
colormap 'jet'
daspect([1 1 1])
c = colorbar;
c.Label.String = 'Stream function (m^2/s)';
set(gca,'FontWeight','bold')
caxis([-floor(grid.rows/2), floor(grid.rows/2)]*grid.h);
title('Stream function');

% Capture the plot as an image
frame = getframe(h);
im = frame2im(frame);

[imind,cm] = rgb2ind(im,256);

% Write to the GIF File
if iter == batch_size
    imwrite(imind,cm,filename,'gif', 'DelayTime', time_step, 'Loopcount',inf);
else
    imwrite(imind,cm,filename,'gif','DelayTime', time_step, ...
        'WriteMode','append');
end

end