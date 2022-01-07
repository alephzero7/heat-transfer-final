%**************************************************************************
% PLOT_VAL.m
% Last edited by: pjh4 Nov 2020
%
% This function plots values provided in 2D array in a filled contour map.
%
% values - 2D array, values to be plotted
% title_str - string, title of the plot
%**************************************************************************

function PLOT_VAL(values, title_str)

% plot a color contour map
figure;
hold on
colormap(gray)
daspect([1 1 1]) % change aspect ratio of the figure
contourf(values, 30);
title(title_str);

end