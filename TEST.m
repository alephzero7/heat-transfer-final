%**************************************************************************
% TEST.m
% Last edited by: pjh4 Nov 2020
%
% This is the file used to initialize the variables used in the simulation
% of 2D convective flow on a smaller grid. (This file is outdated since
% any changes were made on INIT.m
%
%**************************************************************************

clear;
clc;

flow.Re_d = 200;
flow.rho = 1.1774; % density for air (300 K) [kg/m^3]
flow.mu = 1.8462e-5; % dynamic viscosity for air [N s/m^2]
flow.alpha = 0.2216e-4; % thermal diffusivity for air [m^2/s]
flow.Uinf = 1; % flow velocity [m/s]

grid.rows = 51;
grid.cols = 101;
grid.F = 1.5; % over relaxation factor
% convergence criterion (no more than 1% change for all psi)
grid.epsilon = 0.01;

shape.name = 'cylinder'; % type of shape
shape.cx = 21; % center of the shape (x coordinate)
shape.cy = 26; % center (y coordinate)
shape.len = 20; % characteristic length of shape (diam. for cylinder)

% choose spacing for h that will allow for convergence
if shape.len <= 20
    grid.h = flow.Re_d / 50 *(flow.mu/flow.rho) / flow.Uinf ;

else
   grid.h = flow.Re_d / shape.len * (flow.mu/flow.rho) / flow.Uinf;

end
    
grid.dt = 1/4 * grid.h / flow.Uinf;
