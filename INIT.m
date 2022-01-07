%**************************************************************************
% INIT.m
% Last edited by: pjh4 Nov 2020
%
% This is the file used to initialize the variables used in the simulation
% of 2D convective flow on a full 201 by 500 grid.
%
%**************************************************************************
clear;
clc;
% flow.Re_d = 200;
% flow.rho = 1.1774; % density for air (300 K) [kg/m^3]
% flow.mu = 1.8462e-5; % dynamic viscosity for air [N s/m^2]
% flow.alpha = 0.2216e-4; % thermal diffusivity for air [m^2/s]
% flow.Uinf = 1; % flow velocity [m/s]

flow.Re_d = 200;
flow.rho = 1; % density for water (0C) [kg/m^3]
flow.cp = 4.225 * 10^3; % heat capacity [J kg / C]
flow.k = 0.566; % conduction factor [W/(m C)]
flow.alpha = flow.k / (flow.rho*flow.cp); % thermal diffusivity for water [m^2/s]
flow.mu = flow.rho * flow.alpha; % dynamic viscosity 
% modified from water's value to meet constraints

flow.Uinf = 1; % flow velocity [m/s]

grid.rows = 201;
grid.cols = 500;
grid.F = 1.5; % over relaxation factor
% convergence criterion (no more than 1% change for all psi)
grid.epsilon = 0.05;

shape.name = 'cylinder'; % type of shape
shape.cx = 101; % center of the shape (x coordinate)
shape.cy = 101; % center (y coordinate)
shape.len = 50; % characteristic length of shape

% choose spacing for h that will allow for convergence
if shape.len <= 20
    grid.h = flow.Re_d / 50 *(flow.mu/flow.rho) / flow.Uinf ;

else
   grid.h = flow.Re_d / shape.len * (flow.mu/flow.rho) / flow.Uinf;

end

% CFL condition
grid.dt = 1/4 * grid.h / flow.Uinf;
