clc;
fprintf('===========================[ SingleFocusRectArray.m ]===========================\n\n');
fprintf('This script uses the Fast Nearfield Method to calculate the CW pressure field of\n');
fprintf('an array of 128 rectangular elements focused at a single point. The script\n');
fprintf('outputs the pressure field and a diagram of the array.\n\n');

tic
% demo file for the multiple focus array simulation
% constant parameters
f0 = 1.0; % excitation frequency,Hz
soundspeed = 2*pi; % m/s
lambda = soundspeed / f0; % wavelength, m

%define a transducer structure/array
nelex = 1;
neley = 1;
kerf = 5.0e-4;

d = nelex * (width+kerf)

xdcr_array = create_rect_planar_array(1,1,1.0,1.0,0.0,0,[0 0 0])

% create the data structure that specifies the attenuation value, etc.
lossless = set_medium('lossless');

% define the computational grid
ymin = 0;
ymax = 0;
zmin = 0.1;
zmax = 2*d;

nx = 170;
nz = 250;
dx = 0.05;
xmin = -(nx-1.0)/2.0 * dx;
xmax = (nx-1 - (nx-1.0)/2.0) * dx;

dz = 0.0;
x = xmin:dx:xmax;
z = zmin:dz:zmax;

ps = set_coordinate_grid([dx 0 dz],xmin,xmax,ymin,ymax,zmin,zmax);

ndiv = 20;

% generate the pressure field
tic
pressure = cw_pressure(xdcr_array,ps,lossless,ndiv,f0,'fnm');
toc

