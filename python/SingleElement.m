clc;
fprintf('===========================[ SingleFocusRectArray.m ]===========================\n\n');
fprintf('This script uses the Fast Nearfield Method to calculate the CW pressure field of\n');
fprintf('an array of 128 rectangular elements focused at a single point. The script\n');
fprintf('outputs the pressure field and a diagram of the array.\n\n');

tic
% demo file for the multiple focus array simulation
% constant parameters
f0 = 1e6; % excitation frequency,Hz
soundspeed = 1500; % m/s
lambda = soundspeed / f0; % wavelength, m

%define a transducer structure/array
nelex = 1;
neley = 1;
kerf = 5.0e-4;

width = 3e-3; % transducer width, m
height = 50e-3; % transducer height, m

d = nelex * (width+kerf);

xdcr_array = create_rect_planar_array(nelex, neley, width, height, kerf, 0, [0 0 0]);

% create the data structure that specifies the attenuation value, etc.
lossless = set_medium('lossless');

% define the computational grid
xmin = -1.5 * d/2;
xmax = 1.5 * d/2;
ymin = 0;
ymax = 0;
zmin = 0.0;
zmax = 2*d;

nx = 50;
nz = 70;

dx = (xmax - xmin) / nx;
dz = (zmax - zmin) / nz;
x = xmin:dx:xmax;
z = zmin:dz:zmax;

ps = set_coordinate_grid([dx 0 dz],xmin,xmax,ymin,ymax,zmin,zmax);

xmtx = 0;
ymtx = 0;
zmtx = d;

ndiv = 4;

% generate the pressure field
tic
pressure = cw_pressure(xdcr_array,ps,lossless,ndiv,f0,'fnm');
toc

figure(1);
pcolor(z*100, x*100, squeeze(abs(pressure)));
shading flat;
title('Pressure at y = 0');
xlabel('z (cm)');
ylabel('x (cm)');

figure(2);
draw_array(xdcr_array);
