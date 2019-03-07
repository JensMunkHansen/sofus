%
% @file   matlab/singlefocusrectarray.m
% @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
% @date   Tue Jun  7 23:52:37 2018
%
% Copyright 2018, Jens Munk Hansen

%{
//! [LinearArrayMatlab example] 
%}
fnmroot = [pwd filesep '..'];
% Build directory
builddir = [fnmroot filesep 'release'];

% Shared object name and header
if isunix; libname = 'libfnm.so'; else; libname = 'fnm.dll'; end;
libnamefull = [builddir filesep 'fnm' filesep libname];
libheader   = [fnmroot filesep 'fnm' filesep 'if_matlab.h'];

% Load library using alias fnm - only one extra header is allowed
[notfound, warnings] = loadlibrary(libnamefull, libheader,...
                                   'addheader', 'fnm_types.h',...
                                   'addheader', 'if_fnm.h',...
                                   'alias', 'fnm',...
                                   'includepath', builddir,...
                                   'includepath', fnmroot,...
                                   'mfilename', 'fnmM');

[methodinfo, structs, enuminfo, ThunkLibName] = fnmM;

c  = single(1500);
f0 = single(1.0e6);

nElements = 128;
kerf      = 5.0e-4;
width     = 3e-3;
height    = 50e-3;

nDiv = 20;

nx = 170;
nz = 250;

d = (width + kerf) * nElements;
dx = (1.5 * d) / nx;
dz = (2.0 * d) / nz;

focus = zeros(1, 3, 'single');
focus(3) = d;

wx = (nx-1.0) / 2;
xs = ([0:nx-1] - wx)*dx;
zs = [0:nz-1] * dz;

% libpointer for aperture
pAperture = libpointer;

% Create an aperture
calllib('fnm','ApertureLinearCreate',...
        pAperture, nElements, width, kerf, height);

% Set parameters
pf0 = libpointer('singlePtr', f0);
calllib('fnm', 'ApertureRwFloatParamSet', pAperture, int32(9), pf0, ...
        uint64(0), uint64(1));
pc = libpointer('singlePtr', c);
calllib('fnm', 'ApertureRwFloatParamSet', pAperture, int32(12), pc, ...
        uint64(0), uint64(1));
calllib('fnm','ApertureNDivHSet',pAperture, nDiv);
calllib('fnm','ApertureNDivWSet',pAperture, nDiv);
calllib('fnm','ApertureNThreadsSet',pAperture, uint64(4));

nDim = 3;
nPos = nx * nz;

[zs, xs] = meshgrid(zs,xs);
pos = [xs(:), zeros(nPos,1), zs(:)];
pos = single(pos');

nOut = uint64(0);
pdataOut = libpointer;
pnOut = libpointer('uint64Ptr', nOut);


pFocus = libpointer('singlePtr', focus);

% Set focus point
calllib('fnm', 'ApertureRwFloatParamSet', pAperture, int32(5), pFocus, ...
        uint64(1), uint64(3));
% Set focusing type to Rayleigh
calllib('fnm', 'ApertureFocusingTypeSet', pAperture, int32(0));

tic
% Compute CW field
calllib('fnm', 'ApertureCalcCwFieldFast', pAperture, pos, uint64(nPos), uint64(nDim),...
        pdataOut, pnOut);
toc

% Specify pointer type and dimension for output - two samples per complex
setdatatype(pdataOut,'singlePtr', 2, pnOut.Value);

% Convert and reshape
c = pdataOut.Value;
c = complex(c(1,:),c(2,:));
c = reshape(c, nx, nz);

% Plot absolute value
figure(1)
imagesc(abs(c));

% Free libpointers and library
calllib('fnm', 'ApertureDestroy', pAperture);
clear pAperture;
calllib('fnm', 'FreeCArray', pdataOut);
clear pc;
clear pf0;
clear pnOut;

unloadlibrary('fnm');
%{
//! [LinearArrayMatlab @example] 
%}

