function [epsilon, x, dx] = build_GRIN(lambda, Nx, spatial_window, radius, extra_params)
% build_GRIN  Build a GRIN refractive index profile
% lambda - the wavelength, in um
% Nx - the number of spatial points in each dimension
% spatial_window - the total size of space in each dimension, in um
% radius - the radius of the GRIN fiber, in um
%
% extra_params.ncore_diff - the amount to add to the Sellmeier result to get n at the core
% extra_params.alpha - the shape parameter for the GRIN fiber


% Using the Sellmeier equation to generate n(lambda)
a1=0.6961663;
a2=0.4079426;
a3=0.8974794;
b1= 0.0684043;
b2=0.1162414;
b3=9.896161;

nsi=(1+a1*(lambda.^2)./(lambda.^2 - b1^2)+a2*(lambda.^2)./(lambda.^2 - b2^2)+a3*(lambda.^2)./(lambda.^2 - b3^2)).^(0.5);

% The cladding is taken as undoped, and the core is modified
nco = nsi + extra_params.ncore_diff; % core index
ncl = nsi; % cladding index

dx = spatial_window/Nx; % um

x = (-Nx/2:Nx/2-1)*dx;
[X, Y] = meshgrid(x, x);

% GRIN profile
epsilon = (max(ncl, nco - (extra_params.ncore_diff)*(sqrt(X.^2+Y.^2)/radius).^extra_params.alpha)).^2;

end