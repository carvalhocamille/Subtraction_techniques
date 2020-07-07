%% ComputeHarmonicFunction.m
%
% Computes the function and its gradient for use in validating 
% the evaluation of double- and single-layer potentials.
%
% Written by C. Carvalho on 4/2/2020
% License MIT (see Readme.md)

function [ u, gradu_x, gradu_y, gradu_z ] = ComputeFunction( x, y, z, kw )

x0 = 0.1;
y0 = 0.2;
z0 = 0.3;

% compute the harmonic function

ulen = sqrt( ( x - x0 ).^2 + ( y - y0 ).^2 + ( z - z0 ).^2 );
coeff = 0.25 * pi;
u    = coeff * exp(1i * kw * ulen) ./ ulen;

% compute the components of its gradient

gradu_x = -(1./ulen - 1i *kw).*( x - x0 ) ./ ulen .* u;
gradu_y = -(1./ulen - 1i *kw).*( y - y0 ) ./ ulen .* u;
gradu_z = -(1./ulen - 1i *kw).*( z - z0 ) ./ ulen .* u;

return