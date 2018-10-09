% CALCB Volume integral of the 3D phi distribution to produce 2D array of 
%       surface elevations
%
% Governing equation presented in Equation (11) of Jones & Detwiler, 2018
% 
% INPUTS:
%   - Phi = 3D function that contains the Euclidean distance from each node
%           to the nearest interface location
%
%   - dx, dy, dz = grid spacing. Current implementation only accounts for
%                  dx = dy = dz (uniform grid)
%
% OUTPUT:
%
%   - b = 2D array of surface elevations
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%

function [b] = calcB(phi,dx,dy,dz)

V = heaviArea(phi,1.5*dx);
% optimized:
b = V.*dz;

% long form version
% V = V.*(dx*dy*dz);
% 
% b = V./(dx*dy);
        


