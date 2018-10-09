% CALCA_3D approximate surface integral of the 3D phi distribution to generate 2D reactive
% surface area array
%
% Governing equation presented in Equation (12) of Jones & Detwiler, 2018
% 
% INPUTS:
%   - Phi = 3D function that contains the Euclidean distance from each node
%           to the nearest interface location
%
%   - dx, dy, dz = grid spacing. Current implementation only accounts for
%                  dx = dy = dz (uniform grid)
%   
%   - BC = boundary condition type along boundaries parallel to mean flow
%          direction. This is needed when calculating gradient of phi. 
%          Supported BCs: No Flow or Periodic
% 
%
% OUTPUT:
%
%   - A = 2D array of reactive surface areas
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%

function A = calcA_3D(phi,dx,dy,dz,BC)

eps = 1.5*dx;

delta = (1/2/eps) + (1/2/eps)*cos(pi*phi./eps);
idx = phi < -eps;
delta(idx) = 0;
idx = phi > eps;
delta(idx) = 0;

[PX,PY,PZ] = gradPhi(phi,dx,dy,dz,BC);

Amat = delta.*sqrt( PX.^2 + PY.^2 + PZ.^2).*(dx*dy*dz);

A = squeeze(sum(Amat,1));
