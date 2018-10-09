% GRADPHI Explicit gradient approximation 
%
% grad(Phi) is approximated using first-order difference approximations.
% Forward or backward difference is chosen such that information is propagated
% away from the zero-level-set.
%
% At convergence points (regions where level sets converge), we set
% grad(Phi) = 0
%
%
% INPUTS:
%
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
% OUTPUTS:
%
%   - PX,PY,PZ = approximate local gradient in the x, y, and z direction
%                respectively
%
% Coordinate system:
%
% z = across fracture aperture
% y = parallel to mean flow direction
% x = perpendicular to mean flow direction
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%
function [PX,PY,PZ] = gradPhi(phi,dx,dy,dz,BC)

% only works if grid aspect ratio is 1 x 1 x 1
phi = phi./dx;  % if aspect ratio is not 1 x 1 x 1, work dx into differencing step


[nz,nx,ny] = size(phi);

i = 1 : 1 : nz; ip1 = circshift(i,-1); im1 = circshift(i,1);
j = 1 : 1 : nx; jp1 = circshift(j,-1); jm1 = circshift(j,1);
k = 1 : 1 : ny; kp1 = circshift(k,-1); km1 = circshift(k,1);

% calculate forward differences:
fx = phi(i,jp1,k) - phi(i,j,k);
fy = phi(i,j,kp1) - phi(i,j,k);
fz = phi(ip1,j,k) - phi(i,j,k);

idx = fx < 0;
idy = fy < 0; idy(:,:,1) = 1; idy(:,:,ny) = 0;
idz = fz < 0; idz(1,:,:) = 1; idz(nz,:,:) = 0;

PX = idx.*fx + (1-idx).*fx(:,jm1,:);
PY = idy.*fy + (1-idy).*fy(:,:,km1);
PZ = idz.*fz + (1-idz).*fz(im1,:,:);

% apply boundary conditions:
if strcmp(BC,'No Flow')
    PX(:,1,:) = phi(:,2,:)-phi(:,1,:);
    PX(:,nx,:) = phi(:,nx,:)-phi(:,nx-1,:);    
end

% set grad = 0 @ nodes where level sets converge
convgx = idx(:,j,:) == 0 & idx(:,jm1,:) == 1; 
convgz = idz(i,:,:) == 0 & idz(im1,:,:) == 1; convgz(1,:,:) = 0; convgz(nz,:,:) = 0;
convgy = idy(:,:,k) == 0 & idy(:,:,km1) == 1; convgy(:,:,1) = 0; convgy(:,:,ny) = 0;


if strcmp(BC,'No Flow')
    convgx(:,1,:) = 0; convgx(:,nx,:) = 0;
end

PX(convgx) = 0; PY(convgy) = 0; PZ(convgz) = 0;


