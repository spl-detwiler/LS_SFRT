% DIFFPHI Calculates forward and backward difference approximations to local
% gradient of Phi. 
%
% Gradient approximations are first-order accurate in space. Used in Fext
% and reinitphi3D
%
% INPUTS:
%
%   - phi = 3D array of Euclidean distances from each node to nearest
%           interface
%
%   - dx,dy,dz = grid spacing. Currently only accounts for uniform grid
%                spacing (dx = dy = dz). 
%
% OUPUTS:
%   - fxfw, fyfw, fzfw = forward difference approximations in x, y, and z
%                        direction, respectively
% 
%   - fxbw, fybw, fzbw = backward difference approximations in x, y, and z
%                        direction, respectively
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%

function [fxfw,fxbw,fyfw,fybw,fzfw,fzbw] = diffPhi(phi,dx,dy,dz)

% only works if grid aspect ratio is 1 x 1 x 1
phi = phi./dx;

[nz,nx,ny] = size(phi);

% indexing array
i = 1 : 1 : nz; ip1 = circshift(i,-1); im1 = circshift(i,1);
j = 1 : 1 : nx; jp1 = circshift(j,-1); jm1 = circshift(j,1);
k = 1 : 1 : ny; kp1 = circshift(k,-1); km1 = circshift(k,1);

% calculate forward differences:
fxfw = phi(:,jp1,:) - phi(:,j,:);
fyfw = phi(:,:,kp1) - phi(:,:,k);
fzfw = phi(ip1,:,:) - phi(i,:,:);

% calculate backward differences:
fxbw = phi(:,j,:)-phi(:,jm1,:);
fybw = phi(:,:,k)-phi(:,:,km1);
fzbw = phi(i,:,:)-phi(im1,:,:);

% boundary conditions:
fyfw(:,:,ny) = 0; fybw(:,:,1) = 0;
fzfw(nz,:,:) = 0; fzbw(1,:,:) = 0;
