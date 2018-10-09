% HEAVIAREA Heaviside function that separates domain into solid and 
%           liquid nodes. Segragated domain is summed along the
%           z-direction to generate a 2D array of mineral volumes
%
% INPUTS:
%   - Phi = 3D function that contains the Euclidean distance from each node
%           to the nearest interface location
%
%   - cutff = threshold used to define smearing of heaviside function
%
% OUTPUTS:
%
%   - A = 2D array of volumes
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%
function A = heaviArea(phi,cutff)

idx1 = phi >= cutff;
idx2 = phi <= -cutff;
idx3 = abs(phi) < cutff;

H = zeros(size(phi));

H(idx1) = 0;
H(idx2) = 1;
H(idx3) = 0.5*(1 - phi(idx3)/cutff - (1/pi)*sin(pi*phi(idx3)/cutff));

A = squeeze(sum(H,1));