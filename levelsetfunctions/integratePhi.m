% INTEGRATEPHI Implicit integration of reinitialization PDE following scheme
%              presented byPeng et al (1999)
% 
% INPUTS:
%
%   - Phi = 3D function that contains the Euclidean distance from each node
%           to the nearest interface location
%
%   - dt = time step
%
%   - dx, dy, dz = grid spacing. Current implementation only accounts for
%                  dx = dy = dz (uniform grid)
%   
%   - BC = boundary condition type along boundaries parallel to mean flow
%          direction. This is needed when calculating gradient of phi. 
%          Supported BCs: No Flow or Periodic
%
%   - S = 3D array of that specifies sign of the distance function prior to
%         reinitialization
%
% OUTPUTS:
%
%   - Phi = Updated 3D function that contains the Euclidean distance from 
%           each node to the nearest interface location
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%
function phi = integratePhi(phi,dt,dx,dy,dz,S)

[b,a,f,e,d,c] = diffPhi(phi,dx,dy,dz);

ap = max(a,0).^2;
bp = max(b,0).^2;
cp = max(c,0).^2;
dp = max(d,0).^2;
ep = max(e,0).^2;
fp = max(f,0).^2;

am = min(a,0).^2;
bm = min(b,0).^2;
cm = min(c,0).^2;
dm = min(d,0).^2;
em = min(e,0).^2;
fm = min(f,0).^2;

phi = phi - dt.*max(S,0).*(sqrt(ap + bm + cp + dm + ep + fm)-1)...
      - dt.*min(S,0).*(sqrt(am + bp + cm + dp + em + fp)-1);     