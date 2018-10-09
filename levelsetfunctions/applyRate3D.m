% APPLYRATE3D this function applies the externally calculated reaction rate to the
% fluid-rock interface (e.g. fracture surface)
%
% Input:
%
%   - phi: 3D array of distance function values that allow for
%          identification of solid-liquid interface, Gamma. 
%
%   - R: 2D array of reaction rates. A single profile will be taken and
%        placed on cells that contain Gamma
%
%
% Output:
%
%  - F = 3D matrix with reaction rates defined on fluid-rock interface
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler


function F = applyRate3D(phi,R,dx)

[nx,ny] = size(R); % note these are depth-averaged rates

[nzPhi,nxPhi,nyPhi] = size(phi);
F = zeros(nzPhi,nxPhi,nyPhi);

for y = 1 : ny
    for x = 1 : nx
        v = phi(:,x,y);
        idx = (abs(v)<0.707*dx);
        F(idx,x,y) = R(x,y);
    end
end
