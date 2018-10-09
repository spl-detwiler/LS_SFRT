% REINITPHI3D Reinitialization step that resets phi to a signed distance 
%             function after level set function is advanced.
%
% Governing PDE is given by equation (18) in Jones & Detwiler (2018).
%
%
% INPUTS:
%
%   - Phi = 3D function that contains the Euclidean distance from each node
%           to the nearest interface location
%
%   - msk = a 3D binary field that defines convergence tube that will
%           surround the fracture surface
%
%   - dx, dy, dz = grid spacing. Current implementation only accounts for
%                  dx = dy = dz (uniform grid)
%   
%   - BC = boundary condition type along boundaries parallel to mean flow
%          direction. This is needed when calculating gradient of phi. 
%          Supported BCs: No Flow or Periodic
%
%
% OUTPUTS:
%
%   - Phi = 3D function that has been reset such that each node reflects the
%           Euclidean distance to the nearest interface location
%            
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%
function phi = reinitphi3D(phi,msk,dx,dy,dz,BC)

idx = abs(phi)<=0.707*dx;
phiinit = phi;

msk2 = msk-idx;
msk2 = logical(msk2);

conv = 1000;
tol = 1e-8*(mean(mean(abs(phi(msk)))));
iter = 0;

tic
while conv > tol 

    [PX,PY,PZ] = gradPhi(phi,dx,dy,dz,BC);

    
    phimag = sqrt(PX.^2 + PZ.^2 + PY.^2);   
    S = phi./sqrt((phi.^2) + (phimag.^2)*dx*dx);

% Enquist Osher Integration scheme: found in Peng 1999, from Fronts
% Propagating with curvature dependent speed: Algorithms based on
% Hamilton-Jacobi formulations. J. Comp. Phys (1988)
    
    dt = 0.5*dx;
    phip1 = integratePhi(phi,dt,dx,dy,dz,S);
        
    conv = sum(sum(sum((phip1(msk2)-phi(msk2)).^2)));
    phi = phip1;
    phi(idx) = phiinit(idx);

    iter = iter + 1;
    if rem(iter,50) == 0
        figure(2184); 
        imshow(squeeze(phimag(end,:,:)),[0.9 1.1],'colormap',jet,'initialmagnification',200);
        disp(['Reinitialization convergence = ' num2str(conv)])
    end

    if iter > 1000
        conv = 0;
    end
    
    if conv > 1e10
        error('reinitialization is blowing up!')
    end
end
disp(['Reinitialization convergence = ' num2str(conv)])
disp(['Reinitialization iterations: ' num2str(iter)])
display(['Reinitialization Time: ' num2str(toc)])