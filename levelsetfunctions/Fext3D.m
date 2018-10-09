% FEXT3D Velocity Extension step: takes a pre-defined surface velocity and
% extends in the direction of the local surface normal. 
%
% Velocity extension is run until the rates converge within a
% localized region surrounding the fracture surface. Limiting convergence 
% within a small region surrounding the interface decreases computation
% time. 
%
% Velocity extension PDE is given by Equation (17) from Jones & Detwiler (2018)
%
% INPUTS:
%   - Phi = 3D function that contains the Euclidean distance from each node
%           to the nearest interface location
%
%   - F = 3D matrix with reaction rates defined on fluid-rock interface
%         (output of applyRate3D.m)
%
%   - interface = 3D array with flagged interface nodes (needed for
%                 convergence)
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
%   - F = 3D array of surface velocities that will be used to advect phi
%         during level set integration
%
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%

function F = Fext3D(phi,F,interface,msk,dx,dy,dz,BC)

% define initial interface and region to check convergence
idx = abs(phi)<=0.707*dx;
msk = msk-idx;
msk = logical(msk);

% define convergence parameters
conv = 1000;
Finit = F;
tol = 1e-13*mean(abs(Finit(interface)));
iter = 0;

% define characteristic curve (e.g. normal lines) that velocity will extend
% across
[PX,PY,PZ] = gradPhi(phi,dx,dy,dz,BC);

PX = PX./(sqrt(PX.^2 + PZ.^2 + PY.^2));
PZ = PZ./(sqrt(PX.^2 + PZ.^2 + PY.^2));
PY = PY./(sqrt(PX.^2 + PZ.^2 + PY.^2));
S = phi./sqrt(phi.^2 + dx^2);

tic
while conv > tol
 

    [fx_fw,fx_bw,fy_fw,fy_bw,fz_fw,fz_bw] = diffPhi(F,dx,dy,dz);

    % Godunov Integration scheme: Peng 1999
    dt = 0.5*dx;
    Fp1 = F - dt*(max(S.*PX,0).*fx_bw + min(S.*PX,0).*fx_fw ...
                 + max(S.*PZ,0).*fz_bw + min(S.*PZ,0).*fz_fw...
                 + max(S.*PY,0).*fy_bw + min(S.*PY,0).*fy_fw);

    % reset interface velocity             
    Fp1(interface) = Finit(interface);

    % check convergence
    conv = sum(sum(sum((Fp1(msk) - F(msk)).^2)));
    F = Fp1;

    %      if you'd like to visualize extension step:    
     if rem(iter,100) == 0
        figure(26)
        imshow(squeeze(F(end,:,:)),[0 max(Finit(:))],'colormap',jet,'initialmagnification',200)
        pause(0.1)
        disp(conv); format long
     end
     
% check if velocity extension has become unstable     
     if conv > 1e15
         loc = find(F == max(F(:)));
         disp(['Location of explosion: ' num2str(loc(1))]);
         error('Velocity extension is blowing up!!');
     end
     iter = iter + 1;
     
end
disp(['Fext iterations: ' num2str(iter)])
display(['Fext Time: ' num2str(toc)])