% TRANSPORT_2D Solve steady-state 2D groundwater
% advection-dispersion-reaction equation   
% written by RLD - adapted to a depth averaged solver by TAJ
%
% Governing eqn: div(q*c) - div(Db*grad(c)) = R
% where R = k*C
%
% Integration conventions:
%
%  Coordinate System
%  |---> y
%  |
%  V x
%
%              N 
%              |
%              |  
%              V 
%        |----------|
%        |          |
% W -->  |     P    | <-- E
%        |          | 
%        |----------|
%              ^
%              |   
%              |
%              S
%
% Inputs:                                           Units:      Array Size:
%   - qx = flow in the transverse direction         [L3/T]      [ny+1,nx]
%   - qy = flow in the mean flow direction          [L3/T]      [ny,nx+1]
%   - ap = aperture field                           [L]         [ny,nx]
%   - A = reactive surface area field               [L2]        [ny,nx]
%   - Kb = reaction rate constant                   [L/T]       [1,1]
%   - D1 = molecular diffusion coefficient          [L2/T]      [1,1]
%   - Co = inlet concentration                      [M/L3]      [1,1]
%
% Outputs: 
%   - c = concentration field                       [M/L3]      [ny,nx]
%
% Steady-state solution has been compared to a 1D analytical solution. This
% analytical soln can be found in
% '/levelsetmethods/flowandtransportscripts/cAnalytical.m'
% 
% Script that compares simlution to analytical soln can be found in:
% '/levelsetmethods/OutputFiles/AnalyticalSolnTests/TestAnalytical.m'
% 
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%
function [ c ] = transport_2d(qx, qy, ap, A, Kb, D1, co, cs, BC)
[nx, ny]=size(ap);
b=zeros(nx,ny);
c=zeros(nx,ny);
c(:,1)=co;

ap = padarray(ap,[1 1],'replicate');
D=ones(nx+2,ny+2).*D1; % padded array of dispersion coefficients
if strcmp(BC,'periodic') || strcmp(BC,'Periodic')
    % replicate aperture values from either side
    ap(1,:) = ap(nx+1,:); 
    ap(nx+2,:) = ap(2,:);
elseif strcmp(BC,'No Flow') || strcmp(BC,'no flow')
    % no diffusive flux if BC type is no-flow
    D(1,:)=0; D(nx+2,:)=0;
else
    error('unsupported BC type for transport solver')
end


i=2:nx+1; j=2:ny+1; % define indices for domain (excluding ghost nodes)
% positive Darcy fluxes into volume (i,j) through each face
qxn = qx(i-1,:);
qxp = -qx(i,:);
qyn = qy(:,j-1);
qyp = -qy(:,j);

% flag to denote upwind direction --> if = 1 neighbor is upwind 
%                                     if = 0 (i,j) is upwind
axn = qxn>0;
axp = qxp>0;
ayn = qyn>0;
ayp = qyp>0;

in = D(i-1,j).*((ap(i-1,j)+ap(i,j))./2) + qxn.*axn;
ip = D(i+1,j).*((ap(i+1,j)+ap(i,j))./2) + qxp.*axp;
jn = D(i,j-1).*((ap(i,j-1)+ap(i,j))./2) + qyn.*ayn;
jp = D(i,j+1).*((ap(i,j+1)+ap(i,j))./2) + qyp.*ayp;

    

% diagonal terms; these are just a linear combination of the corresponding
% off-diagonal terms
diag = -D(i-1,j).*((ap(i-1,j)+ap(i,j))./2)... 
    - D(i+1,j).*((ap(i+1,j)+ap(i,j))./2)... 
    - D(i,j-1).*((ap(i,j-1)+ap(i,j))./2)...
    - D(i,j+1).*((ap(i,j+1)+ap(i,j))./2)...
    + qxn.*(1-axn) + qxp.*(1-axp) + qyn.*(1-ayn) + qyp.*(1-ayp) - Kb.*A;

% right-hand side; 
b(:,:) = -Kb.*cs.*A;
b(:,1) = (-qyn(:,1).*ayn(:,1) - D(i,1).*ap(i,1)).*c(:,1); 
diag(:,ny) = diag(:,ny) + D(i,ny).*ap(i,ny);

% assemble equations; rather than building a nx*ny coefficient matrix
% containing mostly zeros, we assemble a sparse matrix (i.e., only store
% the diagonal and 4 off-diagonals using Matlab's 'spdiags' function

if strcmp(BC,'No Flow') || strcmp(BC,'no flow')
    
    B=[circshift(jn(:),-nx) circshift(in(:),-1) diag(:) circshift(ip(:),1) circshift(jp(:),nx)];
    d=[-nx,-1,0,1,nx];
    neq=nx*ny;

elseif strcmp(BC,'Periodic') || strcmp(BC,'periodic')
    northbnd = zeros(size(in)); northbnd(1,:) = in(1,:);
    southbnd = zeros(size(in)); southbnd(end,:) = ip(end,:);
    in(1,:) = 0; ip(end,:) = 0;
    
    B=[circshift(jn(:),-nx) circshift(southbnd(:),-(nx-1)) circshift(in(:),-1) diag(:) circshift(ip(:),1) circshift(northbnd(:),nx-1) circshift(jp(:),nx)];
    d=[-nx,-nx+1,-1,0,1,nx-1,nx];
    neq=nx*ny;
end

% this creates the sparse matrix a; you can view the full matrix that it
% represents using full(a) (don't do this with a big problem!)
a=spdiags(B,d,neq,neq); 
% a = full(a);
c=a\b(:);
toc % report time required to build matrix

c = reshape(c,nx,ny);  % this is the concentration field for all of the real domain
end

