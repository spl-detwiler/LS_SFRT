% FLOW_2D solve the 2D depth-averaged Stokes equation in a rectangular domain
% with constant head on two opposite sides
% 
% User specifies boundary type on the two boundaries parallel to the mean
% flow direction: BC options are 'Periodic' or 'No Flow'
%
% INPUTS:
%   T = input transmissivity field                          [L2/T]
%   Ho = user specified head at fracture inlet (left)       [L]
%   BC = boundary condition type at boundaries parallel to mean flow. 
%        - supported boundary conditions: No Flow and Periodic
%
% OUTPUTS:
%   H = head field                              
%   QX = Darcy flow in x-direction (nx+1 x ny)   [L3/T]
%   QY = Darcy flow in y-direction (nx x ny+1)   [L3/T]
%
%
% Geometry follows the following conventions:
% Coordinate system
%   --> y,j  
%    _ _ _ _
% | |1|5|9 |_|
% | |2|6|10|_|
% v |3|7|11|_|   
% x |4|8|12|_| 
% i             
% 
% flow is in y direction; equation numbering is as shown
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler

function [ h, qx, qy, r ] = flow_2d( t, ho, BC)
% specify size of real domain (i.e. excluding ghost nodes)
nx=size(t,1);ny=size(t,2); 

% initialize vectors
b=zeros(size(t)); % right-hand-side vector
h=b;              % solution vector and we use this to initialize constant h BCs

t = padarray(t,[1 1],'replicate');   % add ghost noded to impose no flux BCs
t(1,:) = 0; t(end,:) = 0;            % create no flow along domain boundaries - this will be over-written in BCs are periodic
    
% specify constant head boundaries by specifying head at locations in
% domain
h(:,1)=ho;

tic
i=2:nx+1; j=2:ny+1; % vectors that define the i,j grid 
% build the l.h.s. coefficient matrix
% offdiagonal terms with p in (+) direction and (-) in negative direction
% each offdiagonal term is built as it's own nx*ny matrix
in = 2./(1./t(i-1,j) + 1./t(i,j));
ip = 2./(1./t(i+1,j) + 1./t(i,j));
jn = 2./(1./t(i,j-1) + 1./t(i,j));
jp = 2./(1./t(i,j+1) + 1./t(i,j));

if strcmp(BC,'Periodic') || strcmp(BC,'periodic')
    in(1,:) = 2./((1./t(2,2:ny+1)) + (1./t(nx+1,2:ny+1)));
    ip(end,:) = 2./((1./t(2,2:ny+1)) + (1./t(nx+1,2:ny+1)));
    t(1,:) = t(nx+1,:);
    t(end,:) = t(2,:);
end

% diagonal terms; these are just a linear combination of the corresponding
% off-diagonal terms!
diag=-in(:)-jn(:)-jp(:)-ip(:);

% if you have periodic boundary conditions, we will now break symmetry:
% 

% right-hand side; 
% b(1,:)=-t(1,j).*h(1,:);  b(nx-2,:)=-t(nx,j).*h(nx-2,:);
b(:,1)=-t(i,1).*h(:,1);  b(:,ny-2)=-t(i,ny).*h(:,ny-2);

% assemble equations; rather than building a nx*ny coefficient matrix
% containing mostly zeros, we assemble a sparse matrix (i.e., only store
% the diagonal and 4 off-diagonals using Matlab's 'spdiags' function

if strcmp(BC,'No Flow')
    
    B=[circshift(jn(:),-nx) circshift(in(:),-1) diag(:) circshift(ip(:),1) circshift(jp(:),nx)];
    d=[-nx,-1,0,1,nx];
    neq=nx*ny;
    
elseif strcmp(BC,'Periodic')
    
    northbnd = zeros(size(in)); northbnd(1,:) = in(1,:);
    southbnd = zeros(size(ip)); southbnd(end,:) = ip(end,:);
    in(1,:) = 0; ip(end,:) = 0; 
    
    B=[circshift(jn(:),-nx) circshift(southbnd(:),-(nx-1)) circshift(in(:),-1) diag(:) circshift(ip(:),1) circshift(northbnd(:),nx-1) circshift(jp(:),nx)];
    d=[-nx,-(nx-1),-1,0,1,nx-1,nx];
    neq=nx*ny;

end
    
% this creates the sparse matrix a; you can view the full matrix that it
% represents using full(a) (don't do this with a big problem!)
a=spdiags(B,d,neq,neq); 
% a = full(a);
toc % report time required to build matrix

tic
% solve for flow field
h=a\b(:);
toc % report time required to solve nonlinear system

% check sum of residuals - should be close to zero
r=(a*h-b(:)).^2;  
r=sum(r(:));    

h = reshape(h,nx,ny);  % this is the head field for all of the real domain
h_temp=padarray(h,[1 1],0); % include ghost nodes here to calculate qx and qy fields

if strcmp(BC,'Periodic')
    h_temp(1,:) = h_temp(nx+1,:);
    h_temp(nx+2,:) = h_temp(2,:);
end
h_temp(:,1)=ho;
i=2:nx+2; j=2:ny+1;

% now calculate flows: note this assumes that dx = dy for each grid block
qx = 2./(1./t(i-1,j) + 1./t(i,j)).*(h_temp(i-1,j)-h_temp(i,j));
i=2:nx+1; j=2:ny+2;
qy=2./(1./t(i,j-1) + 1./t(i,j)).*(h_temp(i,j-1)-h_temp(i,j));
end

