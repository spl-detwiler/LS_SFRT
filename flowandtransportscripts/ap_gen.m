% AP_GEN Generate a 2D random aperture field
%
% AP_GEN(NX, NY, LX, LY, H, MU, SIGMA) generates a two-dimensional aperture field of
%   dimensions LX x LY with cutoff length scales of LX and LY, a
%   Hurst roughness exponent of H, and meand and standardard deviation of
%   MU and SIGMA
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler

function ap = ap_gen(nx,ny,lx,ly,H,mu,sigma)

n1=nx; n2=ny;
xx=-n1/2:n1/2-1; yy=-n2/2:n2/2-1;
xx=xx*lx/nx; yy=yy*ly/ny;
[x, y]=meshgrid(xx,yy);
r=(x.^2+y.^2)*2*pi; %define grid of wavenumbers
r=fftshift(r);  % transform wavenumbers so peaks are at corners

% generate power spectrum with cutoff
f=2.0*sqrt(1./((1+r).^(1+H)/(n1*n2)));
p=pi*(2*rand(n2,n1)-1);

re=f.*cos(p);
im=f.*sin(p);

ap=real(ifft2(complex(re,im)));

% set mean and stdev to 0 and 1
ap=(ap-mean(ap(:)))/std2(ap)*sigma+mu;

% set all ap values less than eps equal to eps
eps=1e-8;
bin = ap>eps;
ap=ap.*double(bin);
bin=(1-bin)*eps;
ap=ap+bin;
