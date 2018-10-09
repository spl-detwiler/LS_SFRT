clear all; close all; clc

addpath ../../levelsetfunctions/
addpath ../../flowandtransportscripts/

savedirr1 = 'LevelSetSolns/';
savedirr2 = 'PillarSolns/';
savedirr3 = 'TransportFields/';

% read input file
In = readInput('Input.txt');

% simulation parameters
% nx = In.nx;
% ny = In.ny; 
nt = In.nt;
dx = In.dx;                     % [m]
dz = In.dx;
dy = In.dx;
BC = 'Periodic';

% transport parameters              units
gr = In.g;                          % m/s2
nu = In.nu;                         % m2/s
K = In.K;                           % m/s 
Deff = In.D;                        % m2/s 
rho = In.rho;                       % kg/m3
co = In.Cin;                        % [-]
cs = In.Cs;                         % [kg/m3]
Sh = 4.86;                          % [-]
ho = 1;                             % m
Qin = (In.Qin);                     % m3/s
co = co*cs;

% surface map:
se = fitsread('apLARGE.fits');

% load initial reaction site distributionsp
phi = fitsread('PhiFlat.fits');
[nz,nx,ny] = size(phi);


% area of reaction sites:
A = calcA_3D(phi,dx,dy,dz,BC);

% cutoff values for localization of level set evolution
beta = -8*dx;
g = -12*dx;
ct = 0.707*dx;

% initialize distance fxn and time vector
t = 0;

% now calculate the fracture aperture:
b = se - calcB(phi,dx,dy,dz);

mnb = min(b(:));
if mnb < 0; elevation = elevation - mnb + 10e-6; end


tv = zeros(1,nt);
fitswrite(se,[savedirr1 'ElevationMap.fits']) % initial surface distribution (w/o rxn sites)
fitswrite(phi,[savedirr1 'phi_0000.fits']);    % initial fracture surface distribution (w/ rxn sites)
fitswrite(b,[savedirr1 'ap_0000.fits']);              % initial fracture aperture

n = 1;
Hdlss = 1; Hdlss_o = 1; Hdlss_ref = 1;

while Hdlss/Hdlss_ref < 100
    
    % intialize 2D solution
    if n == 1 
        bref = b;
        Apill = A>0;
        Apill = double(squeeze(Apill))*(dx^2);
    end

    % level set flow field
    T = (b.^3)*gr./12./nu;
    [h,qx,qy] = flow_2d(T,ho, BC);
    
    % scale for constant flow simulation
    scl = Qin./(sum(qy(:,1)));
    qx = qx*scl;
    qy = qy*scl;
    h = h*scl;
  
    if n == 1
        Hdlss_ref = mean(h(:,1));
        Hdlss = mean(h(:,1));
        Hdlss_o = Hdlss;
    else
        Hdlss = mean(h(:,1));
    end
    
    % pillar soln flow field
    Tp = (bref.^3)*gr./12./nu;
    [hp,qxp,qyp] = flow_2d(Tp,ho, BC);
    
    % scale for constant flow simulation (pillar)
    scl = Qin./(sum(qyp(:,1)));
    qxp = qxp*scl;
    qyp = qyp*scl;
    hp = hp*scl;
    
    Keff = K./(1+((2*K.*b)./(Deff*Sh)));
    Keffp = K./(1+((2*K.*bref)./(Deff*Sh)));
    
    % calculate initial Da #
    if n == 1
        Ar = b(1,1)*(dx);
        v = qy(1,1)/Ar;
        L = dx*nx;
        Da = Keff*L/(v*b(1,1)); Da = mean(Da(:))        
        Pe = (v*b(1,1))/Deff; Pe = mean(Pe(:))
    end
    
    % transport solve and rate constant
    conc = transport_2d(qx,qy,b,A,Keff,Deff,co,cs,BC);
    R = ((Keff.*(conc-cs))).*double(A>0);
    R = (R*(100.09/40.08))./rho;
    
    % transport solve and rate constant
    concp = transport_2d(qxp,qyp,bref,Apill,Keffp,Deff,co,cs,BC);
    Rp = ((Keffp.*(concp-cs))).*(double(Apill>0));
    Rp = (Rp*(100.09/40.08))./rho;
    
%     % apply rate to upper surface
    Fold = applyRate3D(phi,R,dx);
    int = abs(phi)<=ct;
    msk = abs(phi)<=-g;
    F = Fext3D(phi,Fold,int,msk,dx,dy,dz,BC);

    
    % generate cutoff mask for upper
    c = zeros(size(phi));
    cidx1 = (phi) >= beta;
    cidx2 = (phi) < g;
    cidx3 = (phi) < beta & (phi) >= g;

    c(cidx1) = 1;
    c(cidx2) = 0;
    c(cidx3) =(((abs(phi(cidx3))-abs(g)).^2).*(2*abs(phi(cidx3)) + abs(g) - 3*abs(beta)))./((abs(g)-abs(beta)).^3); 
        
    
    % choose the time step
    beff = ((sum(qy(:,1))*12*nu)/(mean(h(:,1))-mean(h(:,end)))/gr)^(1/3);
    dt1 = 0.05*beff/max(abs(R(:)));
    
    % rank and choose 5% higher than fastest seal time 
    dtv = b./abs(R);
    idt = dtv == Inf;
    dtv(idt) = [];
    dtv = sort(dtv(:));
    
    dt2 = dtv(1);
    dt3 = 0.1*dx/(max(R(:)));
    
    dt4 = min(dt1,dt2);
    dt = min(dt3,dt4)
    if dt == dt2; disp('Using local closure time'); end
    
    % calculate level set solution advection for upper surface
    [PX,PY,PZ] = gradPhi(phi,dx,dy,dz,BC);
  
    phi = phi - dt.*c.*F.*sqrt(PX.^2 + PZ.^2 + PY.^2);
    
    % calculate pillar growth solution:
    bref = bref - dt*Rp;
    bindx = bref <=1e-6; bref(bindx) = 1e-6; 
    Apill(bindx) = 0;
    
    % now calculate b from level set advection:
    b = se - calcB(phi,dx,dy,dz);
    bindx = b <=1e-6; b(bindx) = 1e-6;
    
%     figure(10)
%     subplot(2,1,1)
%     cimshow(A,[0 dx*dx])
%     subplot(2,1,2)
%     cimshow(Apill,[0 dx*dx])
%     title(n)
% 
%     figure(11)
%     subplot(2,1,1)
%     cimshow(conc,[cs co])
%     subplot(2,1,2)
%     cimshow(concp,[cs co])
%     title(n)
% 
%     figure(12)
%     subplot(2,1,1)
%     cimshow(b,[100 350].*1e-6)
%     subplot(2,1,2)
%     cimshow(bref,[100 350].*1e-6)
%     title(n)
%   
    
  
    if Hdlss/Hdlss_o > 1.01 || rem(n,25) == 0
        Hdlss_o = Hdlss;
        filename = sprintf('%04d',n);
        fitswrite(phi,[savedirr1 'phi_' filename '.fits']); 


        fitswrite(bref,[savedirr2 'ap_' filename '.fits']);
        fitswrite(qxp,[savedirr2 'qx_' filename '.fits']);
        fitswrite(qyp,[savedirr2 'qy_' filename '.fits']);
        fitswrite(concp,[savedirr2 'conc_' filename '.fits']);
        fitswrite(Rp,[savedirr2 'rate_' filename '.fits']);
        fitswrite(hp,[savedirr2 'head_' filename '.fits']);

        fitswrite(qx,[savedirr3 'qx_' filename '.fits']);
        fitswrite(qy,[savedirr3 'qy_' filename '.fits']);
        fitswrite(conc,[savedirr3 'conc_' filename '.fits']);
        fitswrite(R,[savedirr3 'rate_' filename '.fits']);
        fitswrite(h,[savedirr3 'head_' filename '.fits']);
        fitswrite(b,[savedirr3 'ap_' filename '.fits']);

        t = t + dt;
        tv(n) = t;

        fitswrite(tv,[savedirr1 'Time.fits']);
    else
        t = t + dt;
    end
    
    phi = reinitphi3D(phi,abs(phi)<=-g,dx,dy,dz,BC);
    
    % area of reaction sites on upper surface:
    A = calcA_3D(phi,dx,dy,dz,BC);    
    n = n + 1;
    
end
disp('Finished Simulation!')




