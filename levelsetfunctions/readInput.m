% READINPUT reads textfile with simulation parameters
%
% INPUT:
%   - textfile = textfile with pre-defined simulation parameters (see
%                samples in run directories)
%
% OUTPUT (in the following order):
%
%   (1) dx = grid spacing
%
%   (2) g = gravitational acceleration constant
%
%   (3) nu = fluid viscosity [L2/T]
%
%   (4) K = reaction rate constant [L/T]
%
%   (5) D = molecular diffusion coefficient [L2/T]
%
%   (6) rho = mineral solid density [M/L3]
%
%   (7) Cin = inlet concentration (fixed) [M/L3]
%
%   (8) Cs = Dissolved ion solubility concentration [M/L3]
%
%   (9) Qin = inlet flow rate [M/L3]
%
% Copyright (c) 2018 Trevor Jones and Russ Detwiler
%
function Params = readInput(textfile)

fileID = fopen(textfile,'r');
data = textscan(fileID,'%f %s %s %s %s %s %s %s %s\n');
fclose(fileID);

inputdata = data{1};
Params = struct('dx',inputdata(1),'g',inputdata(2),'nu',inputdata(3),...
    'K',inputdata(4),'D',inputdata(5),'rho',inputdata(6),...
    'Cin',inputdata(7),'Cs',inputdata(8),'Qin',inputdata(9));