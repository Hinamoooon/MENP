% Multipole Expansion for NanoPhotonics (MENP)
% Exact multipole expansion (beyond long-wavelength approximation)
% is calculated from current density distribution.[ref.1]
%
% Input properties
% ENxyzf.mat: Electric fields (Ex,Ey,Ez), refractive indices (n_x,n_y,n_z),
%             position vector to each mesh grid point (x,y,z), and
%             frequency (f) exported from EM simulation software such as
%             Lumerical (FDTD), Comsol (FEM), etc.
%
% Output properties
% Cp,Cm,CQe,CQm: Scattering cross section from each of multipole moment.
% Cp: electric dipole
% Cm: magnetic dipole
% CQe: electric quadrupole
% CQm: magnetic quadrupole
%
% References
% 1. "An electromagnetic multipole expansion beyond the long-wavelength approximation"
%    (http://dx.doi.org/10.1016/j.optcom.2017.08.064)
%
% MENP (Multipolar Expansion for NanoPhotonics)
% T. Hinamoto (Kobe University, Japan)

clear all;
close all;

%% loading
addpath('../MENP');
PhysConst;
load('ENxyzf.mat');
[Cp,Cm,CQe,CQm,Csum] = exactME(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z);

%% plot
lbdp = c./f*1e9;
result = [Cp,Cm,CQe,CQm,Csum];
fig = figure();
plot(lbdp,result(:,5)*1e18,'--k'); hold on; % sum
plot(lbdp,result(:,1:4)*1e18,'-o');
xlabel('Wavelength (nm)');
ylabel('Scattering cross section (nm^2)');
xlim([min(lbdp),max(lbdp)]);

%% save
out = [lbdp,result];
outfilename = [mfilename,'.csv'];
csvwrite(outfilename,out);