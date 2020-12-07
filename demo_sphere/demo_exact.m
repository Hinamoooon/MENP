% Multipole Expansion for NanoPhotonics (MENP)
% demo_exact.m is a demo file for calculation of exact multipole expansion
% from current density distribution.[ref.1]
%
% Input data
% ENxyzf.mat: Electric fields (Ex,Ey,Ez), refractive indices (n_x,n_y,n_z),
%             position vector to each mesh grid point (x,y,z), and
%             frequency (f) exported from EM simulation software such as
%             Lumerical (FDTD), Comsol (FEM), etc.
%
% Output properties
% Cp,Cm,CQe,CQm,Csum: Partial and total scattering cross section
% Cp: electric dipole
% Cm: magnetic dipole
% CQe: electric quadrupole
% CQm: magnetic quadrupole
% Csum: total (sum of above four quantities)
%
% References
% 1. "An electromagnetic multipole expansion beyond the long-wavelength approximation"
%    (http://dx.doi.org/10.1016/j.optcom.2017.08.064)
%
% MENP (Multipole Expansion for NanoPhotonics)
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
plot(lbdp,result(:,1:4)*1e18,'-o'); hold on;
plot(lbdp,result(:,5)*1e18,'--k'); % sum
title('Scattering spectra of a silicon nanosphere (D = 200 nm)', 'FontSize',12);
xlabel('Wavelength (nm)', 'FontSize',12);
ylabel('Scattering cross section (nm^2)', 'FontSize',12);
legend('C^{p}_{\rm sca}','C^{m}_{\rm sca}', ...
       'C^{Q^e}_{\rm sca}','C^{Q^m}_{\rm sca}', ...
       'C^{Total}_{\rm sca}','FontSize', 12);
xlim([min(lbdp),max(lbdp)]);

%% save
print(fig, mfilename, '-dpng');

out = [lbdp,result];
outfilename = [mfilename,'.csv'];
csvwrite(outfilename,out);