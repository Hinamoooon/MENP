% demo_toroidal.m show how to compute multipoles including toroidal
% dipole moment. For comparison with a following reference, scattering
% properties of a silicon nanodisk with D=310 and H=50 nm are calculated.
% "Nonradiating anapole modes in dielectric nanoparticles"
% (http://www.nature.com/doifinder/10.1038/ncomms9069)

clear all;
close all;

%% loading
addpath('../MENP');
PhysConst;
load('ENxyzf.mat');
[Cp,CT,Cm,CQe,CQm,Csum] = toroidalME(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z);
[arg_px,~,~,arg_ikTx,~,~] = toroidalME_phase(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z);

%% plot
lbdp = c./f*1e9;
result = [Cp,CT,Cm,CQe,CQm,Csum];
result_phase = [arg_px,arg_ikTx];

fig1 = figure();
plot(lbdp,result(:,1:5)*1e18,'o-'); hold on;
plot(lbdp,result(:,6)*1e18,'k--');  % Csum
% legend('ED','TD','MD','QE','QM','Total');
legend('C^{p}_{\rm sca}', 'C^{T}_{\rm sca}', 'C^{m}_{\rm sca}', ...
       'C^{Q^e}_{\rm sca}','C^{Q^m}_{\rm sca}', ...
       'C^{Total}_{\rm sca}','FontSize', 12);
title('Scattering spectra of a silicon nanodisk', 'FontSize',12);
xlabel('Wavelength (nm)', 'FontSize',12);
ylabel('Scattering cross section (nm^2)', 'FontSize',12);
xlim([min(lbdp),max(lbdp)]);

fig2 = figure();
plot(lbdp,result_phase(:,1:2)); hold on;
legend('p_x','-ikT_x', 'FontSize',12);
title('Phase of electric and toroidal dipole moments, a silicon nanodisk', 'FontSize',12);
xlabel('Wavelength (nm)', 'FontSize',12);
ylabel('Phase (rad)', 'FontSize',12);
xlim([min(lbdp),max(lbdp)]);

%% save
out1 = [lbdp,result];
outfilename = [mfilename,'.csv'];
csvwrite(outfilename,out1);

out2 = [lbdp,result_phase];
outfilename = [mfilename,'_phase.csv'];
csvwrite(outfilename,out2);