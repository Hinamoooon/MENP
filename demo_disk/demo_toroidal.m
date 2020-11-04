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
result = [Cp,CT,Cm,CQe,CQm,Csum,arg_px,arg_ikTx];

fig1 = figure();
plot(lbdp,result(:,6)*1e18,'k--'); hold on; % Csum
plot(lbdp,result(:,1:5)*1e18,'o-');
legend('Total','ED','TD','MD','QE','QM');
xlabel('Wavelength (nm)');
ylabel('Scattering cross section (nm^2)');
xlim([550,800]);

fig2 = figure();
plot(lbdp,result_phase(:,7:8)); hold on;
legend('p_x','-ikT_x');
xlabel('Wavelength (nm)');
ylabel('Phase (rad)');
xlim([550,800]);

%% save
out1 = [lbdp,result];
outfilename = [mfilename,'.csv'];
csvwrite(outfilename,out1);

out2 = [lbdp,result_phase];
outfilename = [mfilename,'_phase.csv'];
csvwrite(outfilename,out2);