% For comparison with a reference.
% "Nonradiating anapole modes in dielectric nanoparticles"
% (http://www.nature.com/doifinder/10.1038/ncomms9069)

clear all;
close all;

%% loading
addpath('../MENP');
PhysConst;
load('ENxyzf.mat');
[Cp,CT,Cm,CQe,CQm,Csum] = toroidME(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z);

%% plot
lbdp = c./f*1e9;
result = [Cp,CT,Cm,CQe,CQm,Csum];
fig = figure();
plot(lbdp,result(:,6)*1e18,'k--'); hold on; % Csum
plot(lbdp,result(:,1:5)*1e18,'o-');

xlabel('Wavelength (nm)');
ylabel('Scattering cross section (nm^2)');
xlim([550,800]);

%% save
out = [lbdp,result];
outfilename = [mfilename,'.csv'];
csvwrite(outfilename,out);