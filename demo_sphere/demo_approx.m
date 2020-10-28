clear all;
close all;

%% loading
addpath('../MENP');
PhysConst;
load('ENxyzf.mat');
[Cp,Cm,CQe,CQm,Csum] = approxME(x,y,z,f,Ex,Ey,Ez,n_x,n_y,n_z);

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