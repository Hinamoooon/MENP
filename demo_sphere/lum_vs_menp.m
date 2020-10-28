%% loading
clear all;
close all;
lum_exact = dlmread('MENP_SiNP_r100_result_exact.txt',',',2,0);
lum_approx = dlmread('MENP_SiNP_r100_result_approx.txt',',',2,0);
menp = csvread('demo_exact.csv');
menp_approx = csvread('demo_approx.csv');

%% plot
% exact
lbdp = lum_exact(:,1);
check = abs(lbdp(:,1) - menp(:,1)) > 0.1;  % 1 if values are different
fig1 = figure();
plot(lbdp,lum_exact(:,2),'k'); hold on;
plot(lbdp,lum_exact(:,3:6));
plot(lbdp,menp(:,2:5)*1e18,'.','MarkerSize',15);
plot(lbdp,menp(:,6)*1e18,'k.','MarkerSize',15);
xlabel('Wavelength (nm)');
ylabel('Scattering cross section (nm^2)');
xlim([min(lbdp),max(lbdp)]);

% approx (vs exact lumerical)
lbdp = lum_approx(:,1);
check2 = abs(lbdp(:,1) - menp_approx(:,1)) > 0.1;  % 1 if values are different
fig2 = figure();
plot(lbdp,lum_approx(:,2),'k'); hold on;
plot(lbdp,lum_approx(:,3:6));
plot(lbdp,menp_approx(:,2:5)*1e18,'.','MarkerSize',15);
plot(lbdp,menp_approx(:,6)*1e18,'k.','MarkerSize',15);
xlabel('Wavelength (nm)');
ylabel('Scattering cross section (nm^2)');
xlim([min(lbdp),max(lbdp)]);