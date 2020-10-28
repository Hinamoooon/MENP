%% loading
clear all;
close all;
lum = dlmread('MENP_SiND_D310H50_result.txt',',',2,0);
menp = csvread('demo_toroidal.csv');

%% plot
% exact
lbdp = lum(:,1);
check = abs(lbdp(:,1) - menp(:,1)) > 0.1;  % 1 if values are different
fig1 = figure();
plot(lbdp,lum(:,2),'k--'); hold on;
plot(lbdp,menp(:,7)*1e18,'k.');
plot(lbdp,menp(:,2:6)*1e18);
legend('Total','ED','TD','MD','EQ','MQ');
xlabel('Wavelength (nm)');
ylabel('Scattering cross section (nm^2)');
xlim([min(lbdp),max(lbdp)]);

fig2 = figure();
plot(lbdp,lum(:,3),'Color',[0 0.4470 0.7410]);hold on;  % p
plot(lbdp,lum(:,5),'Color',[0.8500 0.3250 0.0980]);  % T
plot(lbdp,lum(:,4),'Color',[0.9290 0.6940 0.1250]);  % m
legend('ED','TD','MD');
xlabel('Wavelength (nm)');
ylabel('Scattering intensity (arb. units)');
xlim([min(lbdp),max(lbdp)]);