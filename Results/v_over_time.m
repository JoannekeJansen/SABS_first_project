close all
clear all

load('recorded_times.txt')

load('cai_centre.txt')
load('cai_corner.txt')
load('v_corner.txt')
load('v_centre.txt')

load('cai_centre_gna_90.txt')
load('cai_corner_gna_90.txt')
load('v_corner_gna_90.txt')
load('v_centre_gna_90.txt')

load('cai_centre_gna_110.txt')
load('cai_corner_gna_110.txt')
load('v_corner_gna_110.txt')
load('v_centre_gna_110.txt')

load('cai_centre_sigma_l_90.txt')
load('cai_corner_sigma_l_90.txt')
load('v_corner_sigma_l_90.txt')
load('v_centre_sigma_l_90.txt')

load('cai_centre_sigma_l_110.txt')
load('cai_corner_sigma_l_110.txt')
load('v_corner_sigma_l_110.txt')
load('v_centre_sigma_l_110.txt')

hold on

%plot(recorded_times,cai_centre, 'Linewidth', 3)
plot(recorded_times,cai_corner, 'Linewidth', 3, 'Color', [0.4940    0.1840    0.5560])

%plot(recorded_times,cai_centre_gna_90, 'Linewidth', 3)
plot(recorded_times,cai_corner_gna_90, ':', 'Linewidth', 3, 'Color', [0.6350    0.0780    0.1840])
%plot(recorded_times,cai_centre_gna_110, 'Linewidth', 3)
plot(recorded_times,cai_corner_gna_110,'--', 'Linewidth', 3, 'Color', [0.6350    0.0780    0.1840])

%plot(recorded_times,cai_centre_sigma_l_90, 'Linewidth', 3)
plot(recorded_times,cai_corner_sigma_l_90, ':', 'Linewidth', 3, 'Color', [ 0    0.4470    0.7410])
%plot(recorded_times,cai_centre_sigma_l_110, 'Linewidth', 3)
plot(recorded_times,cai_corner_sigma_l_110,'--', 'Linewidth', 3, 'Color', [ 0    0.4470    0.7410])

xlim([80 86])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

figure
hold on

%plot(recorded_times,v_centre, 'Linewidth', 3)
plot(recorded_times,v_corner, 'Linewidth', 3, 'Color', [0.4940    0.1840    0.5560])

%plot(recorded_times,v_centre_gna_90, 'Linewidth', 3)
plot(recorded_times,v_corner_gna_90,':','Linewidth', 3, 'Color', [0.6350    0.0780    0.1840])
%plot(recorded_times,v_centre_gna_110, 'Linewidth', 3)
plot(recorded_times,v_corner_gna_110,'--','Linewidth', 3, 'Color', [0.6350    0.0780    0.1840])

%plot(recorded_times,v_centre_sigma_l_90, 'Linewidth', 3)
plot(recorded_times,v_corner_sigma_l_90,':', 'Linewidth', 3, 'Color', [ 0    0.4470    0.7410])
%plot(recorded_times,v_centre_sigma_l_110, 'Linewidth', 3)
plot(recorded_times,v_corner_sigma_l_110,'--', 'Linewidth', 3, 'Color', [0    0.4470    0.7410])

xlim([55 65])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

figure
hold on
plot(recorded_times,cai_centre, 'Linewidth', 3)
plot(recorded_times,cai_corner, 'Linewidth', 3)

plot(recorded_times,cai_centre_gna_90, 'Linewidth', 3)
plot(recorded_times,cai_corner_gna_90, 'Linewidth', 3)
plot(recorded_times,cai_centre_gna_110, 'Linewidth', 3)
plot(recorded_times,cai_corner_gna_110, 'Linewidth', 3)

plot(recorded_times,cai_centre_sigma_l_90, 'Linewidth', 3)
plot(recorded_times,cai_corner_sigma_l_90, 'Linewidth', 3)
plot(recorded_times,cai_centre_sigma_l_110, 'Linewidth', 3)
plot(recorded_times,cai_corner_sigma_l_110, 'Linewidth', 3)

%xlim([80 86])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

figure
hold on

plot(recorded_times,v_centre, 'Linewidth', 3)
plot(recorded_times,v_corner, 'Linewidth', 3)

plot(recorded_times,v_centre_gna_90, 'Linewidth', 3)
plot(recorded_times,v_corner_gna_90,'Linewidth', 3)
plot(recorded_times,v_centre_gna_110, 'Linewidth', 3)
plot(recorded_times,v_corner_gna_110,'Linewidth', 3)

plot(recorded_times,v_centre_sigma_l_90, 'Linewidth', 3)
plot(recorded_times,v_corner_sigma_l_90, 'Linewidth', 3)
plot(recorded_times,v_centre_sigma_l_110, 'Linewidth', 3)
plot(recorded_times,v_corner_sigma_l_110, 'Linewidth', 3)
h= legend('$GNa=23$ (centre)','$GNa=23$ (corner)','$GNa=0.9*23$ (centre)', ...
'$GNa=0.9*23$ (corner)','$GNa=1.1*23$ (centre)','$GNa=1.1*23$ (corner)', ...
'$\sigma_l=0.9*0.15$ (centre)','$\sigma_l=0.9*0.15$ (corner)', ...
'$\sigma_l=1.1*0.15$ (centre)','$\sigma_l=1.1*0.15$ (corner)')
set(h, 'FontSize',20,'Interpreter','Latex')
%xlim([55 65])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
shg

shg


