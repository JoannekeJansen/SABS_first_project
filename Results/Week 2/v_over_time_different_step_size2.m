close all
clear all

load('values0.001.txt')
load('values_0.002.txt')
load('values_0.003.txt')
load('values_0.004.txt')
load('values_0.005.txt')
load('values_0.006.txt')
load('values_0.007.txt')
load('values_0.008.txt')
load('values_0.009.txt')
load('values0.01.txt')

hold on
plot(values0_01(:,1),values0_01(:,3), 'Linewidth', 3)
plot(values_0_009(:,1),values_0_009(:,3), 'Linewidth', 3)
plot(values_0_008(:,1),values_0_008(:,3), 'Linewidth', 3)
plot(values_0_007(:,1),values_0_007(:,3), 'Linewidth', 3)
plot(values_0_006(:,1),values_0_006(:,3), 'Linewidth', 3)
plot(values_0_005(:,1),values_0_005(:,3), 'Linewidth', 3)
plot(values_0_004(:,1),values_0_004(:,3), 'Linewidth', 3)
plot(values_0_003(:,1),values_0_003(:,3), 'Linewidth', 3)
plot(values_0_002(:,1),values_0_002(:,3), 'Linewidth', 3)
plot(values0_001(:,1),values0_001(:,3), 'Linewidth', 3)
xlim([0 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h=legend('$\Delta t=0.01$','$\Delta t=0.009$','$\Delta t=0.008$','$\Delta t=0.007$','$\Delta t=0.006$','$\Delta t=0.005$','$\Delta t=0.004$','$\Delta t=0.003$','$\Delta t=0.002$','$\Delta t=0.001$');
set(h,'FontSize',10,'Interpreter', 'Latex')
figure
hold on
plot(values0_01(:,1),values0_01(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_009(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_008(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_007(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_006(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_005(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_004(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_003(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values_0_002(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values0_001(:,2), 'Linewidth', 3)
xlim([1 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')

ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h=legend('$\Delta t=0.01$','$\Delta t=0.009$','$\Delta t=0.008$','$\Delta t=0.007$','$\Delta t=0.006$','$\Delta t=0.005$','$\Delta t=0.004$','$\Delta t=0.003$','$\Delta t=0.002$','$\Delta t=0.001$');
set(h,'FontSize',10,'Interpreter', 'Latex')

% 
shg