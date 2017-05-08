close all
clear all

load('values0.0001.txt')
load('values0.001.txt')
load('values_0.005.txt')
load('values0.01.txt')
load('values0.05.txt')
load('values0.1.txt')
load('values1.0.txt')
load('values2.0.txt')
hold on
plot(values0_01(:,1),values1_0(:,3), 'Linewidth', 3)
plot(values0_01(:,1),values0_1(:,3), 'Linewidth', 3)
plot(values0_01(:,1),values0_05(:,3), 'Linewidth', 3)
plot(values0_01(:,1),values0_01(:,3), 'Linewidth', 3)
plot(values_0_005(:,1),values_0_005(:,3), 'Linewidth', 3)
plot(values0_01(:,1),values0_001(:,3), 'Linewidth', 3)
plot(values0_01(:,1),values0_0001(:,3), 'Linewidth', 3)
xlim([1 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h=legend('$\Delta t=1.0$','$\Delta t=0.1$','$\Delta t=0.05$','$\Delta t=0.01$','$\Delta t=0.005$','$\Delta t=0.001$','$\Delta t=0.0001$');
set(h,'FontSize',10,'Interpreter', 'Latex')
figure
hold on
plot(values0_01(:,1),values1_0(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values0_1(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values0_05(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values0_01(:,2), 'Linewidth', 3)
plot(values_0_005(:,1),values_0_005(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values0_001(:,2), 'Linewidth', 3)
plot(values0_01(:,1),values0_0001(:,2), 'Linewidth', 3)
xlim([1 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')

ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h=legend('$\Delta t=1.0$','$\Delta t=0.1$','$\Delta t=0.05$','$\Delta t=0.01$','$\Delta t=0.005$','$\Delta t=0.001$','$\Delta t=0.0001$');
set(h,'FontSize',10,'Interpreter', 'Latex')

% 
shg