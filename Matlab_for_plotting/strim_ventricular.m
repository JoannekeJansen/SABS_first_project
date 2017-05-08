close all
clear all

load('recorded_times.txt')
load('v_left_g_Na_1.0.txt')
load('v_right_g_Na_1.0.txt')
load('cai_left_g_Na_1.0.txt')
load('cai_right_g_Na_1.0.txt')

subplot(1,2,1)
hold on
plot(recorded_times,cai_left_g_Na_1_0, 'Linewidth', 3)
plot(recorded_times,cai_right_g_Na_1_0, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

%xlim([0 1])
subplot(1,2,2)
hold on
plot(recorded_times,1000*v_left_g_Na_1_0, 'Linewidth', 3)
plot(recorded_times,1000*v_right_g_Na_1_0, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
%xlim([0 1])

h= legend('Left','Right');
set(h, 'FontSize',20,'Interpreter','Latex')

h1=suptitle('Ventricular');
set(h1,'FontSize',30,'Interpreter', 'Latex');

figure

load('recorded_times2.txt')
load('v_left_g_Na_1.02.txt')
load('v_right_g_Na_1.02.txt')
load('cai_left_g_Na_1.02.txt')
load('cai_right_g_Na_1.02.txt')

subplot(1,2,1)
hold on
plot(recorded_times2,cai_left_g_Na_1_02, 'Linewidth', 3)
plot(recorded_times2,cai_right_g_Na_1_02, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

%xlim([0 1])
subplot(1,2,2)
hold on
plot(recorded_times2,1000*v_left_g_Na_1_02, 'Linewidth', 3)
plot(recorded_times2,1000*v_right_g_Na_1_02, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
%xlim([0 1])

h= legend('Left','Right');
set(h, 'FontSize',20,'Interpreter','Latex')

h1=suptitle('Ventricular');
set(h1,'FontSize',30,'Interpreter', 'Latex');

shg
