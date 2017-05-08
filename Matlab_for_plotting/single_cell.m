close all
clear all

load('single_cell_times_atrial.txt')
load('single_cell_cai_atrial.txt')
load('single_cell_v_atrial.txt')

subplot(1,2,1)

plot(single_cell_times_atrial,single_cell_cai_atrial, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,2,2)

%plot(single_cell_times_atrial(1:10000),1000*single_cell_v_atrial(1:10000), 'Linewidth', 3)
plot(single_cell_times_atrial,1000*single_cell_v_atrial, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('Atrial, spontaneous');
set(h1,'FontSize',30,'Interpreter', 'Latex');

figure
load('single_cell_times_ventricular.txt')
load('single_cell_cai_ventricular.txt')
load('single_cell_v_ventricular.txt')

subplot(1,2,1)

plot(single_cell_times_ventricular,single_cell_cai_ventricular, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,2,2)

%plot(single_cell_times_ventricular(1:10000),1000*single_cell_v_ventricular(1:10000), 'Linewidth', 3)
plot(single_cell_times_ventricular,1000*single_cell_v_ventricular, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('Ventricular, spontaneous');
set(h1,'FontSize',30,'Interpreter', 'Latex');


shg
