close all
clear all

load('single_cell_times_atrial_stim.txt')
load('single_cell_cai_atrial_stim.txt')
load('single_cell_v_atrial_stim.txt')

subplot(1,2,1)

plot(single_cell_times_atrial_stim,single_cell_cai_atrial_stim, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
xlim([0 5])
subplot(1,2,2)

%plot(single_cell_times_atrial_stim(1:10000),1000*single_cell_v_atrial_stim(1:10000), 'Linewidth', 3)
plot(single_cell_times_atrial_stim,1000*single_cell_v_atrial_stim, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 5])

h1=suptitle('Atrial, stimulus');
set(h1,'FontSize',30,'Interpreter', 'Latex');

figure
load('single_cell_times_ventricular_stim.txt')
load('single_cell_cai_ventricular_stim.txt')
load('single_cell_v_ventricular_stim.txt')

subplot(1,2,1)

plot(single_cell_times_ventricular_stim,single_cell_cai_ventricular_stim, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
xlim([0 5])

subplot(1,2,2)

%plot(single_cell_times_ventricular_stim(1:10000),1000*single_cell_v_ventricular_stim(1:10000), 'Linewidth', 3)
plot(single_cell_times_ventricular_stim,1000*single_cell_v_ventricular_stim, 'Linewidth', 3)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 5])

h1=suptitle('Ventricular, stimulus');
set(h1,'FontSize',30,'Interpreter', 'Latex');


shg
