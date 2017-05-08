close all
clear all

load('values_GNA_11.5.txt')
load('values_GNA_17.25.txt')
load('values_GNA_23.txt')
load('values_GNA_28.75.txt')
load('values_GNA_34.5.txt')

hold on
plot(values_GNA_11_5(:,1),values_GNA_11_5(:,3), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_17_25(:,3), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_23(:,3), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_28_75(:,3), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_34_5(:,3), 'Linewidth', 3)

xlim([1 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h=legend('$g\_Na = 0.5 * 23$','$g\_Na = 0.75 * 23$','$g\_Na = 23$','$g\_Na = 1.25 * 23$','$g\_Na = 1.5 * 23$');
set(h,'FontSize',10,'Interpreter', 'Latex')
figure
hold on
hold on
plot(values_GNA_11_5(:,1),values_GNA_11_5(:,2), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_17_25(:,2), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_23(:,2), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_28_75(:,2), 'Linewidth', 3)
plot(values_GNA_11_5(:,1),values_GNA_34_5(:,2), 'Linewidth', 3)

xlim([1 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')

ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h=legend('$g\_Na = 0.5 * 23$','$g\_Na = 0.75 * 23$','$g\_Na = 23$','$g\_Na = 1.25 * 23$','$g\_Na = 1.5 * 23$');
set(h,'FontSize',10,'Interpreter', 'Latex')

% 
shg