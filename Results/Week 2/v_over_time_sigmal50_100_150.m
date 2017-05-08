close all
clear all

recorded_times2=(0:1:500);

load('cai_centre_adjusted_sigmal150.txt')
load('cai_corner_adjusted_sigmal150.txt')
load('v_corner_adjusted_sigmal150.txt')
load('v_centre_adjusted_sigmal150.txt')

load('cai_centre_adjusted_Cal100.txt')
load('cai_corner_adjusted_Cal100.txt')
load('v_corner_adjusted_Cal100.txt')
load('v_centre_adjusted_Cal100.txt')

load('cai_centre_adjusted_sigmal50.txt')
load('cai_corner_adjusted_sigmal50.txt')
load('v_corner_adjusted_sigmal50.txt')
load('v_centre_adjusted_sigmal50.txt')

subplot(1,2,1)
hold on

plot(recorded_times2,cai_centre_adjusted_sigmal50, 'Linewidth', 3)
plot(recorded_times2,cai_centre_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,cai_centre_adjusted_sigmal150, 'Linewidth', 3)
plot(recorded_times2,cai_corner_adjusted_sigmal50, 'Linewidth', 3)
plot(recorded_times2,cai_corner_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,cai_corner_adjusted_sigmal150, 'Linewidth', 3)

%xlim([80 86])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
% h= legend('$10$ (centre)','$10$ (corner)','$20$ (centre)', ...
%     '$20$ (corner)','$50$ (centre)','$50$ (corner)','$010$ (centre)','$100$ (corner)');
% set(h, 'FontSize',20,'Interpreter','Latex')
subplot(1,2,2)
hold on

plot(recorded_times2,v_centre_adjusted_sigmal50, 'Linewidth', 3)
plot(recorded_times2,v_centre_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,v_centre_adjusted_sigmal150, 'Linewidth', 3)
plot(recorded_times2,v_corner_adjusted_sigmal50, 'Linewidth', 3)
plot(recorded_times2,v_corner_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,v_corner_adjusted_sigmal150, 'Linewidth', 3)

h= legend('$50\%$ (centre)','$100\%$ (centre)','$150\%$ (centre)', ...
    '$50\%$ (corner)','$100\%$ (corner)','$150\%$ (corner)');
set(h, 'FontSize',20,'Interpreter','Latex')
%xlim([55 65])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h1=suptitle('$\sigma_l$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex')
shg

