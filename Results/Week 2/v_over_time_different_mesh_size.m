close all
clear all

load('recorded_times.txt')
recorded_times2=(0:1:500);

load('cai_centre120.txt')
load('cai_corner120.txt')
load('v_corner120.txt')
load('v_centre120.txt')

load('cai_centre.txt')
load('cai_corner.txt')
load('v_corner.txt')
load('v_centre.txt')

load('cai_centre50.txt')
load('cai_corner50.txt')
load('v_corner50.txt')
load('v_centre50.txt')

load('cai_centre20.txt')
load('cai_corner20.txt')
load('v_corner20.txt')
load('v_centre20.txt')

load('cai_centre_adjusted_Cal100.txt')
load('cai_corner_adjusted_Cal100.txt')
load('v_corner_adjusted_Cal100.txt')
load('v_centre_adjusted_Cal100.txt')

subplot(1,2,1)
hold on

plot(recorded_times2,cai_centre50, 'Linewidth', 3)
plot(recorded_times2,cai_centre_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,cai_centre120, 'Linewidth', 3)

plot(recorded_times2,cai_corner50, 'Linewidth', 3)
plot(recorded_times2,cai_corner_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,cai_corner120, 'Linewidth', 3)

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

plot(recorded_times2,v_centre50, 'Linewidth', 3)
plot(recorded_times2,v_centre_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,v_centre120, 'Linewidth', 3)
plot(recorded_times2,v_corner50,'Linewidth', 3)
plot(recorded_times2,v_corner_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,v_corner120, 'Linewidth', 3)


h= legend('$N=50$ (centre)','$N=100$ (centre)','$N=120$ (centre)', ...
    '$N=50$ (corner)','$N=100$ (corner)','$N=120$ (corner)');
set(h, 'FontSize',20,'Interpreter','Latex')
%xlim([55 65])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h1=suptitle('Uniform triangular mesh consisting of $2\times N\times N$ triangles');
set(h1,'FontSize',30,'Interpreter', 'Latex')
shg

