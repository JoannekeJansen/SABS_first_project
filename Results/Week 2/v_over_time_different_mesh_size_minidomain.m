close all
clear all

load('recorded_times_adjusted_Ca100_mini_dm.txt')
recorded_times2=(0:1:500);


load('v_dlc_adjusted_Cal100.txt')
load('v_trc_adjusted_Cal100.txt')
load('cai_dlc_adjusted_Cal100.txt')
load('cai_trc_adjusted_Cal100.txt')

load('v_dlc_adjusted_Cal100_dm.txt')
load('v_trc_adjusted_Cal100_dm.txt')
load('cai_dlc_adjusted_Cal100_dm.txt')
load('cai_trc_adjusted_Cal100_dm.txt')

subplot(1,2,1)
hold on

plot(recorded_times2,cai_trc_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,cai_trc_adjusted_Cal100_dm, 'Linewidth', 3)

plot(recorded_times2,cai_dlc_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,cai_dlc_adjusted_Cal100_dm, 'Linewidth', 3)


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

plot(recorded_times2,v_trc_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,v_trc_adjusted_Cal100_dm, 'Linewidth', 3)

plot(recorded_times2,v_dlc_adjusted_Cal100, 'Linewidth', 3)
plot(recorded_times2,v_dlc_adjusted_Cal100_dm, 'Linewidth', 3)



h= legend('$N=5$ (bottom left)','$N=10$ (bottom left)', ...
    '$N=5$ (top right)','$N=10$ (top right)');
set(h, 'FontSize',20,'Interpreter','Latex')
%xlim([55 65])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h1=suptitle('Minidomain: Uniform triangular mesh consisting of $2\times N\times N$ triangles');
set(h1,'FontSize',30,'Interpreter', 'Latex')
shg

