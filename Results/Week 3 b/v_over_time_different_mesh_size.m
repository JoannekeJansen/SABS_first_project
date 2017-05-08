close all
clear all

load('recorded_times.txt')
recorded_times2=(0:1:500);

load('cai_left_mesh_120.txt')
load('cai_right_mesh_120.txt')
load('v_left_mesh_120.txt')
load('v_right_mesh_120.txt')

load('cai_left_mesh_100.txt')
load('cai_right_mesh_100.txt')
load('v_left_mesh_100.txt')
load('v_right_mesh_100.txt')

load('cai_left_mesh_80.txt')
load('cai_right_mesh_80.txt')
load('v_left_mesh_80.txt')
load('v_right_mesh_80.txt')

subplot(1,2,1)
hold on

plot(recorded_times2,cai_left_mesh_80, 'Linewidth', 3)
plot(recorded_times2,cai_left_mesh_100, 'Linewidth', 3)
plot(recorded_times2,cai_left_mesh_120, 'Linewidth', 3)

plot(recorded_times2,cai_right_mesh_80, 'Linewidth', 3)
plot(recorded_times2,cai_right_mesh_100, 'Linewidth', 3)
plot(recorded_times2,cai_right_mesh_120, 'Linewidth', 3)

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

plot(recorded_times2,v_left_mesh_80, 'Linewidth', 3)
plot(recorded_times2,v_left_mesh_100, 'Linewidth', 3)
plot(recorded_times2,v_left_mesh_120, 'Linewidth', 3)

plot(recorded_times2,v_right_mesh_80, 'Linewidth', 3)
plot(recorded_times2,v_right_mesh_100, 'Linewidth', 3)
plot(recorded_times2,v_right_mesh_120, 'Linewidth', 3)

h= legend('$N=80$ (left)','$N=100$ (left)','$N=120$ (left)', ...
    '$N=80$ (right)','$N=100$ (right)','$N=120$ (right)');
set(h, 'FontSize',20,'Interpreter','Latex')
%xlim([55 65])

%set(gca, 'XTick', [0:2:35])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h1=suptitle('Uniform triangular mesh consisting of $2\times N\times 1$ triangles on a $5.0\times 0.01$ mm$^2$ domain');
set(h1,'FontSize',30,'Interpreter', 'Latex')
shg

