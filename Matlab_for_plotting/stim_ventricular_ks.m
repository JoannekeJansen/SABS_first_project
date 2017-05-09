close all
clear all
 
FigHandle = figure('Position', [206   475   714   223]);

var='K1';
rtsm=load(char(strcat('recorded_times_g_', var, '_1.5.txt')));
vlsm=load(char(strcat('v_left_g_', var, '_1.5.txt')));
vrsm=load(char(strcat('v_right_g_', var, '_1.5.txt')));
clsm=load(char(strcat('cai_left_g_', var, '_1.5.txt')));
crsm=load(char(strcat('cai_right_g_', var, '_1.5.txt')));

rtla=load(char(strcat('recorded_times_g_', var, '_0.5.txt')));
vlla=load(char(strcat('v_left_g_', var, '_0.5.txt')));
vrla=load(char(strcat('v_right_g_', var, '_0.5.txt')));
clla=load(char(strcat('cai_left_g_', var, '_0.5.txt')));
crla=load(char(strcat('cai_right_g_', var, '_0.5.txt')));

load('recorded_times_g_Na_1.0.txt');
load('v_left_g_Na_1.0.txt');
load('v_right_g_Na_1.0.txt');
load('cai_left_g_Na_1.0.txt');
load('cai_right_g_Na_1.0.txt');

subplot(1,3,1)
hold on
colormap(parula(5))

plot(rtsm,clsm, 'Linewidth', 3, 'Color',[127,205,187]./255);
plot(recorded_times_g_Na_1_0,cai_left_g_Na_1_0, 'Linewidth', 3, 'Color',[29,145,192]./255);
plot(rtla,clla, 'Linewidth', 3, 'Color',[8,29,88]./255);

plot(rtsm,crsm, 'Linewidth', 3, 'Color',[253,141,60]./255)
plot(recorded_times_g_Na_1_0,cai_right_g_Na_1_0, 'Linewidth', 3, 'Color',[227,26,28]./255)
plot(rtla,crla, 'Linewidth', 3, 'Color',[128,0,38]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

xlim([0 1])
subplot(1,3,2)
hold on
plot(rtsm,1000*vlsm, 'Linewidth', 3, 'Color',[127,205,187]./255)
plot(recorded_times_g_Na_1_0,1000*v_left_g_Na_1_0, 'Linewidth', 3, 'Color',[29,145,192]./255)
plot(rtla,1000*vlla, 'Linewidth', 3,'Color',[8,29,88]./255)

plot(rtsm,1000*vlsm, 'Linewidth', 3, 'Color',[253,141,60]./255)
plot(recorded_times_g_Na_1_0,1000*v_right_g_Na_1_0, 'Linewidth', 3, 'Color',[227,26,28]./255)
plot(rtla,1000*vlla, 'Linewidth', 3, 'Color',[128,0,38]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 1])

h= legend('$50\%$, Left','$100\%$, Left','$150\%$, Left', ...
    '$50\%$, Right','$100\%$, Right', ...
    '$150\%$, Right');
set(h, 'FontSize',18,'Interpreter','Latex')
fig_pos = get(gca, 'position');
pos = get(h, 'position');
set(h, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

%h1=suptitle('Ventricular');
%set(h1,'FontSize',30,'Interpreter', 'Latex');

shg
