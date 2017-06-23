close all
clear all
 
FigHandle = figure('Position', [206   475   714   223]);

recorded_times=load('recorded_times_12mm_strip_no_block.txt');
v_left=load('v_left_12mm_strip_no_block.txt');
v_right=load('v_right_12mm_strip_no_block.txt');
cai_left=load('cai_left_12mm_strip_no_block.txt');
cai_right=load('cai_right_12mm_strip_no_block.txt');

brecorded_times=load('recorded_times_12mm_strip_gkr_block.txt');
bv_left=load('v_left_12mm_strip_gkr_block.txt');
bv_right=load('v_right_12mm_strip_gkr_block.txt');
bcai_left=load('cai_left_12mm_strip_gkr_block.txt');
bcai_right=load('cai_right_12mm_strip_gkr_block.txt');

subplot(1,2,1)
hold on
plot(recorded_times,cai_left, 'Linewidth', 3, 'Color',[29,145,192]./255)
plot(recorded_times,cai_right, 'Linewidth', 3, 'Color',[227,26,28]./255)
plot(brecorded_times,bcai_left,':', 'Linewidth', 3, 'Color',[29,145,192]./255)
plot(brecorded_times,bcai_right,':', 'Linewidth', 3, 'Color',[227,26,28]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
%[a,b1]=max(v);
%xlim([0 1])
subplot(1,2,2)
hold on
plot(recorded_times,1000*v_left, 'Linewidth', 3, 'Color',[29,145,192]./255)
plot(recorded_times,1000*v_right, 'Linewidth', 3, 'Color',[227,26,28]./255)
plot(brecorded_times,1000*bv_left,':', 'Linewidth', 3, 'Color',[29,145,192]./255)
plot(brecorded_times,1000*bv_right,':', 'Linewidth', 3, 'Color',[227,26,28]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
%xlim([0 1])

h= legend('Left','Right');
set(h, 'FontSize',18,'Interpreter','Latex')

shg
