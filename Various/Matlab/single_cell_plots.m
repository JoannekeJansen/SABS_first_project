close all
clear all

varsc={'g_{Na}','g_{CaL}','g_{K1}','g_{Kr}','g_{to}','g_{f}','g_{Ks}'};

for i =1:length(varsc)
var=varsc(i);
vars={'1.0', '1.0', '1.0', '1.0', '1.0', '1.0', '1.0'};
varsb={'1.0', '1.0', '1.0', '1.0', '1.0', '1.0', '1.0'};
vars{i}='0.5';
varsb{i}='1.5';
FigHandle = figure('Position', [206   475   714   223]);

%FigHandle = figure('Position', [206   475   900   223]);
rtsm=load(char(strcat('recorded_times_single_cell__', vars(1), '_', vars(2), '_', vars(3), '_', vars(4), '_', vars(5), '_', vars(6), '_', vars(7), '.txt')));
vlsm=load(char(strcat('recorded_v_single_cell__', vars(1), '_', vars(2), '_', vars(3), '_', vars(4), '_', vars(5), '_', vars(6), '_', vars(7), '.txt')));
clsm=load(char(strcat('recorded_cai_single_cell__', vars(1), '_', vars(2), '_', vars(3), '_', vars(4), '_', vars(5), '_', vars(6), '_', vars(7), '.txt')));

rtla=load(char(strcat('recorded_times_single_cell__', varsb(1), '_', varsb(2), '_', varsb(3), '_', varsb(4), '_', varsb(5), '_', varsb(6), '_', varsb(7), '.txt')));
vlla=load(char(strcat('recorded_v_single_cell__', varsb(1), '_', varsb(2), '_', varsb(3), '_', varsb(4), '_', varsb(5), '_', varsb(6), '_', varsb(7), '.txt')));
clla=load(char(strcat('recorded_cai_single_cell__', varsb(1), '_', varsb(2), '_', varsb(3), '_', varsb(4), '_', varsb(5), '_', varsb(6), '_', varsb(7), '.txt')));

recorded_times=load('recorded_times_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
v_left=load('recorded_v_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
cai_left=load('recorded_cai_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');

subplot(1,2,1)
hold on

plot(rtsm(1:10000),clsm(1:10000), 'Linewidth', 3, 'Color',[127,205,187]./255);
plot(recorded_times(1:10000),cai_left(1:10000), 'Linewidth', 3, 'Color',[29,145,192]./255);
plot(rtla(1:10000),clla(1:10000), 'Linewidth', 3, 'Color',[8,29,88]./255);

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

xlim([0 1])
subplot(1,2,2)
hold on
plot(rtsm(1:10000),1000*vlsm(1:10000), 'Linewidth', 3, 'Color',[127,205,187]./255)
plot(recorded_times(1:10000),1000*v_left(1:10000), 'Linewidth', 3, 'Color',[29,145,192]./255)
plot(rtla(1:10000),1000*vlla(1:10000), 'Linewidth', 3,'Color',[8,29,88]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 1])

a=char(strcat('$50\%$ of$\: ', var,'$'));
a2=char(strcat('$100\%$ of$\:', var,'$'));
a3=char(strcat('$150\%$ of$\:', var,'$'));

h= legend(a,a2,a3);
set(h, 'FontSize',18,'Interpreter','Latex')
%fig_pos = get(gca, 'position');
%pos = get(h, 'position');
%set(h, 'position', [pos(1)+1.2*fig_pos(3) pos(2)+0*fig_pos(4) pos(3) pos(4)]);

end
shg
