close all
clear all
 
FigHandle = figure('Position', [206   475   714   223]);

v=load('recorded_single_cell_v_values.txt');
t=load('recorded_single_cell_times.txt');
na=load('recorded_single_cell_Ina.txt');
cal=load('recorded_single_cell_Ical.txt');
kr=load('recorded_single_cell_Ikr.txt');
k1=load('recorded_single_cell_Ik1.txt');
ks=load('recorded_single_cell_Iks.txt');
i_f=load('recorded_single_cell_If.txt');
to=load('recorded_single_cell_Ito.txt');
t=linspace(0,1,10001);

subplot(1,2,2)
%hold on
plot(t,cal, 'Linewidth', 3, 'Color',[227,26,28]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$I_{\mathrm{CaL}}$ (pA/pF)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,2,1)
%hold on
plot(t,na, 'Linewidth', 3, 'Color',[227,26,28]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$I_{\mathrm{Na}}$ (pA/pF)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
%[a,b1]=max(v);
%xlim([0.0 1.0])

FigHandle = figure('Position', [206   475   714   223]);

subplot(1,2,2)
%hold on
plot(t,[k1(10)' k1(2:end-1)' k1(end-1)'], 'Linewidth', 3,'Color',[77,175,74]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$I_{\mathrm{K1}}$ (pA/pF)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,2,1)
%hold on
plot(t,kr, 'Linewidth', 3,'Color',[77,175,74]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$I_{\mathrm{Kr}}$ (pA/pF)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

FigHandle = figure('Position', [206   475   714   223]);

subplot(1,2,2)
%hold on
plot(t,to, 'Linewidth', 3,'Color',[77,175,74]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$I_{\mathrm{to}}$ (pA/pF)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')


subplot(1,2,1)
%hold on
plot(t,ks, 'Linewidth', 3,'Color',[77,175,74]./255)

xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$I_{\mathrm{Ks}}$ (pA/pF)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

FigHandle = figure('Position', [206   475   714   223]);

subplot(1,2,1)
hold on
plot(t,[i_f(1:end-1)' i_f(end-1)'], 'Linewidth', 3,'Color',[227,26,28]./255)
plot(t(i_f>0),i_f(i_f>0), 'Linewidth', 3,'Color',[77,175,74]./255)
xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex')
ylabel('$I_{\mathrm{f}}$ (pA/pF)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,2,2)
%hold on
plot(t,1000*[v(1:10000)', v(10000)],'Linewidth', 3, 'Color',[29,145,192]./255)
xlabel('$t$ (s)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
%xlim([0.0 1.0])

shg
