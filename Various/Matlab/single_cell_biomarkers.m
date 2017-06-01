close all
clear all
 
varsc={'g_{Na}','g_{CaL}','g_{K1}','g_{Kr}','g_{to}','g_{f}','g_{Ks}'};

for i =1:length(varsc)
var=varsc(i);
vars={'1.0', '1.0', '1.0', '1.0', '1.0', '1.0', '1.0'};
varsb={'1.0', '1.0', '1.0', '1.0', '1.0', '1.0', '1.0'};
vars{i}='0.5';
varsb{i}='1.5';

rtsm=load(char(strcat('recorded_times_single_cell__', vars(1), '_', vars(2), '_', vars(3), '_', vars(4), '_', vars(5), '_', vars(6), '_', vars(7), '.txt')));
vlsm=load(char(strcat('recorded_v_single_cell__', vars(1), '_', vars(2), '_', vars(3), '_', vars(4), '_', vars(5), '_', vars(6), '_', vars(7), '.txt')));
clsm=load(char(strcat('recorded_cai_single_cell__', vars(1), '_', vars(2), '_', vars(3), '_', vars(4), '_', vars(5), '_', vars(6), '_', vars(7), '.txt')));

rtla=load(char(strcat('recorded_times_single_cell__', varsb(1), '_', varsb(2), '_', varsb(3), '_', varsb(4), '_', varsb(5), '_', varsb(6), '_', varsb(7), '.txt')));
vlla=load(char(strcat('recorded_v_single_cell__', varsb(1), '_', varsb(2), '_', varsb(3), '_', varsb(4), '_', varsb(5), '_', varsb(6), '_', varsb(7), '.txt')));
clla=load(char(strcat('recorded_cai_single_cell__', varsb(1), '_', varsb(2), '_', varsb(3), '_', varsb(4), '_', varsb(5), '_', varsb(6), '_', varsb(7), '.txt')));

recorded_times=load('recorded_times_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
v_left=load('recorded_v_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
cai_left=load('recorded_cai_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');

rtsm=rtsm(1:10000);
vlsm=vlsm(1:10000);
clsm=clsm(1:10000);

rtla=rtla(1:10000);
vlla=vlla(1:10000);
clla=clla(1:10000);

%%
v_ampli_50(i)=max(vlsm)-min(vlsm) %V
cai_ampli_50(i)=max(clsm)-min(clsm) %mM

v_ampli_150(i)=max(vlla)-min(vlla) %V
cai_ampli_150(i)=max(clla)-min(clla) %mM
%%
[a,b]=max(vlsm);
v_V_max_50(i)=max(vlsm(2:b)-vlsm(1:b-1))*10000 %V/s

[a,b]=max(clsm);
cai_V_max_50(i)=max(clsm(2:b)-clsm(1:b-1))*10000 %mM/s

[a,b]=max(vlla);
v_V_max_150(i)=max(vlla(2:b)-vlla(1:b-1))*10000 %V/s

[a,b]=max(clla);
cai_V_max_150(i)=max(clla(2:b)-clla(1:b-1))*10000 %mM/s
%%
a=find(vlsm>(min(vlsm)+0.1*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.1*v_ampli_50(i)), 1, 'first');
v_90_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.1*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.1*cai_ampli_50(i)), 1, 'first');
cai_90_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.1*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.1*v_ampli_150(i)), 1, 'first');
v_90_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.1*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.1*cai_ampli_150(i)), 1, 'first');
cai_90_150(i)=(b)/10 %ms
%%
a=find(vlsm>(min(vlsm)+0.3*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.3*v_ampli_50(i)), 1, 'first');
v_70_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.3*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.3*cai_ampli_50(i)), 1, 'first');
cai_70_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.3*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.3*v_ampli_150(i)), 1, 'first');
v_70_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.3*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.3*cai_ampli_150(i)), 1, 'first');
cai_70_150(i)=(b)/10 %ms
%%
a=find(vlsm>(min(vlsm)+0.5*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.5*v_ampli_50(i)), 1, 'first');
v_50_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.5*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.5*cai_ampli_50(i)), 1, 'first');
cai_50_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.5*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.5*v_ampli_150(i)), 1, 'first');
v_50_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.5*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.5*cai_ampli_150(i)), 1, 'first');
cai_50_150(i)=(b)/10 %ms
%%
a=find(vlsm>(min(vlsm)+0.7*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.7*v_ampli_50(i)), 1, 'first');
v_30_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.7*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.7*cai_ampli_50(i)), 1, 'first');
cai_30_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.7*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.7*v_ampli_150(i)), 1, 'first');
v_30_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.7*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.7*cai_ampli_150(i)), 1, 'first');
cai_30_150(i)=(b)/10 %ms
% %%
% [a,b]=max(vlsm);
% [c,d]=max(vrsm);
% v_V_C_50(i)=50/(d-b) %mm/ms
% 
% [a,b]=max(clsm);
% [c,d]=max(crsm);
% cai_V_C_50(i)=50/(d-b) %mm/ms
% 
% [a,b]=max(vlla);
% [c,d]=max(vrla);
% v_V_C_150(i)=50/(d-b) %mm/ms
% 
% [a,b]=max(clla);
% [c,d]=max(crla);
% cai_V_C_150(i)=50/(d-b) %mm/ms

end
%%
i=9;

rtsm=load('recorded_times_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
vlsm=load('recorded_v_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
clsm_left=load('recorded_cai_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');

rtla=load('recorded_times_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
vlla=load('recorded_v_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');
clla_left=load('recorded_cai_single_cell__1.0_1.0_1.0_1.0_1.0_1.0_1.0.txt');

% 
% var = 'g_Na';
% 
% rtsm=load(char(strcat('recorded_times_', var, '_1.0.txt')));
% vlsm=load(char(strcat('v_left_', var, '_1.0.txt')));
% vrsm=load(char(strcat('v_right_', var, '_1.0.txt')));
% clsm=load(char(strcat('cai_left_', var, '_1.0.txt')));
% crsm=load(char(strcat('cai_right_', var, '_1.0.txt')));

rtsm=rtsm(1:10000);
vlsm=vlsm(1:10000);
clsm=clsm(1:10000);

% rtla=load(char(strcat('recorded_times_', var, '_1.0.txt')));
% vlla=load(char(strcat('v_left_', var, '_1.0.txt')));
% vrla=load(char(strcat('v_right_', var, '_1.0.txt')));
% clla=load(char(strcat('cai_left_', var, '_1.0.txt')));
% crla=load(char(strcat('cai_right_', var, '_1.0.txt')));

rtla=rtla(1:10000);
vlla=vlla(1:10000);
clla=clla(1:10000);

%%
v_ampli_50(i)=max(vlsm)-min(vlsm) %V
cai_ampli_50(i)=max(clsm)-min(clsm) %mM

v_ampli_150(i)=max(vlla)-min(vlla) %V
cai_ampli_150(i)=max(clla)-min(clla) %mM
%%
[a,b]=max(vlsm);
v_V_max_50(i)=max(vlsm(2:b)-vlsm(1:b-1))*10000 %V/s

[a,b]=max(clsm);
cai_V_max_50(i)=max(clsm(2:b)-clsm(1:b-1))*10000 %mM/s

[a,b]=max(vlla);
v_V_max_150(i)=max(vlla(2:b)-vlla(1:b-1))*10000 %V/s

[a,b]=max(clla);
cai_V_max_150(i)=max(clla(2:b)-clla(1:b-1))*10000 %mM/s
%%
a=find(vlsm>(min(vlsm)+0.1*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.1*v_ampli_50(i)), 1, 'first');
v_90_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.1*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.1*cai_ampli_50(i)), 1, 'first');
cai_90_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.1*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.1*v_ampli_150(i)), 1, 'first');
v_90_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.1*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.1*cai_ampli_150(i)), 1, 'first');
cai_90_150(i)=(b)/10 %ms
%%
a=find(vlsm>(min(vlsm)+0.3*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.3*v_ampli_50(i)), 1, 'first');
v_70_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.3*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.3*cai_ampli_50(i)), 1, 'first');
cai_70_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.3*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.3*v_ampli_150(i)), 1, 'first');
v_70_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.3*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.3*cai_ampli_150(i)), 1, 'first');
cai_70_150(i)=(b)/10 %ms
%%
a=find(vlsm>(min(vlsm)+0.5*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.5*v_ampli_50(i)), 1, 'first');
v_50_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.5*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.5*cai_ampli_50(i)), 1, 'first');
cai_50_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.5*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.5*v_ampli_150(i)), 1, 'first');
v_50_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.5*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.5*cai_ampli_150(i)), 1, 'first');
cai_50_150(i)=(b)/10 %ms
%%
a=find(vlsm>(min(vlsm)+0.7*v_ampli_50(i)), 1, 'first');
b=find(vlsm(a:end)<(min(vlsm)+0.7*v_ampli_50(i)), 1, 'first');
v_30_50(i)=(b)/10 %ms

a=find(clsm>(min(clsm)+0.7*cai_ampli_50(i)), 1, 'first');
b=find(clsm(a:end)<(min(clsm)+0.7*cai_ampli_50(i)), 1, 'first');
cai_30_50(i)=(b)/10 %ms

a=find(vlla>(min(vlla)+0.7*v_ampli_150(i)), 1, 'first');
b=find(vlla(a:end)<(min(vlla)+0.7*v_ampli_150(i)), 1, 'first');
v_30_150(i)=(b)/10 %ms

a=find(clla>(min(clla)+0.7*cai_ampli_150(i)), 1, 'first');
b=find(clla(a:end)<(min(clla)+0.7*cai_ampli_150(i)), 1, 'first');
cai_30_150(i)=(b)/10 %ms
%%
% [a,b]=max(vlsm);
% [c,d]=max(vrsm);
% v_V_C_50(i)=50/(d-b) %mm/ms
% 
% [a,b]=max(clsm);
% [c,d]=max(crsm);
% cai_V_C_50(i)=50/(d-b) %mm/ms
% 
% [a,b]=max(vlla);
% [c,d]=max(vrla);
% v_V_C_150(i)=50/(d-b) %mm/ms
% 
% [a,b]=max(clla);
% [c,d]=max(crla);
% cai_V_C_150(i)=50/(d-b) %mm/ms

%% Plot amplitude
close all
FigHandle = figure('Position', [206   475   900   223]);

subplot(1,2,2)
y=1000*[v_ampli_50(1:7);v_ampli_50(9)*ones(1,7);v_ampli_150(1:7)]';
a=bar(y);
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ax = gca;
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

subplot(1,2,1)
y=[cai_ampli_50(1:7);cai_ampli_50(9)*ones(1,7);cai_ampli_150(1:7)]';
a=bar(y)
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
ax = gca;
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

%% Plot Vmax
close all
FigHandle = figure('Position', [206   475   900   223]);

subplot(1,2,2)
y=1000*[v_V_max_50(1:7);v_V_max_50(9)*ones(1,7);v_V_max_150(1:7)]';
a=bar(y);
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ax = gca;
ylabel('(mV/s)','FontSize',20,'Interpreter','Latex');
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

subplot(1,2,1)
y=[cai_V_max_50(1:7);cai_V_max_50(9)*ones(1,7);cai_V_max_150(1:7)]';
a=bar(y)
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ylabel('(mM/s)','FontSize',20,'Interpreter','Latex')
ax = gca;
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

% %% Plot Vc
% close all
% FigHandle = figure('Position', [206   475   900   223]);
% 
% subplot(1,2,2)
% y=[v_V_C_50(1:7);v_V_C_50(9)*ones(1,7);v_V_C_150(1:7)]';
% a=bar(y);
% a(1).FaceColor = [127,205,187]./255;
% a(2).FaceColor = [29,145,192]./255;
% a(3).FaceColor = [8,29,88]./255;
% ax = gca;
% ylabel('(m/s)','FontSize',20,'Interpreter','Latex');
% ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}','\sigma_t'});
% set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
% xlim([0 8])
% 
% subplot(1,2,1)
% y=[cai_V_C_50(1:7);cai_V_C_50(9)*ones(1,7);cai_V_C_150(1:7)]';
% a=bar(y)
% a(1).FaceColor = [127,205,187]./255;
% a(2).FaceColor = [29,145,192]./255;
% a(3).FaceColor = [8,29,88]./255;
% ylabel('(m/s)','FontSize',20,'Interpreter','Latex')
% ax = gca;
% ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}','\sigma_t'});
% set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
% xlim([0 8])

%% Plot 90%
close all
FigHandle = figure('Position', [206   475   900   223]);

subplot(1,2,2)
y=[v_90_50(1:7);v_90_50(9)*ones(1,7);v_90_150(1:7)]';
a=bar(y);
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ax = gca;
ylabel('(ms)','FontSize',20,'Interpreter','Latex');
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

subplot(1,2,1)
y=[cai_90_50(1:7);cai_90_50(9)*ones(1,7);cai_90_150(1:7)]';
a=bar(y)
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ylabel('(ms)','FontSize',20,'Interpreter','Latex')
ax = gca;
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

%% Plot 70%
close all
FigHandle = figure('Position', [206   475   900   223]);

subplot(1,2,2)
y=[v_70_50(1:7);v_70_50(9)*ones(1,7);v_70_150(1:7)]';
a=bar(y);
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ax = gca;
ylabel('(ms)','FontSize',20,'Interpreter','Latex');
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

subplot(1,2,1)
y=[cai_70_50(1:7);cai_70_50(9)*ones(1,7);cai_70_150(1:7)]';
a=bar(y)
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ylabel('(ms)','FontSize',20,'Interpreter','Latex')
ax = gca;
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

%% Plot 50%
close all
FigHandle = figure('Position', [206   475   900   223]);

subplot(1,2,2)
y=[v_50_50(1:7);v_50_50(9)*ones(1,7);v_50_150(1:7)]';
a=bar(y);
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ax = gca;
ylabel('(ms)','FontSize',20,'Interpreter','Latex');
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

subplot(1,2,1)
y=[cai_50_50(1:7);cai_50_50(9)*ones(1,7);cai_50_150(1:7)]';
a=bar(y)
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ylabel('(ms)','FontSize',20,'Interpreter','Latex')
ax = gca;
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

%% Plot 30%
close all
FigHandle = figure('Position', [206   475   900   223]);

subplot(1,2,2)
y=[v_30_50(1:7);v_30_50(9)*ones(1,7);v_30_150(1:7)]';
a=bar(y);
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ax = gca;
ylabel('(ms)','FontSize',20,'Interpreter','Latex');
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

subplot(1,2,1)
y=[cai_30_50(1:7);cai_30_50(9)*ones(1,7);cai_30_150(1:7)]';
a=bar(y)
a(1).FaceColor = [127,205,187]./255;
a(2).FaceColor = [29,145,192]./255;
a(3).FaceColor = [8,29,88]./255;
ylabel('(ms)','FontSize',20,'Interpreter','Latex')
ax = gca;
ax.XTickLabel = ({'g_{Na}','g_{CaL}','g_{Kr}','g_{K1}','g_{Ks}','g_{f}','g_{to}'});
set(ax,'FontSize',15,'TickLabelInterpreter', 'tex');
xlim([0 8])

