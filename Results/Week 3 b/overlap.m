close all
clear all

%load('recorded_times.txt')
factor=(0.1:0.2:1.9);

recorded_times2=(0:0.01:500);
incr=zeros(size(factor,2),9);

%% g_Cal
figure('pos',[10 10 1200 500])

subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_CaL_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_CaL_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,1)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_CaL_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_CaL_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('$g\_CaL$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_K1

figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_K1_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_K1_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,2)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_K1_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_K1_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('$g\_K1$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_Ks

figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_Ks_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_Ks_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,3)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_Ks_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_Ks_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('$g\_Ks$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_Kr

figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_Kr_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_Kr_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,4)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_Kr_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_Kr_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('$g\_Kr$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_Na
figure('pos',[10 10 1000 800]) 

figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_Na_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_Na_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,5)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_Na_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_Na_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('$g\_Na$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_to
figure('pos',[10 10 1000 800]) 

figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_to_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_to_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,6)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_to_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_to_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');


h1=suptitle('$g\_to$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');


%% SR_Ca_release_ks 
figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,7)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('$ks$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');


%% Sigma_l
figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_sigma_l_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_sigma_l_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,8)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_sigma_l_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_sigma_l_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

h1=suptitle('$\sigma_l$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% Sigma_t
figure('pos',[10 10 1200 500])

factor=(0.1:0.2:1.9);
subplot(1,3,1)
hold on

velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
for i = (1:1:length(factor))
    a = char(strcat('cai_left_sigma_t_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_sigma_t_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    [incra,idxa] = max(A);
    [incrb, idxb] = max(B);
    incr(i,9)=100*(incrb-incra)/incra;
    plot(recorded_times2(idxb-idxa+1:end),A(1:idxa+(length(A)-idxb)), 'Linewidth', 3)
    plot(recorded_times2(idxb-idxa+1:end),B(idxb-idxa+1:end),':', 'Linewidth', 3)
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(1,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_sigma_t_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_sigma_t_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),C(1:velo_l(i)+(length(C)-velo_r(i))), 'Linewidth', 3)
    plot(recorded_times2(velo_r(i)-velo_l(i)+1:end),D(velo_r(i)-velo_l(i)+1:end),':', 'Linewidth', 3)
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)'};
[l,h]=columnlegend(1,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.1*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');


h1=suptitle('$\sigma_t$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');
%%

figure
hold on
plot(factor,incr(:,1))
plot(factor,incr(:,2))
plot(factor,incr(:,3))
plot(factor,incr(:,4))
plot(factor,incr(:,5))
plot(factor,incr(:,6))
plot(factor,incr(:,7))
plot(factor,incr(:,8))
plot(factor,incr(:,9))
legend('g_cal','g_K1','g_Ks','g_Kr','g_Na','g_to','ks','sigma_l','sigma_t')
shg
