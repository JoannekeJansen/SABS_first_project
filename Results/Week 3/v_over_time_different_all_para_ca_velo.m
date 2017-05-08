close all
clear all

load('recorded_times.txt')
recorded_times2=(0:1:500);

%% g_Cal
figure('pos',[10 10 1000 800])

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_CaL_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_CaL_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_CaL_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_CaL_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$g\_CaL$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$g\_CaL$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$g\_CaL$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$g\_CaL$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_K1

figure('pos',[10 10 1000 800])

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_K1_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_K1_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_K1_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_K1_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$g\_K1$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$g\_K1$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$g\_K1$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$g\_K1$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_Ks

figure('pos',[10 10 1000 800]) 

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_Ks_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_Ks_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_Ks_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_Ks_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$g\_Ks$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$g\_Ks$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$g\_Ks$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$g\_Ks$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_Kr

figure('pos',[10 10 1000 800]) 

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_Kr_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_Kr_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_Kr_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_Kr_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$g\_Kr$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$g\_Kr$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$g\_Kr$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$g\_Kr$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_Na
figure('pos',[10 10 1000 800]) 

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_Na_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_Na_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_Na_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_Na_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$g\_Na$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$g\_Na$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$g\_Na$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$g\_Na$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% g_to
figure('pos',[10 10 1000 800]) 

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_g_to_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_g_to_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_g_to_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_g_to_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$g\_to$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$g\_to$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$g\_to$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$g\_to$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');


%% SR_Ca_release_ks
figure('pos',[10 10 1000 800]) 

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_SR_Ca_release_ks_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$ks$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$ks$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$ks$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$ks$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');


%% Sigma_l
figure('pos',[10 10 1000 800]) 

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_sigma_l_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_sigma_l_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_sigma_l_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_sigma_l_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$\sigma_l$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$\sigma_l$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$\sigma_l$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

h1=suptitle('$\sigma_l$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');

%% Sigma_t
figure('pos',[10 10 1000 800]) 

factor=(0.1:0.1:2.0);
velo_l=zeros(size(factor));
velo_r=zeros(size(factor));
velo_l2=zeros(size(factor));
velo_r2=zeros(size(factor));
caimax=zeros(size(factor));
subplot(2,3,1)
hold on

for i = (1:1:length(factor))
    a = char(strcat('cai_left_sigma_t_', num2str(factor(i)), '.txt'));
    b = char(strcat('cai_right_sigma_t_', num2str(factor(i)), '.txt'));
    A=load(a);
    B=load(b);
    plot(recorded_times2,A, 'Linewidth', 3)
    plot(recorded_times2,B, 'Linewidth', 3)
    [camaxv,idxa] = max(A);
    [~, idxb] = max(B);
    caimax(i)=5/(idxb-idxa);
end
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

subplot(2,3,2)
hold on

for i = (1:1:length(factor))
    c = char(strcat('v_left_sigma_t_', num2str(factor(i)), '.txt'));
    d = char(strcat('v_right_sigma_t_', num2str(factor(i)), '.txt'));
    C=load(c);
    D=load(d);
    plot(recorded_times2,C, 'Linewidth', 3)
    plot(recorded_times2,D, 'Linewidth', 3)
    velo_l(i)=find(C>0,1,'first');
    velo_r(i)=find(D>0,1,'first');
    if (D(end) < -30)
    velo_l2(i)=velo_l(i)+find(C(velo_l(i)+1:end)<-30,1,'first');
    velo_r2(i)=velo_r(i)+find(-D(velo_r(i)+1:end)>30,1,'first');
    end
end
legend_str = {'$10\%$ (left)','$10\%$ (right)', ...
    '$20\%$ (left)','$20\%$ (right)', ...
    '$30\%$ (left)','$30\%$ (right)', ...
    '$40\%$ (left)','$40\%$ (right)', ...
    '$50\%$ (left)','$50\%$ (right)', ...
    '$60\%$ (left)','$60\%$ (right)', ...
    '$70\%$ (left)','$70\%$ (right)', ...
    '$80\%$ (left)','$80\%$ (right)', ...
    '$90\%$ (left)','$90\%$ (right)', ...
    '$100\%$ (left)','$100\%$ (right)', ...
    '$110\%$ (left)','$110\%$ (right)', ...
    '$120\%$ (left)','$120\%$ (right)', ...
    '$130\%$ (left)','$130\%$ (right)', ...
    '$140\%$ (left)','$140\%$ (right)', ...
    '$150\%$ (left)','$150\%$ (right)', ...
    '$160\%$ (left)','$160\%$ (right)', ...
    '$170\%$ (left)','$170\%$ (right)', ...
    '$180\%$ (left)','$180\%$ (right)', ...
    '$190\%$ (left)','$190\%$ (right)', ...
    '$200\%$ (left)','$200\%$ (right)'};
[l,h]=columnlegend(3,legend_str,'location','east');
set(findall(h,'type', 'text'), 'interpreter', 'Latex');
fig_pos = get(gca, 'position');
pos = get(l, 'position');
set(l, 'position', [pos(1)+1.1*fig_pos(3) pos(2)-0.5*fig_pos(4) pos(3) pos(4)]);

xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex');
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');

subplot(2,3,4)
hold on
velo = 5./(velo_r-velo_l); %/ms
plot(factor,velo, 'Linewidth', 3)
xlabel('$\sigma_t$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity up (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,5)
hold on
velo2 = 5./(velo_r2-velo_l2); %/ms
plot(factor,velo2, 'Linewidth', 3)
xlabel('$\sigma_t$ factor','FontSize',20,'Interpreter','Latex');
ylabel('Conduction velocity down (mm/ms)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0.0 0.2])

subplot(2,3,6)
hold on
plot(factor,caimax, 'Linewidth', 3)
xlabel('$\sigma_t$ factor','FontSize',20,'Interpreter','Latex');
ylabel('$[Ca]_i$ peak velocity (mM)','FontSize',20,'Interpreter','Latex');
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex');
ylim([0 0.2])

h1=suptitle('$\sigma_t$ factor');
set(h1,'FontSize',30,'Interpreter', 'Latex');
%%
shg
