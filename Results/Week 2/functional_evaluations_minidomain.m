close all
clear all
% g_K1
% Control values and J value 0.8 0.0243219721615
% Control values and J value 0.9 0.00631999426715
% Control values and J value 1.0 8.61952650229e-06
% Control values and J value 1.1 0.00823217488743
% Control values and J value 1.2 0.0343262662775
% 
% g_CaL
% Control values and J value 0.8 155.839512253
% Control values and J value 0.9 36.1586291941
% Control values and J value 1.0 8.61952650229e-06
% Control values and J value 1.1 29.9848434951
% Control values and J value 1.2 108.734512786
% 
% g_Kr
% Control values and J value 0.8 1.55027073199
% Control values and J value 0.9 0.394283404582
% Control values and J value 1.0 8.61952650229e-06
% Control values and J value 1.1 0.403457678855
% Control values and J value 1.2 1.64252579319

A1=[0.0243219721615; 0.00631999426715; 8.61952650229e-06; 0.00823217488743; 0.0343262662775];
A2=[155.839512253; 36.1586291941; 8.61952650229e-06; 29.9848434951; 108.734512786];
A3=[1.55027073199; 0.394283404582; 8.61952650229e-06; 0.403457678855; 1.64252579319];
B3= (0.8 : 0.1 : 1.2);

figure
plot(B3,A1,'*', 'Linewidth', 4)
hold on
%xlim([0 1.5])
idx = find(min(abs(B3-1))==abs(B3-1));
plot(B3(idx),A2(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%set(gca, 'XTick', linspace(0, 1.5, 16))
h1=title('Minidomain');
set(h1,'FontSize',25,'Interpreter', 'Latex')
xlabel('$g\_K1$ factor','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_K1)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
figure
plot(B3,A2,'*', 'Linewidth', 4)
hold on
%xlim([0 1.5])
idx = find(min(abs(B3-1))==abs(B3-1));
plot(B3(idx),A2(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%set(gca, 'XTick', linspace(0, 1.5, 16))
h1=title('Minidomain');
set(h1,'FontSize',25,'Interpreter', 'Latex')
xlabel('$g\_Cal$ factor','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_CaL)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
figure
plot(B3,A3,'*', 'Linewidth', 4)
hold on
%xlim([0 1.5])
idx = find(min(abs(B3-1))==abs(B3-1));
plot(B3(idx),A3(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%set(gca, 'XTick', linspace(0, 1.5, 16))
h1=title('Minidomain');
set(h1,'FontSize',25,'Interpreter', 'Latex')
xlabel('$g\_Kr$ factor','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_Kr)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
shg