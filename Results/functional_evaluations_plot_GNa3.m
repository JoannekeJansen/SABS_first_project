close all
clear all

load('J_results_gna.txt')
A=J_results_gna(:,1);
B=J_results_gna(:,2);
plot(B,A,'*', 'Linewidth', 4)
hold on
%xlim([0 35])
idx = find(min(abs(B-23))==abs(B-23));
plot(B(idx),A(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%set(gca, 'XTick', [0:2:35])
xlabel('$g\_Na$','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_Na)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
figure
load('J_results_sigma.txt')
A=J_results_sigma(:,1);
B=J_results_sigma(:,2);
plot(B,A,'*', 'Linewidth', 4)
hold on
%xlim([0 35])
idx = find(min(abs(B-0.15))==abs(B-0.15));
plot(B(idx),A(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
set(gca, 'XTick', [0:2:35])
xlabel('$g\_Na$','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_Na)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')

shg