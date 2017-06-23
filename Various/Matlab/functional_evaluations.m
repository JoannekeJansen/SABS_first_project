close all
clear all
FigHandle = figure('Position', [206    70   714   900]);

A=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gna.txt');
B=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gcal.txt');
C=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gkr.txt');
D=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gk1.txt');
subplot(2,1,2)
colorord = get(gca, 'ColorOrder');
%semilogy(A(:,4),A(:,2),'*', 'Linewidth', 4)
plot(A(:,4),A(:,2),'-o', 'MarkerSize', 7, 'MarkerFaceColor', colorord(1,:))
hold on
plot(B(:,4),B(:,2),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(2,:))
plot(C(:,4),C(:,2),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(3,:))
plot(D(:,4),D(:,2),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(4,:))
%plot(E(7,4),E(7,2),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%ylim([0 0.5])
xlim([0.75 1.25])
xlabel('Percentage of default value','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}$ (right side)','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [0.5 0.6 0.7 0.8 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.4 1.5])
%set(gca, 'XTickLabel', {'$50\%$','$60\%$','$70\%$','$80\%$','$90\%$','$95\%$','$100\%$','$105\%$','$110\%$','$120\%$','$130\%$','$140\%$','$150\%$'})
set(gca, 'XTickLabel', {'$50$','$60$','$70$','$80$','$90$','$95$','$100$','$105$','$110$','$120$','$130$','$140$','$150$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
subplot(2,1,1)
%semilogy(A(:,4),A(:,3),'*', 'Linewidth', 4)
plot(A(:,4),A(:,3),'-o', 'MarkerSize', 7, 'MarkerFaceColor', colorord(1,:))
hold on
plot(B(:,4),B(:,3),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(2,:))
plot(C(:,4),C(:,3),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(3,:))
plot(D(:,4),D(:,3),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(4,:))
%ylim([0 0.1])
xlim([0.75 1.25])

%plot(E(7,4),E(7,3),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
xlabel('Percentage of default value','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}$ (left side)','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [0.5 0.6 0.7 0.8 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.4 1.5])
%set(gca, 'XTickLabel', {'$50\%$','$60\%$','$70\%$','$80\%$','$90\%$','$95\%$','$100\%$','$105\%$','$110\%$','$120\%$','$130\%$','$140\%$','$150\%$'})
set(gca, 'XTickLabel', {'$50$','$60$','$70$','$80$','$90$','$95$','$100$','$105$','$110$','$120$','$130$','$140$','$150$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
FigHandle = figure('Position', [302   306   714   399]);

subplot(2,1,1:2)
%semilogy(A(:,4),A(:,1),'*', 'Linewidth', 4)
plot(A(:,4),A(:,1),'-o', 'MarkerSize', 7, 'MarkerFaceColor', colorord(1,:))
hold on
plot(B(:,4),B(:,1),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(2,:))
plot(C(:,4),C(:,1),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(3,:))
plot(D(:,4),D(:,1),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(4,:))

ylim([0 0.09])
xlim([0.75 1.25])
%plot(A(7,4),A(7,1),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
xlabel('Percentage of default value','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}$ (total)','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [0.5 0.6 0.7 0.8 0.9 0.95 1.0 1.05 1.1 1.2 1.3 1.4 1.5])
%set(gca, 'XTickLabel', {'$50\%$','$60\%$','$70\%$','$80\%$','$90\%$','$95\%$','$100\%$','$105\%$','$110\%$','$120\%$','$130\%$','$140\%$','$150\%$'})
set(gca, 'XTickLabel', {'$50$','$60$','$70$','$80$','$90$','$95$','$100$','$105$','$110$','$120$','$130$','$140$','$150$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
l = legend('$g_{Na}$','$g_{CaL}$','$g_{Kr}$','$g_{K1}$','$\sigma_t$');
l.Orientation = 'horizontal';
l.Location = 'SouthOutside';
l.Interpreter = 'Latex';
l.FontSize = 15;


shg