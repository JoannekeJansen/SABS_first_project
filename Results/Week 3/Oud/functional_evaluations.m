close all
clear all

A1=[0.10023123904; 0.0605024672409; 0.03515795203968196531; 0.0178271723453; 5.868856901723199133e-03; 3.006519911435719518e-06];
B3= (0.5 : 0.1 : 1.0);

figure
plot(B3,A1,'*', 'Linewidth', 4)
hold on
idx = find(min(abs(B3-1))==abs(B3-1));
plot(B3(idx),A1(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%h1=title('Minidomain');
%set(h1,'FontSize',25,'Interpreter', 'Latex')
xlabel('$g\_Na$ factor','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_Na)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
shg

% 0.10023123904 0.5
% 0.0605024672409 0.6
% 0.03515795203968196531 0.7
% 0.0178271723453 0.8
% 5.868856901723199133e-03 0.9
% 3.006519911435719518e-06 1.0

