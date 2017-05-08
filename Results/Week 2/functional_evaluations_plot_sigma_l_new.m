close all
clear all
%Compute_sigmal_J.py
%Control values and J value 0.195 0.00147188155416
%Control values and J value 0.18 0.000903548457085
%Control values and J value 0.165  0.000348335550311
%Control values and J value 0.15 3.06419440913e-05
%Control values and J value 0.135 0.000395252549906
%Control values and J value 0.12 0.00151639422458
%Control values and J value 0.105 0.00346284763428
B=[0.105, 0.12, 0.135, 0.15, 0.165, 0.18, 0.195];
A=[0.00346284763428, 0.00151639422458, 0.000395252549906, 3.06419440913e-05, 0.000348335550311, 0.000903548457085, 0.00147188155416 ];

% Control values and J value 1.3 0.549511990443
% Control values and J value 1.2 0.268457106259
% Control values and J value 1.1 0.0735259191771
% Control values and J value 1.0 6.28776737841e-05
% Control values and J value 0.9 0.0865321405431
% Control values and J value 0.8 0.355748107536
% Control values and J value 0.7 0.797752642545
B2=[1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6];
A2=[0.549511990443, 0.268457106259, 0.0735259191771, 6.28776737841e-05, 0.0865321405431, 0.355748107536, 0.797752642545, 1.35548600711
];
%single_cell_inverse
%Control values and J value 0.5 356.429544852
%Control values and J value 0.6 246.234454886
%Control values and J value 0.7 139.013550523
%Control values and J value 0.8 59.3960886327
%Control values and J value 0.9 14.4811542074
%Control values and J value 1.0 0.99635199126
%Control values and J value 1.1 12.0896942824
%Control values and J value 1.2 41.1345373958
%Control values and J value 1.3 82.9126031408
%Control values and J value 1.4 133.602686677
%Control values and J value 1.5 190.470758379
A3= [356.429544852; 246.234454886; 139.013550523;59.3960886327;14.4811542074;0.99635199126;12.0896942824;41.1345373958;82.9126031408;133.602686677;190.470758379];
B3= (0.5 : 0.1 : 1.5);

plot(B,A,'*', 'Linewidth', 4)
hold on
min(abs(B-0.15))
%xlim([0 1.5])
idx = find(min(abs(B-0.15))==abs(B-0.15));
plot(B(idx),A(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%set(gca, 'XTick', linspace(0, 1.5, 16))
xlabel('$\sigma_l$','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, \sigma_l)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
figure
plot(B2,A2,'*', 'Linewidth', 4)
hold on
%xlim([0 1.5])
idx = find(min(abs(B2-1))==abs(B2-1));
plot(B2(idx),A2(idx),'-o','markers',9,'MarkerEdgeColor','k','MarkerFaceColor',[0.6350    0.0780    0.1840])
%set(gca, 'XTick', linspace(0, 1.5, 16))
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
h1=title('Single Cell Solver');
set(h1,'FontSize',25,'Interpreter', 'Latex')
xlabel('$g\_Cal$ factor','FontSize',20,'Interpreter','Latex')
ylabel('$\mathcal{J}(v, [Ca]_i, g\_CaL)$','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
shg