close all
clear all

A=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gna_gcal.txt');
A1=reshape(A(1:2601,1),51,51);
A2=reshape(A(1:2601,2),51,51);
A3=reshape(A(1:2601,3),51,51);

FigHandle = figure('Position', [1 496 1027 209]);
subplot(1,3,1)
image(A1,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,2)
image(A2,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,3)
image(A3,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
%%
A=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gna_gk1.txt');
A1=reshape(A(1:2601,1),51,51);
A2=reshape(A(1:2601,2),51,51);
A3=reshape(A(1:2601,3),51,51);

FigHandle = figure('Position', [1 496 1027 209]);
subplot(1,3,1)
image(A1,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,2)
image(A2,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,3)
image(A3,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
%%
A=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gna_gkr.txt');
A1=reshape(A(1:2601,1),51,51);
A2=reshape(A(1:2601,2),51,51);
A3=reshape(A(1:2601,3),51,51);

FigHandle = figure('Position', [1 496 1027 209]);
subplot(1,3,1)
image(A1,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,2)
image(A2,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,3)
image(A3,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Na}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
%%
A=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gkr_gcal.txt');
A1=reshape(A(1:2601,1),51,51);
A2=reshape(A(1:2601,2),51,51);
A3=reshape(A(1:2601,3),51,51);

FigHandle = figure('Position', [1 496 1027 209]);
subplot(1,3,1)
image(A1,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,2)
image(A2,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,3)
image(A3,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
%%
A=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gk1_gcal.txt');
A1=reshape(A(1:2601,1),51,51);
A2=reshape(A(1:2601,2),51,51);
A3=reshape(A(1:2601,3),51,51);

FigHandle = figure('Position', [1 496 1027 209]);
subplot(1,3,1)
image(A1,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,2)
image(A2,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,3)
image(A3,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{CaL}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
%%
A=load('J_values_1.0_1.0_1.0_1.0_12mm_strip_gk1_gkr.txt');
A1=reshape(A(1:2601,1),51,51);
A2=reshape(A(1:2601,2),51,51);
A3=reshape(A(1:2601,3),51,51);

FigHandle = figure('Position', [1 496 1027 209]);
subplot(1,3,1)
image(A1,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,2)
image(A2,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
subplot(1,3,3)
image(A3,'CDataMapping','scaled')
set(gca,'YDir','normal')
xlabel('$g_{K1}$','FontSize',20,'Interpreter','Latex')
ylabel('$g_{Kr}$','FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [1 26 51])
set(gca, 'XTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca, 'YTick', [1 26 51])
set(gca, 'YTickLabel', {'$75\%$','$100\%$', '$125\%$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')
colorbar
%%
shg