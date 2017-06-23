close all
clear all

A=load('J2_values_1.0_1.0_1.0_1.0_12mm_strip_gna.txt');
B=load('J2_values_1.0_1.0_1.0_1.0_12mm_strip_gcal.txt');
C=load('J2_values_1.0_1.0_1.0_1.0_12mm_strip_gkr.txt');
D=load('J2_values_1.0_1.0_1.0_1.0_12mm_strip_gk12.txt');
vars={'(total)','(left)','(right)','\mbox{wave speed}','\mbox{amplitude $v$}','\mbox{amplitude $[Ca]_i$}', '\mbox{upstroke velocity $[Ca]_i$}', '\mbox{upstroke velocity $v$}', 'v \:30 \%', ' [Ca]_i \:30 \%', 'v \:50 \%', ' [Ca]_i \:50_p', 'v \:70 \%', ' [Ca]_i \:70 \%', 'v\: 90 \%', ' [Ca]_i\: 90 \%'};
for i=4:16
FigHandle = figure('Position', [198   448   696   257]);
colorord = get(gca, 'ColorOrder');
hold on
plot(A(:,17),A(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor', colorord(1,:))
plot(B(:,17),B(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(2,:))
plot(C(:,17),C(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(3,:))
plot(D(:,17),D(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(4,:))

xlim([0.75 1.25])
a=char(strcat('$\mathcal{J}_{', vars(i),'}$'));
xlabel('Percentage of default value','FontSize',20,'Interpreter','Latex')
ylabel(a, 'FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25])
set(gca, 'XTickLabel', {'$0.75$', '$0.8$', '$0.85$', '$0.9$', '$0.95$', '$1.0$', '$1.05$', '$1.1$', '$1.15$', '$1.2$', '$1.25$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')

l = legend('$g_{Na}$','$g_{CaL}$','$g_{Kr}$','$g_{K1}$');
l.Orientation = 'horizontal';
l.Location = 'North';
l.Interpreter = 'Latex';
l.FontSize = 15;
end

for i=1:3
FigHandle = figure('Position', [198   448   696   257]);
colorord = get(gca, 'ColorOrder');
hold on
plot(A(:,17),A(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor', colorord(1,:))
plot(B(:,17),B(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(2,:))
plot(C(:,17),C(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(3,:))
plot(D(:,17),D(:,i),'-o', 'MarkerSize', 7, 'MarkerFaceColor',colorord(4,:))

xlim([0.75 1.25])
a=char(strcat('$\mathcal{J}$ ', vars(i),''));
xlabel('Percentage of default value','FontSize',20,'Interpreter','Latex')
ylabel(a, 'FontSize',20,'Interpreter','Latex')
set(gca, 'XTick', [0.75 0.8 0.85 0.9 0.95 1.0 1.05 1.1 1.15 1.2 1.25])
set(gca, 'XTickLabel', {'$0.75$', '$0.8$', '$0.85$', '$0.9$', '$0.95$', '$1.0$', '$1.05$', '$1.1$', '$1.15$', '$1.2$', '$1.25$'})
set(gca,'FontSize',15,'TickLabelInterpreter', 'Latex')

l = legend('$g_{Na}$','$g_{CaL}$','$g_{Kr}$','$g_{K1}$');
l.Orientation = 'horizontal';
l.Location = 'North';
l.Interpreter = 'Latex';
l.FontSize = 15;
end
shg