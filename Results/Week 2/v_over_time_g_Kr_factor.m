close all
clear all

load('values_g_Kr_factor_0.1.txt')
load('values_g_Kr_factor_0.2.txt')
load('values_g_Kr_factor_0.3.txt')
load('values_g_Kr_factor_0.4.txt')
load('values_g_Kr_factor_0.5.txt')
load('values_g_Kr_factor_0.6.txt')
load('values_g_Kr_factor_0.7.txt')
load('values_g_Kr_factor_0.8.txt')
load('values_g_Kr_factor_0.9.txt')
load('values_g_Kr_factor_1.txt')
load('values_g_Kr_factor_0.txt')

subplot(1,2,1)
hold on
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_1(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_9(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_8(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_7(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_6(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_5(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_4(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_3(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_2(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_1(1:end-1,3), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0(1:end-1,3), 'Linewidth', 3)

xlim([1 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$[Ca]_i$ (mM)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
subplot(1,2,2)
hold on
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_1(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_9(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_8(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_7(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_6(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_5(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_4(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_3(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_2(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0_1(1:end-1,2), 'Linewidth', 3)
plot(values_g_Kr_factor_1(1:end-1,1),values_g_Kr_factor_0(1:end-1,2), 'Linewidth', 3)

xlim([1 499])
xlabel('$t$ (ms)','FontSize',20,'Interpreter','Latex')
ylabel('$v$ (mV)','FontSize',20,'Interpreter','Latex')
set(gca,'FontSize',15,'TickLabelInterpreter', 'tex')
h=legend('$100\%$','$90\%$','$80\%$','$70\%$','$60\%$','$50\%$','$40\%$','$30\%$','$20\%$','$10\%$','$0\%$');
set(h,'FontSize',10,'Interpreter', 'Latex')
h1=suptitle('$g\_Kr$');
set(h1,'FontSize',30,'Interpreter', 'Latex')
shg