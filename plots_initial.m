clc
clear 
close all
format long e
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextInterpreter', 'tex'); 
set(0,'DefaultTextFontSize',18);
addpath('./ORIGIN/')

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\1_5000_5comp_withoutNO_ex3.mat')

n_N2_1_5 = n_N2;
n_O2_1_5 = n_O2;
n_NO_1_5 = n_NO;
n_N_1_5 = n_N;
n_O_1_5 = n_O;
T_1_5 = T;
v_1_5 = v;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_5000_5comp_withoutNO_ex3.mat')

n_N2_100_5 = n_N2;
n_O2_100_5 = n_O2;
n_NO_100_5 = n_NO;
n_N_100_5 = n_N;
n_O_100_5 = n_O;
T_100_5 = T;
v_100_5 = v;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\1_7000_5comp_withoutNO_ex3.mat')

n_N2_1_7 = n_N2;
n_O2_1_7 = n_O2;
n_NO_1_7 = n_NO;
n_N_1_7 = n_N;
n_O_1_7 = n_O;
T_1_7 = T;
v_1_7 = v;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex3.mat')

n_N2_100_7 = n_N2;
n_O2_100_7 = n_O2;
n_NO_100_7 = n_NO;
n_N_100_7 = n_N;
n_O_100_7 = n_O;
T_100_7 = T;
v_100_7 = v;

figure1 = figure(1);
plot(X,T_1_7,'r--','linewidth',2), hold on
plot(X,T_100_7,'r-','linewidth',2)
plot(X,T_1_5,'b--','linewidth',2)
plot(X,T_100_5,'b-','linewidth',2)
ylabel('T, K');
xlabel('x/r^*');
legend({'T^* = 7000 K, p^* = 1 atm','T^* = 7000 K, p^* = 100 atm',...
        'T^* = 5000 K, p^* = 1 atm','T^* = 5000 K, p^* = 100 atm'},'location','northeast');
set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure1,'./ORIGIN/PLOTS/T_initial','-dpdf','-r0');
hold off

figure2 = figure(2);
plot(X,v_1_7,'r--','linewidth',2), hold on
plot(X,v_100_7,'r-','linewidth',2)
plot(X,v_1_5,'b--','linewidth',2)
plot(X,v_100_5,'b-','linewidth',2)
ylabel('v, m/s');
xlabel('x/r^*');
legend({'T^* = 7000 K, p^* = 1 atm','T^* = 7000 K, p^* = 100 atm',...
        'T^* = 5000 K, p^* = 1 atm','T^* = 5000 K, p^* = 100 atm'},'location','southeast');
set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure2,'./ORIGIN/PLOTS/v_initial','-dpdf','-r0');
hold off

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\10_7000_5comp_withoutNO_ex3.mat')

n_N2_10_7 = n_N2;
n_O2_10_7 = n_O2;
n_NO_10_7 = n_NO;
n_N_10_7 = n_N;
n_O_10_7 = n_O;
T_10_7 = T;
v_10_7 = v;

figure3 = figure(3);
semilogy(X, n_N2_1_7, 'k', X, n_N2_10_7, 'b', X, n_N2_100_7, 'r','linewidth',2)
title('n_{N_2}/n')
ylabel('n_{N_2}/n')
xlabel('x/r^*')
legend({'1 atm','10 atm','100 atm'})
xlim([0,5]);
set(figure3,'Units','Inches');
pos = get(figure3,'Position');
set(figure3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure3,'./ORIGIN/PLOTS/N2_initial','-dpdf','-r0');

figure4 = figure(4);
semilogy(X, n_O2_1_7, 'k', X, n_O2_10_7, 'b', X, n_O2_100_7, 'r','linewidth',2)
title('n_{O_2}/n')
legend({'1 atm','10 atm','100 atm'})
ylabel('n_{O_2}/n')
xlim([0,5]);
xlabel('x/r^*')
set(figure4,'Units','Inches');
pos = get(figure4,'Position');
set(figure4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure4,'./ORIGIN/PLOTS/O2_initial','-dpdf','-r0');


load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\10_7000_5comp_withoutNO_ex1.mat')

n_N2_10_7 = n_N2;
n_O2_10_7 = n_O2;
n_NO_10_7 = n_NO;
n_N_10_7 = n_N;
n_O_10_7 = n_O;
T_10_7 = T;
v_10_7 = v;


load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex1.mat')

n_N2_100_7 = n_N2;
n_O2_100_7 = n_O2;
n_NO_100_7 = n_NO;
n_N_100_7 = n_N;
n_O_100_7 = n_O;
T_100_7 = T;
v_100_7 = v;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\1_7000_5comp_withoutNO_ex1.mat')

n_N2_1_7 = n_N2;
n_O2_1_7 = n_O2;
n_NO_1_7 = n_NO;
n_N_1_7 = n_N;
n_O_1_7 = n_O;
T_1_7 = T;
v_1_7 = v;

figure5 = figure(5);
semilogy(X, n_N2_1_7, 'k', X, n_N2_10_7, 'b', X, n_N2_100_7, 'r','linewidth',2)
title('n_{N_2}/n')
ylabel('n_{N_2}/n')
xlim([0,5]);
xlabel('x/r^*')
legend({'1 atm','10 atm','100 atm'})
set(figure5,'Units','Inches');
pos = get(figure5,'Position');
set(figure5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure5,'./ORIGIN/PLOTS/N2_initial_ex1','-dpdf','-r0');

figure6 = figure(6);
semilogy(X, n_O2_1_7, 'k', X, n_O2_10_7, 'b', X, n_O2_100_7, 'r','linewidth',2)
title('n_{O_2}/n')
legend({'1 atm','10 atm','100 atm'})
ylabel('n_{O_2}/n')
xlim([0,5]);
xlabel('x/r^*')
set(figure6,'Units','Inches');
pos = get(figure6,'Position');
set(figure6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure6,'./ORIGIN/PLOTS/O2_initial_ex1','-dpdf','-r0');

