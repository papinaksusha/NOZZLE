clc
clear 
close all
format long e
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextInterpreter', 'tex'); 
set(0,'DefaultTextFontSize',18);

sw_o = 2;                                                                  % switch on oscillator, 1 -  harmonic oscillator; 2 -  anharmonic oscillator 3 - STELLAR
switch sw_o 
    case 1
        wx = [0 0 0];   
    case 2
        wx = [1432 1198 1407.5];                                           % m^-1, spectroscopic constant
    case 3
        wx = [1432 1198 1407.5];   
end
                                                 
% Max vibrational levels                                                 
% I(1,:) - harmonic oscillator, I(2,:) - anharmonic oscillator, I(3,:) - STELLAR distributions

I = [33 26 27; 47 36 39; 60, 45, 47]; 


load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex5.mat')

n_i_N2_ex = n_i_N2;
n_i_O2_ex = n_i_O2;
n_N2_ex = n_N2;
n_O2_ex = n_O2;
n_NO_ex = n_NO;
n_N_ex = n_N;
n_O_ex = n_O;
T_ex = T;
v_ex = v;
X_ex = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex3_rec0.mat')

n_i_N2_rec = n_i_N2;
n_i_O2_rec = n_i_O2;
n_N2_rec = n_N2;
n_O2_rec = n_O2;
n_NO_rec = n_NO;
n_N_rec = n_N;
n_O_rec = n_O;
T_rec = T;
v_rec = v;
X_rec = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex1.mat')

n_i_N2_1 = n_i_N2;
n_i_O2_1 = n_i_O2;
n_N2_1 = n_N2;
n_O2_1 = n_O2;
n_NO_1 = n_NO;
n_N_1 = n_N;
n_O_1 = n_O;
T_1 = T;
v_1 = v;
X_1 = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex2.mat')

n_i_N2_2 = n_i_N2;
n_i_O2_2 = n_i_O2;
n_N2_2 = n_N2;
n_O2_2 = n_O2;
n_NO_2 = n_NO;
n_N_2 = n_N;
n_O_2 = n_O;
T_2 = T;
v_2 = v;
X_2 = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex3.mat')

n_i_N2_3 = n_i_N2;
n_i_O2_3 = n_i_O2;
n_N2_3 = n_N2;
n_O2_3 = n_O2;
n_NO_3 = n_NO;
n_N_3 = n_N;
n_O_3 = n_O;
T_3 = T;
v_3 = v;
X_3 = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_N2_N.mat')

X_N2 = X;
n_N2_bin = n_m;
n_N_bin = n_a;
T_N2 = T;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_O2_O.mat')

X_O2 = X;
n_O2_bin = n_m;
n_O_bin = n_a;
T_O2 = T;

figure1 = figure(1);
plot(X_1 , T_1), hold on
plot(X_2 , T_2)
plot(X_3 , T_3) 
ylim([0,7000]);
legend('Русанов','Полак','Варнатц')
hold off

figure2 = figure(2);
plot(X_3, T_3,'g-','linewidth',2), hold on
plot(X_ex , T_ex, 'b:','linewidth',2)
plot(X_rec , T_rec,'r','linewidth',2)
plot(X_N2, T_N2,'b','linewidth',2)
plot(X_O2, T_O2,'k','linewidth',2)
ylim([0,7000]);
ylabel('T, K');
xlabel('x/r^*');
legend('все реакции','без обмена','без рекомбинации','N_2/N','O_2/O');
set(figure2,'Units','Inches');
pos = get(figure2,'Position');
set(figure2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure2,'./PLOTS/T_react','-dpdf','-r0');
hold off

figure3 = figure(3);
plot(X_3, T_3,'g-','linewidth',2), hold on
plot(X_ex , T_ex, 'b:','linewidth',2)
plot(X_rec , T_rec,'r','linewidth',2)
%ylim([3000,7000]);
xlim([0,0.01]);
ylabel('T, K');
xlabel('x/r^*');
legend('все реакции','без обмена','без рекомбинации','location','southwest');
set(figure3,'Units','Inches');
pos = get(figure3,'Position');
set(figure3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
%print(figure3,'./ORIGIN/PLOTS/T_react_close','-dpdf','-r0');
hold off


%%

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex3.mat')

n_i_N2_3 = n_i_N2;
n_i_O2_3 = n_i_O2;
n_N2_3 = n_N2;
n_O2_3 = n_O2;
n_NO_3 = n_NO;
n_N_3 = n_N;
n_O_3 = n_O;
T_3 = T;
v_3 = v;
X_3 = X;
% 
load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex5.mat')

n_i_N2_ex = n_i_N2;
n_i_O2_ex = n_i_O2;
n_N2_ex = n_N2;
n_O2_ex = n_O2;
n_NO_ex = n_NO;
n_N_ex = n_N;
n_O_ex = n_O;
T_ex = T;
v_ex = v;
X_ex = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex3_rec0.mat')

n_i_N2_rec = n_i_N2;
n_i_O2_rec = n_i_O2;
n_N2_rec = n_N2;
n_O2_rec = n_O2;
n_NO_rec = n_NO;
n_N_rec = n_N;
n_O_rec = n_O;
T_rec = T;
v_rec = v;
X_rec = X;

figure4 = figure(4);
semilogy(X_3, n_N2_3,'linewidth',2), hold on
semilogy(X_3, n_O2_3,'linewidth',2)
semilogy(X_3, n_NO_3,'linewidth',2)
semilogy(X_3, n_N_3,'linewidth',2)
semilogy(X_3, n_O_3,'linewidth',2)
title('все реакции')
xlim([0,5]);
ylim([1e-6,1]);
xlabel('x/r^*');
ylabel('n_c/n');
legend('N_2','O_2','NO','N','O');
set(figure4,'Units','Inches');
pos = get(figure4,'Position');
set(figure4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure4,'all_all','-dpdf','-r0');
hold off

figure5 = figure(5);
semilogy(X_ex, n_N2_ex,'linewidth',2), hold on
semilogy(X_ex, n_O2_ex,'linewidth',2)
semilogy(X_ex, n_NO_ex,'linewidth',2)
semilogy(X_ex, n_N_ex,'linewidth',2)
semilogy(X_ex, n_O_ex,'linewidth',2)
title('без обмена')
xlim([0,5]);
ylim([1e-6,1]);
xlabel('x/r^*');
ylabel('n_c/n');
legend('N_2','O_2','NO','N','O');
set(figure5,'Units','Inches');
pos = get(figure5,'Position');
set(figure5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure5,'all_ex','-dpdf','-r0');
hold off

figure6 = figure(6);
semilogy(X_rec, n_N2_rec,'linewidth',2), hold on
semilogy(X_rec, n_O2_rec,'linewidth',2)
semilogy(X_rec, n_NO_rec,'linewidth',2)
semilogy(X_rec, n_N_rec,'linewidth',2)
semilogy(X_rec, n_O_rec,'linewidth',2)
title('без рекомбинации')
xlim([0,5]);
ylim([1e-6,1]);
xlabel('x/r^*');
ylabel('n_c/n');
legend('N_2','O_2','NO','N','O');
hold off




%%

xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

i_N2 = 0 : I(sw_o,1);
i_O2 = 0 : I(sw_o,2);

u_N2_ex = zeros(length(xr),length(i_N2));
u_O2_ex = zeros(length(xr),length(i_O2));

u_N2_3 = zeros(length(xr),length(i_N2));
u_O2_3 = zeros(length(xr),length(i_O2));

figure7 = figure(7);

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2_3(i,g) = interp1q(X_3,n_i_N2_1(:,g),xr(i));
    end
    semilogy(i_N2, u_N2_3(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{N_2i}/n');
xlim([0,I(sw_o,1)]);
ylim([1e-30,1])
xlabel('i');
title('N_2, все реакции');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',10);
set(figure7,'Units','Inches');
pos = get(figure7,'Position');
set(figure7,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure7,'N2_distr_allreact','-dpdf','-r0');
hold off

figure8 = figure(8);

for i = 1 : length(xr)
    for g = 1 : length(i_O2)
        u_O2_3(i,g) = interp1q(X_3,n_i_O2_1(:,g),xr(i));
    end
    semilogy(i_O2, u_O2_3(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{O_2i}/n');
xlim([0,I(sw_o,2)]);
title('O_2, все реакции')
%ylim([1e-15,1])
xlabel('i');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',10);
set(figure8,'Units','Inches');
pos = get(figure8,'Position');
set(figure8,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure8,'O2_distr_allreact','-dpdf','-r0');
   
hold off

figure9 = figure(9);

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2_ex(i,g) = interp1q(X_ex,n_i_N2_ex(:,g),xr(i));
    end
    semilogy(i_N2, u_N2_ex(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{N_2i}/n');
xlim([0,I(sw_o,1)]);
%ylim([1e-15,1])
xlabel('i');
title('N_2, без обмена');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',10);
set(figure9,'Units','Inches');
pos = get(figure9,'Position');
set(figure9,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure9,'N2_distr_ex','-dpdf','-r0');
hold off

figure10 = figure(10);

for i = 1 : length(xr)
    for g = 1 : length(i_O2)
        u_O2_ex(i,g) = interp1q(X_ex,n_i_O2_ex(:,g),xr(i));
    end
    semilogy(i_O2, u_O2_ex(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{O_2i}/n');
xlim([0,I(sw_o,2)]);
title('O_2, без обмена')
%ylim([1e-15,1])
xlabel('i');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',10);
set(figure10,'Units','Inches');
pos = get(figure10,'Position');
set(figure10,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure10,'O2_distr_ex','-dpdf','-r0');
   
hold off


