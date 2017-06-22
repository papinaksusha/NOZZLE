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


f1 = figure(1);
plot(X_3, n_N2_3,'k', X_ex, n_N2_ex,'b', X_rec, n_N2_rec,'r','linewidth',2), hold on
legend({'все реакции','без обмена','без рекомбинации'}, 'Position',[0.438690486248759 0.297619134708057 0.437499989941716 0.215873009912552])
%semilogy(X_3, n_O2_3,'r-', X_ex, n_O2_ex, 'r--', X_rec, n_O2_rec, 'r:','linewidth',2)
title('N_2')
xlim([0,5]);
xlabel('x/r^*')
ylabel('n_{N_2}/n')
set(f1,'Units','Inches');
pos = get(f1,'Position');
set(f1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f1,'N2_react','-dpdf','-r0');
hold off

f2 = figure(2);
semilogy(X_3, n_N_3,'k', X_ex, n_N_ex,'b', X_rec, n_N_rec,'r','linewidth',2), hold on
legend({'все реакции','без обмена','без рекомбинации'})%, 'Position',[0.438690486248759 0.297619134708057 0.437499989941716 0.215873009912552])
title('N')
xlabel('x/r^*')
xlim([0,5]);
ylabel('n_N/n')
set(f2,'Units','Inches');
pos = get(f2,'Position');
set(f2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f2,'N_react','-dpdf','-r0');
hold off

f3 = figure(3);
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
set(f3,'Units','Inches');
pos = get(f3,'Position');
set(f3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f3,'all_all','-dpdf','-r0');
hold off

f4 = figure(4);
semilogy(X_ex, n_N2_ex,'linewidth',2), hold on
semilogy(X_ex, n_O2_ex,'linewidth',2)
semilogy(X_ex, n_NO_ex,'linewidth',2)
semilogy(X_ex, n_N_ex,'linewidth',2)
semilogy(X_ex, n_O_ex,'linewidth',2)
title('без обмена')
xlim([0,5]);
ylim([1e-8,1]);
xlabel('x/r^*');
ylabel('n_c/n');
legend('N_2','O_2','NO','N','O');
set(f4,'Units','Inches');
pos = get(f4,'Position');
set(f4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f4,'all_ex','-dpdf','-r0');
hold off


f4 = figure(4);
plot(X_3, n_O2_3,'k', X_ex, n_O2_ex,'b', X_rec, n_O2_rec,'r','linewidth',2), hold on
legend({'все реакции','без обмена','без рекомбинации'}, 'Position',[0.438690486248759 0.297619134708057 0.437499989941716 0.215873009912552])
%semilogy(X_3, n_O2_3,'r-', X_ex, n_O2_ex, 'r--', X_rec, n_O2_rec, 'r:','linewidth',2)
title('O_2')
xlim([0,5]);
xlabel('x/r^*')
ylabel('n_{O_2}/n')
set(f4,'Units','Inches');
pos = get(f4,'Position');
set(f4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f4,'O2_react','-dpdf','-r0');
hold off