clc
clear 
close all
format long e
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextInterpreter', 'tex'); 
set(0,'DefaultTextFontSize',18);

I = [33 26 27; 47 36 39; 60, 45, 47]; 
sw_o = 2;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\1_7000_5comp_withoutNO_ex3.mat')


n_i_N2_c = n_i_N2;
n_i_O2_c = n_i_O2;
n_N2_c = n_N2;
n_O2_c = n_O2;
n_NO_c = n_NO;
n_N_c = n_N;
n_O_c = n_O;
T_c = T;
v_c = v;
X_c = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\1_7000_5comp_withoutNO_ex3_hyper.mat')


n_i_N2_h = n_i_N2;
n_i_O2_h = n_i_O2;
n_N2_h = n_N2;
n_O2_h = n_O2;
n_NO_h = n_NO;
n_N_h = n_N;
n_O_h = n_O;
T_h = T;
v_h = v;
X_h = X;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\1_7000_5comp_withoutNO_ex3_f.mat')


n_i_N2_f = n_i_N2;
n_i_O2_f = n_i_O2;
n_N2_f = n_N2;
n_O2_f = n_O2;
n_NO_f = n_NO;
n_N_f = n_N;
n_O_f = n_O;
T_f = T;
v_f = v;
X_f = X;


figure1 = figure(1);
plot(X_c, T_c, 'k', X_h, T_h, 'b', X_f, T_f, 'r', 'linewidth', 2)
xlabel('x/r^*')
ylabel('T, K')
legend({'коническое','гиперболическое','F4'})
set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure1,'T_profiles','-dpdf','-r0');
hold off

figure2 = figure(2);
plot(X_c, n_N2_c, 'k', X_h, n_N2_h, 'b', X_f, n_N2_f, 'r', 'linewidth', 2)
xlabel('x/r^*')
ylabel('n_{N_2}/n')
xlim([0,5]);
legend({'коническое','гиперболическое','F4'})
set(figure2,'Units','Inches');
pos = get(figure2,'Position');
set(figure2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure2,'N2_profiles','-dpdf','-r0');
hold off

figure3 = figure(3);
plot(X_c, n_O2_c, 'k', X_h, n_O2_h, 'b', X_f, n_O2_f, 'r', 'linewidth', 2)
xlabel('x/r^*')
xlim([0,5]);
ylabel('n_{O_2}/n')
legend({'коническое','гиперболическое','F4'})
set(figure3,'Units','Inches');
pos = get(figure3,'Position');
set(figure3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure3,'O2_profiles','-dpdf','-r0');
hold off

figure4 = figure(4);
plot(X_c, n_N_c, 'k', X_h, n_N_h, 'b', X_f, n_N_f, 'r', 'linewidth', 2)
xlabel('x/r^*')
xlim([0,5]);
ylabel('n_{N}/n')
legend({'коническое','гиперболическое','F4'})
set(figure4,'Units','Inches');
pos = get(figure4,'Position');
set(figure4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure4,'N_profiles','-dpdf','-r0');
hold off


%%

xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

i_N2 = 0 : I(sw_o,1);
i_O2 = 0 : I(sw_o,2);
u_N2_c = zeros(length(xr),length(i_N2));
u_O2_c = zeros(length(xr),length(i_O2));
u_N2_h = zeros(length(xr),length(i_N2));
u_O2_h = zeros(length(xr),length(i_O2));
u_N2_f = zeros(length(xr),length(i_N2));
u_O2_f = zeros(length(xr),length(i_O2));

figure5 = figure(5);

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2_c(i,g) = interp1q(X_c,n_i_N2_c(:,g),xr(i));
    end
    semilogy(i_N2, u_N2_c(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{N_2i}/n');
xlim([0,I(sw_o,1)]);
%ylim([1e-15,1])
xlabel('i');
title('N_2, коническое');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',12);
set(figure5,'Units','Inches');
pos = get(figure5,'Position');
set(figure5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure5,'N2_distr_с','-dpdf','-r0');
hold off

figure6 = figure(6);

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2_h(i,g) = interp1q(X_h,n_i_N2_h(:,g),xr(i));
    end
    semilogy(i_N2, u_N2_h(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{N_2i}/n');
xlim([0,I(sw_o,1)]);
%ylim([1e-15,1])
xlabel('i');
title('N_2, гиперболическое');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',12);
set(figure6,'Units','Inches');
pos = get(figure6,'Position');
set(figure6,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure6,'N2_distr_h','-dpdf','-r0');
hold off


figure7 = figure(7);

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2_f(i,g) = interp1q(X_f,n_i_N2_f(:,g),xr(i));
    end
    semilogy(i_N2, u_N2_f(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{N_2i}/n');
xlim([0,I(sw_o,1)]);
%ylim([1e-15,1])
xlabel('i');
title('N_2, F4');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',12);
set(figure7,'Units','Inches');
pos = get(figure7,'Position');
set(figure7,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure7,'N2_distr_f','-dpdf','-r0');
hold off


figure8 = figure(8);

for i = 1 : length(xr)
    for g = 1 : length(i_O2)
        u_O2_c(i,g) = interp1q(X_c,n_i_O2_c(:,g),xr(i));
    end
    semilogy(i_O2, u_O2_c(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{O_2i}/n');
xlim([0,I(sw_o,2)]);
title('O_2, коническое')
%ylim([1e-15,1])
xlabel('i');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',12);
set(figure8,'Units','Inches');
pos = get(figure8,'Position');
set(figure8,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure8,'O2_distr_с','-dpdf','-r0');

figure9 = figure(9);

for i = 1 : length(xr)
    for g = 1 : length(i_O2)
        u_O2_h(i,g) = interp1q(X_h,n_i_O2_h(:,g),xr(i));
    end
    semilogy(i_O2, u_O2_h(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{O_2i}/n');
xlim([0,I(sw_o,2)]);
title('O_2, гиперболическое')
%ylim([1e-15,1])
xlabel('i');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',12);
set(figure9,'Units','Inches');
pos = get(figure9,'Position');
set(figure9,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure9,'O2_distr_h','-dpdf','-r0');
hold off

figure10 = figure(10);

for i = 1 : length(xr)
    for g = 1 : length(i_O2)
        u_O2_f(i,g) = interp1q(X_f,n_i_O2_f(:,g),xr(i));
    end
    semilogy(i_O2, u_O2_f(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{O_2i}/n');
xlim([0,I(sw_o,2)]);
title('O_2, F4')
%ylim([1e-15,1])
xlabel('i');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',12);
set(figure10,'Units','Inches');
pos = get(figure10,'Position');
set(figure10,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure10,'O2_distr_f','-dpdf','-r0');