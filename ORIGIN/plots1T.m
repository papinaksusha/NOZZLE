clc
clear 
close all
format long e
set(0,'DefaultAxesFontSize',18);
set(0,'DefaultTextInterpreter', 'tex'); 
set(0,'DefaultTextFontSize',18);

I = [33 26 27; 47 36 39; 60, 45, 47]; 
sw_o = 2;
k = 1.380648e-23; 

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

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\1T_100_7000.mat')

n_N2_1t = nn2(6,:);
n_O2_1t = no2(6,:);
T_1t = T(6,:);
X_1t = X6;


e = e_i_c;
e_i_N2 = cell2mat(e(1));
e_i_O2 = cell2mat(e(2));

n_i_N2_1t = zeros(length(T_1t),I(sw_o,1)+1);
n_i_O2_1t = zeros(length(T_1t),I(sw_o,2)+1);
for i = 1 : length(T_1t)
n_i_N2_1t(i,:) = n_N2_1t(:,i)/sum(exp(-e_i_N2/k/T_1t(i))).*exp(-e_i_N2/k/T_1t(i)); 
n_i_O2_1t(i,:) = n_O2_1t(:,i)/sum(exp(-e_i_O2/k/T_1t(i))).*exp(-e_i_O2/k/T_1t(i)); 
end
xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

i_N2 = 0 : I(sw_o,1);
i_O2 = 0 : I(sw_o,2);

u_N2_3 = zeros(length(xr),length(i_N2));
u_O2_3 = zeros(length(xr),length(i_O2));

u_N2_1t = zeros(length(xr),length(i_N2));
u_O2_1t = zeros(length(xr),length(i_O2));

figure1 = figure(1);

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2_3(i,g) = interp1q(X_3,n_i_N2_3(:,g),xr(i));
    end
    semilogy(i_N2, u_N2_3(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{N_2i}/n');
xlim([0,I(sw_o,1)]);
ylim([1e-30,1])
xlabel('i');
title('N_2');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',10);
for i = 1 : length(xr)
    TT = interp1q(X_1t',T_1t',xr(i));
    u_N2_1t(i,:) = n_N2_1t(:,i)/sum(exp(-e_i_N2/k/TT)).*exp(-e_i_N2/k/TT);    
    semilogy(i_N2,  u_N2_1t(i,:), '--', 'color', colours(i,:)), hold on
end
set(figure1,'Units','Inches');
pos = get(figure1,'Position');
set(figure1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure1,'./PLOTS/N2_distr_1T','-dpdf','-r0');
hold off

figure2 = figure(2);

for i = 1 : length(xr)
    for g = 1 : length(i_O2)
        u_O2_3(i,g) = interp1q(X_3,n_i_O2_3(:,g),xr(i));
    end
    semilogy(i_O2, u_O2_3(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_{O_2i}/n');
xlim([0,I(sw_o,2)]);
title('O_2')
ylim([1e-15,1])
xlabel('i');
legend({'x/{r^*} = 0', 'x/{r^*} = 1', 'x/{r^*} = 2', 'x/{r^*} = 3', 'x/{r^*} = 5', 'x/{r^*} = 10', ...
       'x/{r^*} = 15', 'x/{r^*} = 20', 'x/{r^*} = 25', 'x/{r^*} = 30', 'x/{r^*} = 40', 'x/{r^*} = 50'},'location','southwest','Fontsize',10);
for i = 1 : length(xr)
    TT = interp1q(X_1t',T_1t',xr(i));
    u_O2_1t(i,:) = n_O2_1t(:,i)/sum(exp(-e_i_O2/k/TT)).*exp(-e_i_O2/k/TT);
   semilogy(i_O2,  u_O2_1t(i,:), '--', 'color', colours(i,:)),hold on
end
   
set(figure2,'Units','Inches');
pos = get(figure2,'Position');
set(figure2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure2,'./PLOTS/O2_distr_1T','-dpdf','-r0');

figure3 = figure(3);
hold on
plot(X_1t, T_1t, X_3, T_3,'linewidth',2)
legend('однотемпературное','поуровневое')
xlabel('x/r^*');
ylabel('T, K')
hold off   


i_N2 = [1,3,6,11,21,31,41,46];
colours = colormap(lines(length(i_N2)));

figure4 = figure(4);
for i = 1 : length(i_N2)
    semilogy(X_3, n_i_N2_3(:,i_N2(i)), 'color', colours(i,:),'linewidth',2), hold on
    semilogy(X_1t, n_i_N2_1t(:,i_N2(i)),'--', 'color', colours(i,:),'linewidth',2)
end
title('N_2')
xlabel('x/r^*');
ylabel('n_{N2i}/n');
ylim([1e-30,1])
set(figure4,'Units','Inches');
pos = get(figure4,'Position');
set(figure4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure4,'./PLOTS/N2_distr_1T_alongx','-dpdf','-r0');

hold off


figure5 = figure(5);
semilogy(X_3, n_i_O2_3(:,1), X_3, n_i_O2_3(:,2), X_3, n_i_O2_3(:,6), X_3, n_i_O2_3(:,11), X_3, n_i_O2_3(:,21),X_3, n_i_O2_3(:,31),X_3, n_i_O2_3(:,33),'linewidth',2), hold on
semilogy(X_1t, n_i_O2_1t(:,1),'--', X_1t, n_i_O2_1t(:,2),'--', X_1t, n_i_O2_1t(:,6),'--', X_1t, n_i_O2_1t(:,11),'--', X_1t, n_i_O2_1t(:,21),'--',X_1t, n_i_O2_1t(:,31),'--',X_1t, n_i_O2_1t(:,33),'--','linewidth',2)
title('O_2')
xlabel('x/r^*');
ylabel('n_{O2i}/n');
ylim([1e-10,1])
set(figure5,'Units','Inches');
pos = get(figure5,'Position');
set(figure5,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figure5,'./PLOTS/O2_distr_1T_alongx','-dpdf','-r0');
hold off

save('1T_100_7000_distr.mat', 'u_N2_1t','u_O2_1t','X_1t')