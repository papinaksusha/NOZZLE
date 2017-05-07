clc
clear 
close all
format long e

set(gca,'FontSize',20);
set(0, 'DefaultAxesFontSize',16)

global I sw_o w wx h c k

Na = 6.022141e23;                                                          % mol^-1, Avogadro constant 
k = 1.380648e-23;                                                          % J/K, Boltzmann constant
h = 6.626070e-34;                                                          % J*s, Planck constant
c = 299792458;                                                             % m/s, speed of light 
w = [235857 158019 190420];                                                % m^-1, spectroscopic constant
%wx = [1432 1198 1407.5];                                                   % m^-1, spectroscopic constant
m = 1.6605402e-27*[28.0134 31.99880 30.00610 14.0067 15.9994];             % kg, molecular mass
theta_r = [2.86 2.07 2.42];                                                % K, characteristic rotational emperature
D = h*c.*[7.871e6 4.126e6 5.24e6];                                         % J, dissociation energy

sw_o = 2;

switch sw_o 
    case 1
        wx = [0 0 0];   
    case 2
        wx = [1432 1198 1407.5];   
end

x_N = 0.01;
x = 0 : x_N : 50;

q15 = load('1_5000_war_anhar_full');
q1005 = load('100_5000_war_anhar_full');
q17 = load('1_7000_war_anhar_full');
q1007 = load('100_7000_war_anhar_full');

T15 = q15.Y(:,end).*5000;
T1005 = q1005.Y(:,end).*5000;
T17 = q17.Y(:,end).*7000;
T1007 = q1007.Y(:,end).*7000;

v_cr1 = 1.496343752444559e+03;

v_cr2 = 1.768253770314827e+03;



v15 = q15.Y(:,end-1).*v_cr1;
v1005 = q1005.Y(:,end-1).*v_cr1;
v17 = q17.Y(:,end-1).*v_cr2;
v1007 = q1007.Y(:,end-1).*v_cr2;

figure(1)
plot(q15.X, [T15 T1005 T17 T1007],'Linewidth',2)
xlabel('x/r^*');
ylabel('T');
legend('T^* = 5000 K, p^* = 1 atm','T^* = 5000 K, p^* = 100 atm',...
       'T^* = 7000 K, p^* = 1 atm','T^* = 7000 K, p^* = 100 atm');
   
figure(2)
plot(q15.X, [v15 v1005 v17 v1007],'Linewidth',2)
xlabel('x/r^*');
ylabel('v');
legend('T^* = 5000 K, p^* = 1 atm','T^* = 5000 K, p^* = 100 atm',...
       'T^* = 7000 K, p^* = 1 atm','T^* = 7000 K, p^* = 100 atm','Location','SouthEast');

   
n_i_N2_17 = q17.Y(: , 1 : I(sw_o,1) + 1)./(sum(q17.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,1) + 1));
n_i_O2_17 = q17.Y(: , I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)./(sum(q17.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,2) + 1));   
n_N2_17 = sum(n_i_N2_17,2);
n_O2_17 = sum(n_i_O2_17,2);
n_NO_17 = q17.Y(: , I(sw_o,1) + I(sw_o,2) + 3)./sum(q17.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_N_17 = q17.Y(: , I(sw_o,1) + I(sw_o,2) + 4)./sum(q17.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_O_17 = q17.Y(: , I(sw_o,1) + I(sw_o,2) + 5)./sum(q17.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);


figure(3)
semilogy(q17.X, [n_N2_17 n_O2_17 n_NO_17 n_N_17 n_O_17],'LineWidth',1)
legend('N_2','O_2','NO','N','O');
title('T^* = 7000 K, p^* = 1 atm');
ylabel('n_c/n');
xlabel('x/_r*');
xlim([0,5]);

n_i_N2_1007 = q1007.Y(: , 1 : I(sw_o,1) + 1)./(sum(q1007.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,1) + 1));
n_i_O2_1007 = q1007.Y(: , I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)./(sum(q1007.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,2) + 1));   
n_N2_1007 = sum(n_i_N2_1007,2);
n_O2_1007 = sum(n_i_O2_1007,2);
n_NO_1007 = q1007.Y(: , I(sw_o,1) + I(sw_o,2) + 3)./sum(q1007.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_N_1007 = q1007.Y(: , I(sw_o,1) + I(sw_o,2) + 4)./sum(q1007.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_O_1007 = q1007.Y(: , I(sw_o,1) + I(sw_o,2) + 5)./sum(q1007.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);

xr = [0 1 2 3 5 10 15 20 25 30 40 50];


colours = colormap(jet(length(xr)));

figure(4)

for i = 1 : length(xr)
    semilogy(0 : I(sw_o,1), n_i_N2_1007((xr(i) - x(1))./x_N + 1,:),'color', colours(i,:)), hold on
end
ylabel('n_i/n');
xlim([0,I(sw_o,1)]);
%ylim([1e-15,1])
xlabel('i');
title('N_2');
legend({'x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50'},'Fontsize',10);
hold off

figure(5)

for i = 1 : length(xr)
    semilogy(0 : I(sw_o,2), n_i_O2_1007((xr(i) - x(1))/x_N + 1,:),'color', colours(i,:)), hold on
end
ylabel('n_i/n');
xlim([0,I(sw_o,2)]);
%ylim([1e-15,1])
xlabel('i');
title('O_2');
legend({'x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50'},'Location','SouthWest','Fontsize',10);
%hold off


n_i_N2_15 = q15.Y(: , 1 : I(sw_o,1) + 1)./(sum(q15.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,1) + 1));
n_i_O2_15 = q15.Y(: , I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)./(sum(q15.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,2) + 1));   
n_N2_15 = sum(n_i_N2_15,2);
n_O2_15 = sum(n_i_O2_15,2);
n_NO_15 = q15.Y(: , I(sw_o,1) + I(sw_o,2) + 3)./sum(q15.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_N_15 = q15.Y(: , I(sw_o,1) + I(sw_o,2) + 4)./sum(q15.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_O_15 = q15.Y(: , I(sw_o,1) + I(sw_o,2) + 5)./sum(q15.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);

n_i_N2_1005 = q1005.Y(: , 1 : I(sw_o,1) + 1)./(sum(q1005.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,1) + 1));
n_i_O2_1005 = q1005.Y(: , I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)./(sum(q1005.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,2) + 1));   
n_N2_1005 = sum(n_i_N2_1005,2);
n_O2_1005 = sum(n_i_O2_1005,2);
n_NO_1005 = q1005.Y(: , I(sw_o,1) + I(sw_o,2) + 3)./sum(q1005.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_N_1005 = q1005.Y(: , I(sw_o,1) + I(sw_o,2) + 4)./sum(q1005.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_O_1005 = q1005.Y(: , I(sw_o,1) + I(sw_o,2) + 5)./sum(q1005.Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);

q1t = load('100_7000_u_1');
n_N21t = q1t.Y5(1,:)./sum(q1tY(1:5,:));
n_O21t = q1t.Y5(2,:)./sum(q1tY(1:5,:));
n_N21ti = ;
n_O21ti = ;