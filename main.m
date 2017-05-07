%% 5-components air mixture N2/O2/NO/N/O, state-to-state kinetics in nozzle  

clc
clear 
close all
format long e
global Na k h c w wx m theta_r D I sw_o mol sw_n ex_model

%% CONSTs
% J = kg*m^2/s^2, Pa = kg/m/s^2
% [N2 O2 NO N O]

tic

Na = 6.022141e23;                                                          % mol^-1, Avogadro constant 
k = 1.380648e-23;                                                          % J/K, Boltzmann constant
h = 6.626070e-34;                                                          % J*s, Planck constant
c = 299792458;                                                             % m/s, speed of light 
w = [235857 158019 190420];                                                % m^-1, spectroscopic constant
sw_o = 2;                                                                  % swith on oscillator, 1 -  harmonic oscillator; 2 -  anharmonic oscillator
switch sw_o 
    case 1
        wx = [0 0 0];   
    case 2
        wx = [1432 1198 1407.5];                                           % m^-1, spectroscopic constant
end
molar = 1e-3.*[28.0134 31.99880 30.00610 14.0067 15.9994];                 % kg/mol, molar mass         
m = 1.6605402e-24.*molar;                                                  % kg, molecular mass
theta_r = [2.86 2.07 2.42];                                                % K, characteristic rotational emperature
D = h*c.*[7.871e6 4.126e6 5.24e6];                                         % J, dissociation energy
sw_n = 1; % switcn on nozzle

mol = 3; %1 - бинарная, 3 - пятикомпонентная    надо ли???                                                                       % сделать диалоговое окно?
% Number of vibrational levels                                                 
% I(1,:) - harmonic oscillator, I(2,:) - anharmonic oscillator
ex_model = 3; %выбор модели реакции обмена
I = [33 26 27; 47 36 39]; 

%% Проблемы с диссоциацией О2, VT O2, переписать все через вызов e_i_c, один раз кэйс на свитч осциллятора с wx = [000]?
% переписать детальный балана через уровни
% 1 - 'N2' и тд
fig = 1;

% N_T = 100;
% TT = 2000 : N_T : 14000;
%
% TEST_SSH(fig,N_T,TT);
% TEST_DISS_TM(fig,N_T,TT);
% TEST_EXCHANGE(fig,N_T,TT)

%%
%initial values
x_N = 0.01;
x = 0 : x_N : 50;
T_cr = 7000;
p_cr = 100*101325;
n_cr = p_cr/k/T_cr; 
n_N2_cr = 0.79;
n_O2_cr = 0.21;

e = e_i_c;
e_i_N2 = cell2mat(e(1));
e_i_O2 = cell2mat(e(2));
e_i_NO = cell2mat(e(3));
Z_vibr_N2 = sum(exp(-e_i_N2/k/T_cr));
Z_vibr_O2 = sum(exp(-e_i_O2/k/T_cr));

init = zeros(I(sw_o,1) + 1 + I(sw_o,2) + 1 + 5,1);

init(1 : I(sw_o,1) + 1) = n_N2_cr/Z_vibr_N2.*exp(-e_i_N2./k./T_cr);
init(I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2) = n_O2_cr/Z_vibr_O2.*exp(-e_i_O2./k./T_cr);
init(I(sw_o,1) + I(sw_o,2) + 6) = 1;
init(I(sw_o,1) + I(sw_o,2) + 7) = 1;
v_cr = v_critical([sum(init(1 : I(sw_o,1) + 1)) sum(init(I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)) ...
                   init(I(sw_o,1) + I(sw_o,2) + 4) init(I(sw_o,1) + I(sw_o,2) + 4) init(I(sw_o,1) + I(sw_o,2) + 5)],T_cr);
v_cr = v_cr + v_cr*0.1;

%v_cr = 2000;

options = odeset('AbsTol', 1e-54, 'RelTol', 2.3e-14, 'OutputFcn', @odeplot, 'OutputSel', I(sw_o,1) + I(sw_o,2) + 7);

[X,Y] = Nozzle_5_full(x,init,options,T_cr,p_cr,v_cr);

toc

n_i_N2 = Y(: , 1 : I(sw_o,1) + 1)./(sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,1) + 1));
n_i_O2 = Y(: , I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)./(sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,2) + 1));
n_N2 = sum(n_i_N2 , 2);
n_O2 = sum(n_i_O2 , 2);
n_NO = Y(: , I(sw_o,1) + I(sw_o,2) + 3)./sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_N = Y(: , I(sw_o,1) + I(sw_o,2) + 4)./sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_O = Y(: , I(sw_o,1) + I(sw_o,2) + 5)./sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
v = Y(: , I(sw_o,1) + I(sw_o,2) + 6).*v_cr;
T = Y(: , end).*T_cr;

figure(1)
plot(X , T);
%save('MAT/1_5000_war_anhar_full.mat','X','Y');


%%

xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

figure(2)

for i = 1 : length(xr)
    semilogy(0 : I(sw_o,1), n_i_N2((xr(i) - X(1))/x_N + 1,:), 'color', colours(i,:)), hold on
end

ylabel('n_i/n');
xlim([0,I(sw_o,1)]);
%ylim([1e-15,1])
xlabel('i');
title('N_2');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50');
hold off

figure(3)

for i = 1 : length(xr)
    semilogy(0 : I(sw_o,2), n_i_O2((xr(i) - X(1))/x_N + 1,:), 'color', colours(i,:)), hold on
end
ylabel('n_i/n');
xlim([0,I(sw_o,2)]);
%ylim([1e-15,1])
xlabel('i');
title('O_2');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50');
hold off

figure(4)
semilogy(X, [n_N2 n_O2 n_NO n_N n_O])
legend('N_2','O_2','NO','N','O');
ylabel('n_c/n');
xlabel('x/_r*');
xlim([0,5]);

%% conservation laws

N_init = [sum(init(1 : I(sw_o,1) + 1)) sum(init(I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)) ...
          init(I(sw_o,1) + I(sw_o,2) + 3 : I(sw_o,1) + I(sw_o,2) + 5)'].*n_cr;
n_i_N2_d = Y(: , 1 : I(sw_o,1) + 1).*n_cr;
n_i_O2_d = Y(: , I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2).*n_cr;
n_N2_d = sum(n_i_N2_d , 2);
n_O2_d = sum(n_i_O2_d , 2);
n_NO_d = Y(: , I(sw_o,1) + I(sw_o,2) + 3).*n_cr;
n_N_d = Y(: , I(sw_o,1) + I(sw_o,2) + 4).*n_cr;
n_O_d = Y(: , I(sw_o,1) + I(sw_o,2) + 5).*n_cr;
N = [n_N2_d n_O2_d n_NO_d n_N_d n_O_d];
rho = sum(ones(length(x),1)*m.*N , 2);

%  continuity equation rho*v*S = const (S - conical)

r_cr = 1e-3;
r = r_cr.*(1 + x.*tan(0.117*pi)); 

Q_cr = sum(m.*N_init)*v_cr*r_cr^2;
Q = rho'.*r.^2.*v';
eps1 = max(abs(Q - Q_cr)./Q_cr);

% energy equation

R_c = k*Na./molar;
rho_c = (ones(length(x),1)*m.*N);
Y_c =  rho_c./(rho*ones(1,5));

h_c(:,1) = 3.5*R_c(1).*T + 1./rho_c(:,1).*sum(ones(length(x),1)*(e_i_N2 + h*c*(w(1)/2 - wx(1)/4)).*n_i_N2_d , 2);
h_c(:,2) = 3.5*R_c(2).*T + 1./rho_c(:,2).*sum(ones(length(x),1)*(e_i_O2 + h*c*(w(2)/2 - wx(2)/4)).*n_i_O2_d , 2);
h_c(:,3) = 3.5*R_c(3).*T + h*c*(w(3)/2 - wx(3)/4)/m(3) + (D(1)/2 + D(2)/2 - D(3))/m(3);
h_c(:,4) = 2.5*R_c(4).*T + D(1)/2;
h_c(:,5) = 2.5*R_c(5).*T + D(2)/2;

H = v.^2./2 + sum(Y_c.*h_c , 2);
%eps2 = max(abs(H - H(1))./H(1));