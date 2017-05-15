global k_ex_N2_STELLAR k_ex_O2_STELLAR L

l_N2 = I(sw_o , 1) + 1;
l_O2 = I(sw_o , 2) + 1;
l_NO = I(sw_o , 3) + 1;
l_mol =  l_N2 + l_O2 + l_NO;
l_c = l_mol + 2;
l_v = l_c + 1;
l_T = l_v + 1;

L = [l_N2 l_O2 l_NO l_mol l_c l_v l_T];

x_N = 1e-2;
x = 0 : x_N : 50;
T_cr = 7000;
p_cr = 100*101325;
n_cr = p_cr/k/T_cr; 
n_N2_cr = 0.79;
n_O2_cr = 0.21;
n_NO_cr = 0;

load('./MAT/N2-O-NO-N-QCT-3D-Rates.mat');
load('./MAT/O2-N-NO-O-QCT-3D-Rates.mat');

k_ex_N2_STELLAR = 1.e-6.*N2_O_NO_N_Rates;
k_ex_O2_STELLAR = 1.e-6.*O2_N_NO_O_Rates;

e = e_i_c;
e_i_N2 = cell2mat(e(1))';
e_i_O2 = cell2mat(e(2))';
e_i_NO = cell2mat(e(3))';
Z_vibr_N2 = sum(exp(-e_i_N2/k/T_cr));
Z_vibr_O2 = sum(exp(-e_i_O2/k/T_cr));
Z_vibr_NO = sum(exp(-e_i_NO/k/T_cr));

init = zeros(l_T , 1);

init(1 : l_N2) = n_N2_cr/Z_vibr_N2.*exp(-e_i_N2./k./T_cr);
init(l_N2 + 1 : l_N2 + l_O2) = n_O2_cr/Z_vibr_O2.*exp(-e_i_O2./k./T_cr);
init(l_N2 + l_O2 + 1 : l_mol) = n_NO_cr/Z_vibr_NO.*exp(-e_i_NO./k./T_cr);
init(l_v) = 1;
init(l_T) = 1;
v_cr = v_critical([sum(init(1 : l_N2)) sum(init(l_N2 + 1 : l_N2 + l_O2)) sum(l_N2 + l_O2 + 1 : l_mol)...
                   init(l_mol + 1) init(l_c)] , T_cr);
v_cr = v_cr + v_cr*0.2;

%v_cr = 2000;

%options = odeset('AbsTol', eps, 'RelTol', 2.3e-14, 'OutputFcn', @odeplot, 'OutputSel', l_T);
options = odeset('AbsTol', 1e-54, 'RelTol', 2.3e-14,'Stats', 'on', 'OutputFcn', @odeplot, 'BDF', 'off', 'OutputSel', l_T);

[X,Y] = Nozzle_5_full_STELLAR(x,init,options,T_cr,p_cr,v_cr);

toc
%%
n_i_N2 = Y(: , 1 : l_N2)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_N2));
n_i_O2 = Y(: , l_N2 + 1 : l_N2 + l_O2)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_O2));
n_i_NO = Y(: , l_N2 + l_O2 + 1 : l_mol)./(sum(Y(: , 1 : l_c) , 2)*ones(1 , l_NO));
n_N2 = sum(n_i_N2 , 2);
n_O2 = sum(n_i_O2 , 2);
n_NO = sum(n_i_NO , 2);
n_N = Y(: , l_mol + 1)./sum(Y(: , 1 : l_c) , 2);
n_O = Y(: , l_c)./sum(Y(: , 1 : l_c) , 2);
v = Y(: , l_v).*v_cr;
T = Y(: , l_T).*T_cr;

figure(1)
plot(X , T);
%save('MAT/1_5000_war_anhar_full.mat','X','Y');

xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

i_N2 = 0 : I(sw_o,1);
i_O2 = 0 : I(sw_o,2);
i_NO = 0 : I(sw_o,3);

u_N2 = zeros(length(xr) , length(i_N2));
u_O2 = zeros(length(xr) , length(i_O2));
u_NO = zeros(length(xr) , length(i_NO));

figure(2)

for i = 1 : length(xr)
    for g = 1 : length(i_N2)
        u_N2(i,g) = interp1q(X,n_i_N2(:,g),xr(i));
    end
    semilogy(i_N2, u_N2(i,:), 'color', colours(i,:)), hold on
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
    for g = 1 : length(i_O2)
        u_O2(i,g) = interp1q(X,n_i_O2(:,g),xr(i));
    end
    semilogy(i_O2, u_O2(i,:), 'color', colours(i,:)), hold on
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
ylim([1e-5,1]);


figure(5)

for i = 1 : length(xr)
    for g = 1 : length(i_NO)
        u_NO(i,g) = interp1q(X,n_i_NO(:,g),xr(i));
    end
    semilogy(i_NO, u_NO(i,:), 'color', colours(i,:)), hold on
end

ylabel('n_i/n');
xlim([0,I(sw_o,3)]);
%ylim([1e-15,1])
xlabel('i');
title('NO');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50');
hold off

%% conservation laws

N_init = [sum(init(1 : l_N2)) sum(init(l_N2 + 1 : l_N2 + l_O2)) ...
          sum(init(1 : l_mol)) init(l_mol + 1) init(l_c)].*n_cr; 
n_i_N2_d = Y(: , 1 : l_N2).*n_cr;
n_i_O2_d = Y(: , l_N2 + 1 : l_N2 + l_O2).*n_cr;
n_i_NO_d = Y(: , l_N2 + l_O2 + 1 : l_mol).*n_cr;
n_N2_d = sum(n_i_N2_d , 2);
n_O2_d = sum(n_i_O2_d , 2);
n_NO_d = sum(n_i_NO_d , 2);
n_N_d = Y(: , l_mol + 1).*n_cr;
n_O_d = Y(: , l_c).*n_cr;
N = [n_N2_d n_O2_d n_NO_d n_N_d n_O_d];
rho = sum(ones(length(X),1)*m.*N , 2);

%  continuity equation rho*v*S = const (S - conical)
 
r_cr = 1e-3;
r = r_cr.*(1 + X.*tan(0.117*pi)); 
 
Q_cr = sum(m.*N_init)*v_cr*r_cr^2;
Q = rho.*r.^2.*v;
eps1 = max(abs(Q - Q_cr)./Q_cr);
 
% % energy equation
 
% R_c = k*Na./molar;
% rho_c = (ones(length(X),1)*m.*N);
% Y_c =  rho_c./(rho*ones(1,5));
% 
% h_c(:,1) = 3.5*R_c(1).*T + 1./rho_c(:,1).*sum(ones(length(X),1)*(e_i_N2 + e_0_N2_STELLAR).*n_i_N2_d , 2);
% h_c(:,2) = 3.5*R_c(2).*T + 1./rho_c(:,2).*sum(ones(length(X),1)*(e_i_O2 + e_0_O2_STELLAR).*n_i_O2_d , 2);
% h_c(:,3) = 3.5*R_c(3).*T + 1./rho_c(:,3).*sum(ones(length(X),1)*(e_i_NO + e_0_NO_STELLAR).*n_i_NO_d , 2) +...
%            (D(1)/2 + D(2)/2 - D(3))/m(3);
% h_c(:,4) = 2.5*R_c(4).*T + D(1)/2;
% h_c(:,5) = 2.5*R_c(5).*T + D(2)/2;
%  
% H = v.^2./2 + sum(Y_c.*h_c , 2);
% eps2 = max(abs(H - H(1))./H(1));