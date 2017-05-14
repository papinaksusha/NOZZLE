sp = 1;

p_cr = 100*101325;
T_cr = 8000;
x_N = 0.01;
x = 0 : x_N : 50;
x = [0,50];

e = e_i_c;
e_i = cell2mat(e(sp));
Z_vibr = sum(exp(-e_i/k/T_cr));

%% Бинарная 

init = zeros(I(sw_o,sp) + 1 + 3 , 1);

% Александрова

%n_at = 0.234; % 1 7000 N2
%n_at = 0.0264; % 100 7000 N2
%n_at = 0.29; % 1 4000 O2
%n_at = 0.0339; % 100 4000 O2
%n_at = 0.36; % 100 6000 O2

% Cиницын

% n_at = 0.12659;% 5 4000 O2
% n_at = 0.04192;% 50 4000 O2
% n_at = 0.80778;% 5 6000 O2
% n_at = 0.43712;% 50 6000 O2

% n_at = 0.00044;% 1 4000 N2
% n_at = 0.00004;% 100 4000 N2
% n_at = 0.48739;% 1 8000 N2
 n_at = 0.06579;% 100 8000 N2

n_mol = 1 - n_at;
n = zeros(1 , 5);
n(sp) = n_mol;
n(sp + 3) = n_at;

init(1 : I(sw_o , sp) + 1) = n_mol/Z_vibr.*exp(-e_i./k./T_cr);
init(I(sw_o , sp) + 2) = n_at;
init(I(sw_o , sp) + 3) = 1;
init(I(sw_o , sp) + 4) = 1;
v_cr = v_critical(n , T_cr);
v_cr = v_cr + v_cr*0.05;

%v_cr = 2000;

options = odeset('AbsTol', eps, 'RelTol', 2.22045e-14, 'OutputFcn', @odeplot, 'OutputSel', I(sw_o,sp) + 4);
%options = odeset('OutputFcn', @odeplot, 'OutputSel', I(sw_o,sp) + 4);

[X,Y] = Nozzle_2_full(x,init,options,T_cr,p_cr,v_cr,sp);

toc
%%
T = Y(: , I(sw_o , sp) + 4).*T_cr;
n_m = Y(: , 1 : I(sw_o , sp) + 1)./(sum(Y(:, 1 : I(sw_o , sp) + 2) , 2)*ones(1 , I(sw_o,sp) + 1));
n_a = Y(: , I(sw_o , sp) + 2)/sum(Y(:, 1 : I(sw_o , sp) + 2) , 2);

figure(1)
plot(X, T)

tit = ['N_2', 'O_2'];
xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

figure(2)

i_sp = 0 : I(sw_o,sp);
u = zeros(length(xr),length(i_sp));

for i = 1 : length(xr)
    for g = 1 : length(i_sp)
        u(i,g) = interp1q(X,n_m(:,g),xr(i));
    end
    semilogy(0 : I(sw_o,sp), u(i,:), 'color', colours(i,:)), hold on
end
%colormap jet
ylabel('n_i/n');
xlim([0,I(sw_o,sp)]);
%ylim([1e-15,1])
xlabel('i');
switch sp
    case 1
        title('N_2');
    case 2
        title('O_2');
end
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r =50');
hold off

switch sp
    case 1
        figure(3)
        semilogy(X, [n_m(:,1) n_m(:,2) n_m(:,4) n_m(:,6) n_m(:,11) n_m(:,21) n_m(:,31) n_m(:,41) n_m(:,46)])
        legend('0','1','3','5','10','20','30','40','45')
        title('N_2')
    case 2
        figure(3)
        semilogy(X, [n_m(:,1) n_m(:,2) n_m(:,4) n_m(:,6) n_m(:,11) n_m(:,21) n_m(:,31) n_m(:,33)])
        legend('0','1','3','5','10','20','30','32')
        title('O_2')
end
