e = e_i_c;
e_i_N2 = cell2mat(e(1));
e_i_O2 = cell2mat(e(2));
Z_vibr_N2 = sum(exp(-e_i_N2/k/T_cr));
Z_vibr_O2 = sum(exp(-e_i_O2/k/T_cr));

init = zeros(I(sw_o,1) + 1 + I(sw_o,2) + 1 + 5,1);

n_N2_cr = 0.79;
n_O2_cr = 0.21;

%init(1) = 0.79;
%init(I(sw_o,1) + 2) = 0.21;

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

figure(1)
plot(X,Y(:,end).*T_cr);
%save('MAT/1_5000_war_anhar_full.mat','X','Y');


%%
n_i_N2 = Y(: , 1 : I(sw_o,1) + 1)./(sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,1) + 1));
n_i_O2 = Y(: , I(sw_o,1) + 2 : I(sw_o,1) + I(sw_o,2) + 2)./(sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2)*ones(1,I(sw_o,2) + 1));
xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

figure(2)

for i = 1 : length(xr)
    semilogy(0 : I(sw_o,1), n_i_N2((xr(i) - X(1))/x_N + 1,:), 'color', colours(i,:)), hold on
end
%colormap jet
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

n_N2 = sum(n_i_N2,2);
n_O2 = sum(n_i_O2,2);
n_NO = Y(: , I(sw_o,1) + I(sw_o,2) + 3)./sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_N = Y(: , I(sw_o,1) + I(sw_o,2) + 4)./sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);
n_O = Y(: , I(sw_o,1) + I(sw_o,2) + 5)./sum(Y(: , 1 : I(sw_o,1) + I(sw_o,2) + 5) , 2);

figure(4)
semilogy(X, [n_N2 n_O2 n_NO n_N n_O])
legend('N_2','O_2','NO','N','O');
ylabel('n_c/n');
xlabel('x/_r*');
xlim([0,5]);