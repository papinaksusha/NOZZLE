clc
clear 
close all
format long e

global sw_o h c I w wx

sw_o = 2;                                                                  % switch on oscillator, 1 -  harmonic oscillator; 2 -  anharmonic oscillator 3 - STELLAR
switch sw_o 
    case 1
        wx = [0 0 0];   
    case 2
        wx = [1432 1198 1407.5];                                           % m^-1, spectroscopic constant
    case 3
        wx = [1432 1198 1407.5];   
end
w = [235857 158019 190420]; 
molar = 1e-3.*[28.0134 31.99880 30.00610 14.0067 15.9994];                 % kg/mol, molar mass         
m = 1.6605402e-24.*molar; 
h = 6.626070e-34;                                                          % J*s, Planck constant
c = 299792458;  
I = [33 26 27; 47 36 39; 60 45 47]; 
l_N2 = I(sw_o , 1) + 1;
l_O2 = I(sw_o , 2) + 1;
l_mol =  l_N2 + l_O2 + 1;
l_c = l_mol + 2;
l_v = l_c + 1;
l_T = l_v + 1;

load('C:\Users\Ksusha\Desktop\S-t-s\Matlab\ORIGIN\100_7000_5comp_withoutNO_ex3.mat')
xr = [0 1 2 3 5 10 15 20 25 30 40 50];

colours = colormap(jet(length(xr)));

i_N2 = 0 : I(sw_o,1);
i_O2 = 0 : I(sw_o,2);

u_N2 = zeros(length(xr),length(i_N2));
u_O2 = zeros(length(xr),length(i_O2));

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
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r = 50');
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
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r = 50');
hold off

figure(4)
semilogy(X, [n_N2 n_O2 n_NO n_N n_O])
legend('N_2','O_2','NO','N','O');
ylabel('n_c/n');
xlabel('x/_r*');
xlim([0,5]);
ylim([1e-5,1]);

u_N2 = u_N2';
i_N2 = i_N2';
u_O2 = u_O2';
i_O2 = i_O2';

e = e_i_c();
e_i_N2 = cell2mat(e(1))'.*6.242e18; %1*48
e_i_O2 = cell2mat(e(2))'.*6.242e18;
    
figure(5)

for i = 1 : length(xr)
    semilogy(e_i_N2, u_N2(:,i), 'color', colours(i,:)), hold on
end

ylabel('n_i/n');
xlim([0,e_i_N2(end)]);
%ylim([1e-15,1])
xlabel('\epsilon_i^{N_2}, eV');
title('N_2');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r = 50');
hold off

figure(6)

for i = 1 : length(xr)
    semilogy(e_i_O2, u_O2(:,i), 'color', colours(i,:)), hold on
end

ylabel('n_i/n');
xlim([0,e_i_O2(end)]);
%ylim([1e-15,1])
xlabel({'\epsilon_i^{O_2}, eV'});
title('O_2');
legend('x/r = 0', 'x/r = 1', 'x/r = 2', 'x/r = 3', 'x/r = 5', 'x/r = 10', ...
       'x/r = 15', 'x/r = 20', 'x/r = 25', 'x/r = 30', 'x/r = 40', 'x/r = 50');
hold off

save('100_7000_distrib.mat','e_i_N2','e_i_O2','u_N2','u_O2','X')