function f = k_diss(sp,T)

 %%%!!!!! Не различается для столкновения с атомом или молекулой !!!!!%%%%
 
% Dissociation rate coefficients, Treanor-Marrone model
% INPUT:
% sp - chemical specie
% T - temperature
% OUTPUT:
% f(1) = k_N2i_diss
% f(2) = k_O2i_diss
% f(3) = k_rec_N2i
% f(4) = k_rec_O2i

global h c w wx D k I sw_o Na
% Constants for Arrenius law from
% International Workshop on Radiation of High Temperature Gases in Atmospheric
% Entry. Part II // 30 Sep.- 10 Oct.2004. Porquerolles, France.

%N = [-1 -1 -1];
%A = [4.15e-11 1.51e-11 6.81e-12];
N = [-1.6 -1.5]; % Olya disser
A = [3.7e15 1.99e15]./Na;

i = 0 : I(sw_o,sp);
%f = zeros(I(sw_o,sp));

e = h*c.*(w(sp).*i - wx(sp).*i - wx(sp).*i.^2); 

                       % Vibrational energy of molecules
Z_vibr = sum(exp(-e./k./T));                                           % Equilibrium vibrational partition function                                               
k_c_diss_eq = A(sp)*T^N(sp)*exp(-D(sp)/k/T);                              % Thermal equilibrium dissociation rate coefficient, m^3/s
U = D(sp)/6/k;                                                          % Parameter (U = inf, U = D/6k, U = 3T)
Z_ci = Z_vibr/sum(exp(e./k./U))*exp(e(i+1)/k*(1/T + 1/U));             % Nonequilibrium factor Z_ci
f = Z_ci*k_c_diss_eq;

end

