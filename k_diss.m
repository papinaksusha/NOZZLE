function f = k_diss(sp,T)
% Dissociation rate coefficients, Treanor-Marrone model
% INPUT:
% sp - chemical specie
% T - temperature

global D k I sw_o Na
% Constants for Arrenius law from
% International Workshop on Radiation of High Temperature Gases in Atmospheric
% Entry. Part II // 30 Sep.- 10 Oct.2004. Porquerolles, France.


% как брал Синицын
A = [1.6e16 3.7e15; 
     9.99e15 1.99e15;
     0.3e12 0.41e12]./Na; % A{(sp,:) = [Aat Amol]
N = [-1.6 -1.6;
     -1.5 -1.5;
     0.5 -1];
 
e = e_i_c;
e = cell2mat(e(sp));

Z_vibr = sum(exp(-e./k./T));                                           % Equilibrium vibrational partition function                                               
k_c_diss_eq = A(sp,:).*T.^N(sp,:)*exp(-D(sp)/k/T);                              % Thermal equilibrium dissociation rate coefficient, m^3/s
U = D(sp)/6/k;                                                          % Parameter (U = inf, U = D/6k, U = 3T)
%U = Inf;
Z_ci = Z_vibr/sum(exp(e./k./U)).*exp(e./k*(1/T + 1/U));             % Nonequilibrium factor Z_ci
f = zeros(I(sw_o,sp) + 1, 5);
f(:,1:3) = Z_ci'*k_c_diss_eq(2)*ones(1,3);
f(:,4:5) = Z_ci'*k_c_diss_eq(1)*ones(1,2);

end

