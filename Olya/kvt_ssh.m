function k_down = kvt_ssh(t,sp)
% function for vib. exch. rates, SHH model, Stupochenko
% t - temperature, K
% sw_o - oscillator switcher
% k_down - array of the rate coef-s of VT exchanges A2(i+1)+M = A2(i)+M
% i+1 - vib. level before collision
% i - vib. level after collision

global k c h h_bar om_e om_x_e l sw_o m r0 em re
% constantes
if sw_o == 1 % anh.os. 
    om10 = om_e(sp)-2*om_x_e(sp); % lenear oscillation frequency
elseif sw_o == 2 % har.os.
    om10 = om_e(sp);
end
Lmax = l(sp)-1; % number of maximum level
i = (0:Lmax-1);
om0 = 2*pi*c*om10;% circular oscillation frequency, sec^-1

% reduced masses, kg
mu = m(sp)*m./(m(sp)+m);

% inverse radius in 1st approximation, m^-1
a = 17.5./r0(sp,:);

% collision frequency, m^3/sec
Z = r0(sp,:).^2./sqrt(mu)*sqrt(8*pi*k*t);

% firstly, find rates of transition from 1st lev. to 2nd lev.
chi = (pi^2*om0^2/2/k*(mu./a.^2)*t^(-1)).^(1/3);% dim-less
%
r = r0(sp,:) .* (0.5 * sqrt(1+(chi*t)./em(sp,:)) + 0.5).^(-1/6);
% steric factor
Z_0 = (a .* re(sp)).^2 .* exp(-((3*a .* re(sp).^2)./(8*r)));

Const = 1.294./Z_0*4*pi^2*om0/h_bar*sqrt(4*pi/3);
%
r_c_r_02 = (0.5*sqrt(1+(chi*t)./em(sp,:))+0.5).^(-1/3);
%
P10 = Const.*r_c_r_02.*(1+1.1*em(sp,:)/t).^(-1).*...
    mu./a.^2.*chi.^0.5.*exp(-3*chi+h_bar*om0/(2*k*t)+em(sp,:)/t);
k10 = Z.*P10; % m^3/sec

% secondly, find rates coef-s for transition i+1 -> i
if sw_o == 2
    k_down = k10'*(i+1);
elseif sw_o == 1
    % anharmonicity factor for VT transitions (Gordietz)
    aA = a * 1e-10; % A^-1
    % reduced masses
    mu_amu = mu/1.6605e-27; % amu
    % adiabatic factor for transition i+1 -> i
    % E(i+1) - E(i)
    diffE = (om_e(sp)-2*(i+1)*om_x_e(sp))*h*c; % J
    %
    E1 = (om_e(sp)-2*om_x_e(sp))*1.4388e-2; % K
    dE = om_x_e(sp)*1.4388e-2; % K
    %
    gamma_n = pi ./ a / h_bar .* sqrt(mu/(2*k*t));
    gamma_n = gamma_n' * diffE;
    %
    gamma_0 = 0.32 ./ aA .* sqrt(mu_amu/t)*E1;
    %
    delta = zeros(5,length(i));
    k_down = zeros(5,length(i));
    for iM = 1:5
    delta(iM,:) = (4/3*gamma_0(iM) * dE/E1).^(gamma_n(iM,:) < 20).*...
        (4*(gamma_0(iM))^(2/3) * dE/E1).^(gamma_n(iM,:) >= 20);
    %
    k_down(iM,:) = (i+1) * k10(iM) .* exp(i .* delta(iM,:)) .* exp(-i*h*c*om_x_e(sp) / (k*t));
    end
end