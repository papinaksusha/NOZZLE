function kdown = kvvs_ssh(t,sp1,sp2)
% rate coef-s of VV exchanges
% kdown - array of k_(i->i-1)^(j->j+1)
% t - temperature, K

global c h k m om_e om_x_e l sw_o r0

m1 = m(sp1); % mass of 1st molecule, kg
m2 = m(sp2); % mass of 2nd molecule, kg
mu = m1 * m2 / (m1+m2); % reduced mass
if (sp1 == 1)||(sp1 == 2)
    m_osc = 0.25 * m1; % reduced mass of osc-or
elseif sp1 == 3
    m_osc = m(4)*m(5)/(m(4)+m(5));
end
if sw_o == 1 % number of vibr. levels
    om10_1 = om_e(sp1)-2*om_x_e(sp1); % lenear oscillation frequency, anh.os.
    om10_2 = om_e(sp2)-2*om_x_e(sp2);
elseif sw_o == 2
    om10_1 = om_e(sp1);
    om10_2 = om_e(sp2);
end
Lmax1 = l(sp1)-1; Lmax2 = l(sp2)-1; % max level
om0 = 2*pi*c*[om10_1 om10_2];% circular oscillation frequency, sec^-1
%
a = 17.5 / r0(sp1,sp2); % inverse radius in 1st approximation m^-1
Z = r0(sp1,sp2)^2 / sqrt(mu) * sqrt(8*pi*k*t); % collision frequency, m^3/sec
%
Q10 = 0.5^4 / m_osc * a^2 / om0(1)^2 * 4 * k * t;
k10 = Q10 * Z;

kdown = zeros(Lmax1,Lmax2);
j_up = 0:Lmax2-1;
if sw_o == 2
    for i_down = 1:Lmax1
        kdown(i_down,:) = i_down * (j_up+1) * k10;
    end
elseif sw_o == 1
    % anharmonicity factor for VV transitions
    a_1 = 17.5 / r0(sp1,sp1) * 1e-10; % A^-1
    a_2 = 17.5 / r0(sp2,sp2) * 1e-10; % A^-1
    mu_amu = mu / 1.6605e-27; % amu
    dE_1 = om_x_e(sp1)*1.4388e-2; % K
    dE_2 = om_x_e(sp2)*1.4388e-2; % K
    delta_1 = 0.427 / a_1 * sqrt(mu_amu / t) * dE_1; % Gorgietz book
    delta_2 = 0.427 / a_2 * sqrt(mu_amu / t) * dE_2;
    p = (om_e(sp1)-om_e(sp2)-2*(om_x_e(sp2)-om_x_e(sp1)))/(2*om_x_e(sp1));
    %
    for i_down = 1:Lmax1
        kdown(i_down,:) = i_down * (j_up+1) * k10 .* ...
            exp(-abs(delta_2 * j_up - delta_1 * (i_down-1-p))) * exp(delta_1 * p) .* ...
            exp((j_up * om_x_e(sp1) - (i_down-1) * om_x_e(sp2)) * h * c / (k * t));
    end
end