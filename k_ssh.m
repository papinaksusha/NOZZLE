function [VT,VV_N2, VV_O2, VV_NO] = k_ssh(sp,T)
% посмотреть для гармонического более подробно в лосеве
% Детальный баланс писать в функции правых частей
% ВСЕ ОПИСАНИЕ ПЕРЕПИСАТЬ
% VT exchange rate coefficients k_{c,i+1->i}^M, SSH theory
% INPUT:
% sp - chemical specie sp =1,2,3
% T - temperature
% OUTPUT:
% f(1) =
% f(2) = 
% порядки с Олей не совсем сходятся
global m c k h w wx sw_o I

d = [3.621 3.458 3.47 3.298 2.75].*1e-10; % диаметры молекул m
% wave_vibr.pdf
r0 = [3.621 3.458 3.47 3.298 2.75].*1e-10; % m
r0 = diag(r0);
eps = [97.5 107.4 119 71.4 80]; % eps/k, K
eps = diag(eps);

for y = 1 : 5
    for j = 1 : 5
        eps(y,j) = sqrt(eps(y,y)*eps(j,j));
        r0(y,j) = 0.5*(r0(y,y) + r0(j,j));
    end
end
% VT

Z0 = 3; % В статье есть формула. она есть и в ступоченко, но написано, что для столкновения одинаковых частиц

mu = m(sp).*m./(m(sp) + m); % kg
omega = 2*pi*c*w(sp); % c-1 %в справочнике omega = 2 pi c (w - 2wx)
alpha = 17.5./r0(sp,:); % m-1
hi = ((0.5*pi^2*omega^2/k/T).*mu./(alpha.^2)).^(1/3);
R = (0.5.*sqrt(1 + hi.*T./eps(sp,:)) + 0.5).^(-1/3);                       % (r/r0)^2
P10 = 1.294.*R./Z0./(1 + 1.1.*eps(sp,:)./T).*8.*pi^3.*mu.*omega./alpha.^2./h.* ... 
      sqrt(4*pi/3.*hi).*exp(-3.*hi + h*omega/4/pi/k/T + eps(sp,:)./T);

  % НЕ ВЕРНО!  но у Оли так
%Zn = 2.*(d(sp) + d).^2.*sqrt(2*pi*k*T./mu);
%Zn = (d(sp) + d).^2.*sqrt(2*pi*k*T./mu)./8;
Zn = (d(sp) + d).^2.*sqrt(pi*k*T/2./mu);
k10 = P10.*Zn;

i = {0 : I(sw_o,1) - 1, 0 : I(sw_o,2) - 1, 0 : I(sw_o,3) - 1};
    
switch sw_o
        case 1
            VT = (cell2mat(i(sp)) + 1)'*k10;
        case 2
            E1 = (w(sp) - 2*wx(sp))*0.014;
            gamma0 = 0.32./alpha.*1e10.*sqrt(mu.*6.022e26/T).*E1;
            g_n = gamma_n(sp,T);
            deltaVT = zeros(I(sw_o,sp),5);
            l = find(g_n >= 20);
            ll = find(g_n < 20);
            deltaVT(l) = 4.*gamma0(fix((l - 1)./I(sw_o,sp)) + 1).^(2/3).*wx(sp).*0.014./E1;
            deltaVT(ll) = 4/3.*gamma0(fix((ll - 1)./I(sw_o,sp)) + 1).*wx(sp).*0.014./E1;

            VT = (cell2mat(i(sp)) + 1)'*k10.*exp((ones(1,5)'*cell2mat(i(sp)))'.*deltaVT);
end

% VV & VV'

m_r(1) = m(4)*m(4)/(m(4)+m(4));
m_r(2) = m(5)*m(5)/(m(5)+m(5));
m_r(3) = m(4)*m(5)/(m(4)+m(5));
lambda1 = 0.5;
lambda2 = 0.5;

Q10 = lambda1^2*lambda2^2*4.*alpha(1:3).^2.*k.*T./omega^2./m_r(sp);
k1001 = Q10.*Zn(1:3);

VV = cell(1,3);

switch sw_o
        case 1
            for l = 1:3
                VV(l) = {(cell2mat(i(sp)) + 1)'*(cell2mat(i(l)) + 1).*k1001(l)};
            end
        case 2
            for l = 1:3
                if (l == sp)
                    deltaVV = 0.427/alpha(sp).*1e10*sqrt(mu(sp).*6.0221366516752e26/T)*wx(sp)*0.01438776978978;
                    Abs = zeros(I(sw_o,sp),I(sw_o,sp));
                    Abs(1,:) = 0 : I(sw_o,sp) - 1;
                    Abs(:,1) = 0 : I(sw_o,sp) - 1;
                    for q = 2 : I(sw_o,sp)
                        for g = 2 : I(sw_o,sp)
                            Abs(q,g) = abs(Abs(1,g) - Abs(q,1));
                        end
                    end
                    VV(l) = {(cell2mat(i(sp)) + 1)'*(cell2mat(i(l)) + 1).*k1001(l).* ...
                            exp(-deltaVV.*Abs).*(1.5 - 0.5.*exp(-deltaVV.*Abs))};
                else
                    deltaC = 0.427/alpha(sp)*1e10*sqrt(mu(sp)*6.022e26/T)*wx(sp)*0.014;
                    deltaD = 0.427/alpha(l)*1e10*sqrt(mu(l)*6.022e26/T)*wx(l)*0.014;
                    p = (w(sp) - 2*wx(sp) - w(l) + 2*wx(l))/2/wx(sp);
                    Abs = zeros(I(sw_o,l) + 1, I(sw_o,sp) + 1);
                    Abs(1,2:end) = 0 : I(sw_o,sp) - 1;
                    Abs(2:end,1) = 0 : I(sw_o,l) - 1;
                    for q = 2 : I(sw_o,l) + 1
                        for g = 2 : I(sw_o,sp) + 1
                            Abs(q,g) = abs(deltaD*Abs(q,1) - deltaC*Abs(1,g) + deltaC*p);
                        end
                    end
                    Abs = Abs(2:end,2:end);
                    VV(l) = {(cell2mat(i(sp)) + 1)'*(cell2mat(i(l)) + 1).*k1001(l).* ...
                            exp(-Abs').*exp(deltaC.*p)};
                end
            end
end

VV_N2 = cell2mat(VV(1));
VV_O2 = cell2mat(VV(2));
VV_NO = cell2mat(VV(3));

end