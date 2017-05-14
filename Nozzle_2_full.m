function [f1,f2] = Nozzle_2_full(x,init,options,T_cr,p_cr,v_cr,sp)

global k h c w wx m theta_r D I sw_o sw_n

n_cr = p_cr/k/T_cr;

[f1,f2] = ode15s(@nozzle_2,x,init,options);

    function dy = nozzle_2(x,y)
    % Без возбуждения NO
    
    l = I(sw_o,sp) + 1;
    
   % dy = zeros(l_N2 + 1 + l_O2 + 1 + 5 , 1); %  90 *1
    A = zeros(l + 3 , l + 3); 
    b = zeros(l + 3 , 1);
    S = zeros(2); % S(1) сама площадь, S(2) производная
    
    % геом параметры сопла
    switch sw_n
        case 1
             r_cr = 1e-3;
    end
    
    % безразмерные переменные

    n_mol = y(1 : l);
    n_at = y(l + 1);
    v = y(l + 2);
    T = y(l + 3);
    M_cr = m(sp) + m(sp + 3); 
    M = [m(sp) m(sp + 3)]./M_cr;
    
    n = [n_mol' n_at];
    rho = sum([M(1).*n_mol' M(2).*n_at]);
    
    switch sw_n
        case 1
            alpha = 0.117*pi;
            S(1) = (1+x*tan(alpha))^2;
            S(2) = 2*tan(alpha)*(1+x*tan(alpha));
    end
    
    % вектора колебательных энергий молекул
     e = e_i_c;
     e_i = cell2mat(e(sp))./k./T_cr; %1*48
     e_0 = h*c/k/T_cr*(w(sp)/2 - wx(sp)/4);
    
     % энергии образования
     
     e_f = D(sp)/2/k/T_cr;
     
     %% ЛЧ, матрица А

    for j = 1 : l
        A(j , j) = v;     % уравнения кинетики
        A(j , l + 2) = n(j);     % уравнения кинетики
        A(l + 2 , j) = T;  % уравнение импульса
        A(l + 3 , j) = 2.5*T + e_i(j) + e_0;  % уравнение энергии
    end
    
    A(l + 1 , l + 1) = v; % уравнения кинетики
    A(l + 1 , l + 2) = n(l + 1);  % уравнения кинетики
    A(l + 2 , l + 1) = T; % уравнение импульса
    
    % уравнение импульса
    
    A(l + 2 , l + 2) = M_cr*v_cr^2/k/T_cr * v * rho;
    A(l + 2 , l + 3) = sum(n);
    
    % уравнение энергии
    
    A(l + 3 , l + 1) = 1.5*T + e_f;
    A(l + 3 , l + 2) = 1/v*(sum((3.5*T + e_i + e_0).*n_mol') + (2.5*T + e_f)*n_at);
    A(l + 3 , l + 3) = 2.5*sum(n_mol) + 1.5*n_at;

    AA = sparse(A);
    
    %% ПЧ вектор b
    % коэффициенты
        % размерные переменные
    
    n_mol_d = [0 n_mol'.*n_cr 0]'*ones(1,2); %50*5
    n_at_d = n_at*n_cr;
    T_d = T*T_cr;
    
    % колебательный обмен

    [k_VT, k_VV_N2, k_VV_O2, ~] = k_ssh(sp,T_d);
    
    switch sp
        case 1
            k_VV = k_VV_N2;
        case 2
            k_VV = k_VV_O2;
    end
    %%% VT
    
    k_VT = [k_VT(: , sp) k_VT(: , sp + 3)];
   
    i = 0 : l - 2;

    k_VT_r = k_VT.*exp(-h*c/k/T_cr.*(w(sp) - 2*wx(sp) - 2*wx(sp).*i)'*ones(1,2)./T); %детальный баланс
    
    k_VT = [zeros(1,2); k_VT; zeros(1,2)]; %49*5
    k_VT_r = [zeros(1,2); k_VT_r; zeros(1,2)];
    
    n_c = ones(l , 1)*[sum(n_mol) n_at].*n_cr; % мб убрать n-cr  а потом др коэфф обезраз? короче, подумать как лучше обезраз релакс члены
    %*  48*5
    R_VT = sum(n_c.*(n_mol_d(1 : end - 2,:).*k_VT_r(1 : end - 1,:) - ...
                        n_mol_d(2 : end - 1,:).*k_VT(1 : end - 1,:) + ...
                        n_mol_d(3 : end,:).*k_VT(2 : end,:) - ...
                        n_mol_d(2 : end - 1,:).*k_VT_r(2 : end,:)) , 2); %48*1
                    
    %%% VV
    
    n_mol_d = [0 n_mol'.*n_cr 0]; %1*50
    
    k_VV_r = k_VV.*exp(-2*h*c*wx(sp)/k/T_cr.*(i'*ones(1,l-1) - ... 
                                   ones(l-1,1)*i)./T);

    k_VV = [zeros(1 , l - 1) ; k_VV; zeros(1, l - 1)]; %49*47
    k_VV_r = [zeros(1 , l - 1) ; k_VV_r; zeros(1 , l - 1)]; %49*47
    
    R_VV = sum(n_mol_d(1 : end - 2)'*n_mol_d(3 : end - 1).*k_VV_r(1 : end - 1, :) - ...
                  n_mol_d(2 : end - 1)'*n_mol_d(2 : end - 2).*k_VV(1 : end - 1, :) + ...
                  n_mol_d(3 : end)'*n_mol_d(2 : end - 2).*k_VV(2 : end, :) - ...
                  n_mol_d(2 : end - 1)'*n_mol_d(3 : end - 1).*k_VV_r(2 : end, :) , 2); % 48*1
    
    R_vibr = R_VT + R_VV; %48*1
    
    % диссоциация
    
    n_mol_d = n_mol.*n_cr*ones(1,2); %48*2
    
    k_dis = k_diss(sp,T_d); % 48*5 
    k_dis = [k_dis(:,sp) k_dis(:,sp+3)];
%    k_dis = k_dis'*ones(1,2) ;% потому что не различается для партнера 48*5
   
    k_rec = k_dis.*(m(sp)/m(sp + 3)^2)^(1.5).*h^3.*(2*pi*k*T_d)^(-1.5).*...
                T_d./theta_r(sp).*0.5.*exp(-e_i'*ones(1,2)./T + D(sp)/k/T_d); 

    R_mol_dis = sum(n_c.*(n_at_d^2.*k_rec - n_mol_d.*k_dis) , 2);
    R_at_dis = -2*sum(R_mol_dis);

    R_mol = r_cr/n_cr/v_cr.*(R_vibr + R_mol_dis);
    R_at = r_cr/n_cr/v_cr.*R_at_dis;
   
    b(1 : l) = R_mol - v.*n_mol./S(1).*S(2); %%% далее индексы поправить
    b(l + 1) = R_at - v.*n_at/S(1)*S(2);
    b(l + 2) = 0;
    b(l + 3) = -S(2)/S(1)*(sum((3.5*T + e_i + e_0).*n_mol') + (2.5*T + e_f)*n_at);
    %%
    
    dy = AA\b;
    
    end

end