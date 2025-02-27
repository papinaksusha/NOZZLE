function [f1,f2] = Nozzle_5_full_STELLAR_withoutNO(x,init,options,T_cr,p_cr,v_cr)

global k h c m theta_r D sw_n L sw_o k_ex_N2_STELLAR k_ex_O2_STELLAR Na

n_cr = p_cr/k/T_cr;

[f1,f2] = ode15s(@nozzle_5_STELLAR_NO_0,x,init,options);

    function dy = nozzle_5_STELLAR_NO_0(x,y)
    % STELLAR without NO excitation 
    
    l_N2 = L(1);
    l_O2 = L(2);
    l_mol = L(3);
    l_c = L(4);
    l_v = L(5);
    l_T = L(6);
    
    A = zeros(l_T , l_T); 
    b = zeros(l_T , 1);
    S = zeros(2);                                                          % S(1) - nozzle cross-section, S(2) - its derivative 
    
    % nozzle geometry
    
    switch sw_n
        case 1
             r_cr = 1e-3;
             alpha = 0.117*pi;
             S(1) = (1 + x*tan(alpha))^2;
             S(2) = 2*tan(alpha)*(1 + x*tan(alpha));
    end
    
    % dimensionless variables
    
    n_N2 = y(1 : l_N2);
    n_O2 = y(l_N2 + 1 : l_N2 + l_O2);
    n_NO = y(l_mol);
    n_N = y(l_mol + 1);
    n_O = y(l_c);
    v = y(l_v);
    T = y(l_T);
    M_cr = sum(m); 
    M = m./M_cr;
    
    n = [n_N2' n_O2' n_NO n_N n_O];
    rho = sum([M(1).*n_N2' M(2).*n_O2' M(3)*n_NO M(4)*n_N M(5)*n_O]);
    
    % vibrational energy of molecules
    
     e = e_i_c;
     e_i_N2 = cell2mat(e(1))./k./T_cr; %1*61
     e_i_O2 = cell2mat(e(2))./k./T_cr;
     %e_i_NO = cell2mat(e(3))./k./T_cr;
    
     % vibrational energy of 0th levels
     
     e_0_N2 = h*c/k/T_cr*1e-2*1175.78;
     e_0_O2 = h*c/k/T_cr*1e-2*787.38;
     e_0_NO = h*c/k/T_cr*1e-2*948.646642922263;
     e_0 = [e_0_N2.*ones(1 , l_N2) e_0_O2.*ones(1 , l_O2) e_0_NO 0 0];
     
     e_i = [e_i_N2 e_i_O2 0 0 0];
     
     % formation energy
     
     e_NO = (D(1)/2 + D(2)/2 - D(3))/k/T_cr;
     e_N = D(1)/2/k/T_cr;
     e_O = D(2)/2/k/T_cr;
     
     e_f = [zeros(1 , l_N2 + l_O2) e_NO e_N e_O];
     
    %% Left part, matrix A (A(y)*y = b)
    
    % kinetic equations
    
    A(1 : l_c , 1 : l_c) = diag(v.*ones(1 , l_c)); 
    A(1 : l_c , l_v) = n; 
    
    % momentum equation
    
    A(l_v, 1 : l_c) = T; 
    A(l_v, l_v) = M_cr*v_cr^2/k/T_cr * v *rho;
    A(l_v, l_T) = sum(n);
    
    % energy equation
    
    T_energy1 = [2.5*T.*ones(1 , l_mol) 1.5*T 1.5*T];
    T_energy2 = [3.5*T.*ones(1 , l_mol) 2.5*T 2.5*T];
    
    A(l_T, 1 : l_c) = T_energy1 + e_i + e_0 + e_f; 
    A(l_T, l_v) = 1/v*sum(n.*(T_energy2  + e_i + e_0 + e_f));
    A(l_T, l_T) = 2.5*(sum(n_N2) + sum(n_O2) + n_NO) + 1.5*(n_N + n_O);
   
    AA = sparse(A);
    
    %% Right part, vector b
    
    %%% reaction rate coefficients and source terms
        
    % dimensional variables
    
    n_N2_d = [0 n_N2'.*n_cr 0]'*ones(1 , 5); %50*5
    n_O2_d = [0 n_O2'.*n_cr 0]'*ones(1 , 5);
    n_NO_d = n_NO*n_cr;
    n_N_d = n_N*n_cr;
    n_O_d = n_O*n_cr;
    T_d = T*T_cr;
    
    %%% vibrational energy transitions (SSH)
    
    sw_o = 2;
    
    [k_N2_VT, k_N2_N2_VV, k_N2_O2_VV, ~] = k_ssh(1 , T_d);
    [k_O2_VT, k_O2_N2_VV, k_O2_O2_VV, ~] = k_ssh(2 , T_d);
    
    % interpoltion anharmonical on STELLAR dictribution
    
    e = e_i_c;
    
    e_i_N2_an = cell2mat(e(1))./k./T_cr;
    e_i_O2_an = cell2mat(e(2))./k./T_cr;
    
    k_N2_VT = interp1(e_i_N2_an(1 : end - 1), k_N2_VT, e_i_N2(1 : end - 1), 'spline');
    k_O2_VT = interp1(e_i_O2_an(1 : end - 1), k_O2_VT, e_i_O2(1 : end - 1), 'spline'); 

    k_N2_N2_VV = interp2(e_i_N2_an(1 : end - 1), e_i_N2_an(1 : end - 1)', ...
                 k_N2_N2_VV, e_i_N2(1 : end - 1), e_i_N2(1 : end - 1)', 'spline');
    k_N2_O2_VV = interp2(e_i_O2_an(1 : end - 1), e_i_N2_an(1 : end - 1)', ...
                 k_N2_O2_VV, e_i_O2(1 : end - 1), e_i_N2(1 : end - 1)', 'spline');
    
    k_O2_N2_VV = interp2(e_i_N2_an(1 : end - 1), e_i_O2_an(1 : end - 1)', ...
                 k_O2_N2_VV, e_i_N2(1 : end - 1), e_i_O2(1 : end - 1)', 'spline');
    k_O2_O2_VV = interp2(e_i_O2_an(1 : end - 1), e_i_O2_an(1 : end - 1)', ...
                 k_O2_O2_VV, e_i_O2(1 : end - 1), e_i_O2(1 : end - 1)', 'spline');
   
    sw_o = 3;
    
    % VT

    k_N2_VT_r = k_N2_VT.*exp(-(e_i_N2(2 : end) - e_i_N2(1 : end - 1))'*ones(1,5)/T); %��������� ������
    k_O2_VT_r = k_O2_VT.*exp(-(e_i_O2(2 : end) - e_i_O2(1 : end - 1))'*ones(1,5)/T);

    k_N2_VT = [zeros(1,5); k_N2_VT; zeros(1,5)]; %49*5
    k_O2_VT = [zeros(1,5); k_O2_VT; zeros(1,5)];
    k_N2_VT_r = [zeros(1,5); k_N2_VT_r; zeros(1,5)]; %49*5
    k_O2_VT_r = [zeros(1,5); k_O2_VT_r; zeros(1,5)];
    
    n_c_N2 = ones(l_N2 , 1)*[sum(n_N2) sum(n_O2) n_NO n_N n_O].*n_cr; %  %*  48*5
    
    R_N2_VT = sum(n_c_N2.*(n_N2_d(1 : end - 2,:).*k_N2_VT_r(1 : end - 1,:) - ...
                           n_N2_d(2 : end - 1,:).*k_N2_VT(1 : end - 1,:) + ...
                           n_N2_d(3 : end,:).*k_N2_VT(2 : end,:) - ...
                           n_N2_d(2 : end - 1,:).*k_N2_VT_r(2 : end,:)) , 2); %48*1
                     
    n_c_O2 = ones(l_O2 , 1)*[sum(n_N2) sum(n_O2) n_NO n_N n_O].*n_cr; 
    
    R_O2_VT = sum(n_c_O2.*(n_O2_d(1 : end - 2,:).*k_O2_VT_r(1 : end - 1,:) - ...
                           n_O2_d(2 : end - 1,:).*k_O2_VT(1 : end - 1,:) + ...
                           n_O2_d(3 : end,:).*k_O2_VT(2 : end,:) - ...
                           n_O2_d(2 : end - 1,:).*k_O2_VT_r(2 : end,:)) , 2);
               
    % VV 
    
    n_N2_d = [0 n_N2'.*n_cr 0]; %1*50
    n_O2_d = [0 n_O2'.*n_cr 0];
    
    k_N2_N2_VV_r = k_N2_N2_VV.*exp(-((e_i_N2(2 : end) - e_i_N2(1 : end - 1))'*ones(1 , l_N2 - 1) + ... 
                                   ones(l_N2 - 1 , 1)*(e_i_N2(1 : end - 1) - e_i_N2(2 : end)))./T);
    k_O2_O2_VV_r = k_O2_O2_VV.*exp(-((e_i_O2(2 : end) - e_i_O2(1 : end - 1))'*ones(1 , l_O2 - 1) + ... 
                                   ones(l_O2 - 1 , 1)*(e_i_O2(1 : end - 1) - e_i_O2(2 : end)))./T);

    k_N2_N2_VV = [zeros(1 , l_N2 - 1) ; k_N2_N2_VV; zeros(1 , l_N2 - 1)]; %49*47
    k_O2_O2_VV = [zeros(1 , l_O2 - 1) ; k_O2_O2_VV; zeros(1 , l_O2 - 1)];
    k_N2_N2_VV_r = [zeros(1 , l_N2 - 1) ; k_N2_N2_VV_r; zeros(1 , l_N2 - 1)]; %49*47
    k_O2_O2_VV_r = [zeros(1 , l_O2 - 1) ; k_O2_O2_VV_r; zeros(1 , l_O2 - 1)];
    
    R_N2_VV = sum(n_N2_d(1 : end - 2)'*n_N2_d(3 : end - 1).*k_N2_N2_VV_r(1 : end - 1, :) - ...
                  n_N2_d(2 : end - 1)'*n_N2_d(2 : end - 2).*k_N2_N2_VV(1 : end - 1, :) + ...
                  n_N2_d(3 : end)'*n_N2_d(2 : end - 2).*k_N2_N2_VV(2 : end, :) - ...
                  n_N2_d(2 : end - 1)'*n_N2_d(3 : end - 1).*k_N2_N2_VV_r(2 : end, :) , 2); % 48*1
    
    R_O2_VV = sum(n_O2_d(1 : end - 2)'*n_O2_d(3 : end - 1).*k_O2_O2_VV_r(1 : end - 1, :) - ...
                  n_O2_d(2 : end - 1)'*n_O2_d(2 : end - 2).*k_O2_O2_VV(1 : end - 1, :) + ...
                  n_O2_d(3 : end)'*n_O2_d(2 : end - 2).*k_O2_O2_VV(2 : end, :) - ...
                  n_O2_d(2 : end - 1)'*n_O2_d(3 : end - 1).*k_O2_O2_VV_r(2 : end, :) , 2);
       
    % VV'
              
    k_N2_O2_VV_r =  k_N2_O2_VV.*exp(-((e_i_N2(2 : end) - e_i_N2(1 : end - 1))'*ones(1 , l_O2 - 1) + ...
                                    ones(l_N2 - 1 , 1)*(e_i_O2(1 : end - 1) - e_i_O2(2 : end)))./T);
    k_O2_N2_VV_r =  k_O2_N2_VV.*exp(-((e_i_O2(2 : end) - e_i_O2(1 : end - 1))'*ones(1 , l_N2 - 1) + ...
                                    ones(l_O2 - 1, 1)*(e_i_N2(1 : end - 1) - e_i_N2(2 : end)))./T);                             
    
    k_N2_O2_VV = [zeros(1 , l_O2 - 1);  k_N2_O2_VV; zeros(1 , l_O2 - 1)];
    k_N2_O2_VV_r = [zeros(1 , l_O2 - 1);  k_N2_O2_VV_r; zeros(1 , l_O2 - 1)];
    k_O2_N2_VV = [zeros(1 , l_N2 - 1);  k_O2_N2_VV ; zeros(1 , l_N2 - 1)];
    k_O2_N2_VV_r = [zeros(1 , l_N2 - 1);  k_O2_N2_VV_r; zeros(1 , l_N2 - 1)];
    
    R_N2_VV_s = sum(n_N2_d(1 : end - 2)'*n_O2_d(3 : end - 1).*k_N2_O2_VV_r(1 : end - 1, :) - ...
                  n_N2_d(2 : end - 1)'*n_O2_d(2 : end - 2).*k_N2_O2_VV(1 : end - 1, :) + ...
                  n_N2_d(3 : end)'*n_O2_d(2 : end-2).*k_N2_O2_VV(2 : end, :) - ...
                  n_N2_d(2 : end - 1)'*n_O2_d(3 : end - 1).*k_N2_O2_VV_r(2 : end, :) , 2);
    
    R_O2_VV_s = sum(n_O2_d(1 : end - 2)'*n_N2_d(3 : end - 1).*k_O2_N2_VV_r(1 : end - 1, :) - ...
                  n_O2_d(2 : end - 1)'*n_N2_d(2 : end - 2).*k_O2_N2_VV(1 : end - 1, :) + ...
                  n_O2_d(3 : end)'*n_N2_d(2 : end-2).*k_O2_N2_VV(2 : end, :) - ...
                  n_O2_d(2 : end - 1)'*n_N2_d(3 : end - 1).*k_O2_N2_VV_r(2 : end, :) , 2);
       
    R_N2_vibr = R_N2_VT + R_N2_VV + R_N2_VV_s; 
    R_O2_vibr = R_O2_VT + R_O2_VV + R_O2_VV_s;
    
    % Chemical exchange
    
    n_N2_d = n_N2.*n_cr;
    n_O2_d = n_O2.*n_cr;
  
    TT = 100 : 100 : 100000;
    
    k_ex_N2 = interp2(TT , [0 : l_N2 - 1]',  squeeze(k_ex_N2_STELLAR(: , 1, :)) , T_d, [0 : l_N2 - 1]', 'spline');
    k_ex_O2 = interp2(TT , [0 : l_O2 - 1]' , squeeze(k_ex_O2_STELLAR(: , 1, :)) , T_d, [0 : l_O2 - 1]', 'spline');
%     
%     k_ex_N2loop = zeros(l_N2 , l_NO);
%     k_ex_O2loop = zeros(l_O2 , l_NO);
%     for q = 1 : l_NO
%         for v1 = 1 : l_N2
%             k_vector1 = reshape(k_ex_N2_STELLAR(v1,q,:),[1 1e3]);
%             k_ex_N2loop(v1,q) = interp1(TT,k_vector1,T_d,'spline');
%         end
%         for v2 = 1 : l_O2
%             k_vector2 = reshape(k_ex_O2_STELLAR(v2,q,:),[1 1e3]);
%             k_ex_O2loop(v2,q) = interp1(TT,k_vector2,T_d,'spline');
%         end
%     end

    k_ex_N2_r = k_ex_N2.*(m(1)*m(5)/m(3)/m(4))^(1.5).*theta_r(3)/theta_r(1)* ...
                0.5.*exp(-e_i_N2'./T + D(1)/k/T_d - D(3)/k/T_d);
    
    k_ex_O2_r = k_ex_O2.*(m(2)*m(4)/m(3)/m(5))^(1.5).*theta_r(3)/theta_r(2)* ...
                0.5.*exp(-e_i_O2'./T + D(2)/k/T_d - D(3)/k/T_d);
    
    R_N2_ex = n_NO_d*n_N_d.*k_ex_N2_r - n_N2_d.*n_O_d.*k_ex_N2; %48*1
    R_O2_ex = n_NO_d*n_O_d.*k_ex_O2_r - n_O2_d.*n_N_d.*k_ex_O2;
    R_NO_ex = - sum(R_N2_ex) - sum(R_O2_ex);
    R_N_ex = - sum(R_N2_ex) + sum(R_O2_ex);
    R_O_ex = sum(R_N2_ex) - sum(R_O2_ex);
    
    % Dissociation
    
    n_N2_d = n_N2.*n_cr*ones(1 , 5); 
    n_O2_d = n_O2.*n_cr*ones(1 , 5);
    
    k_diss_N2 = k_diss(1 , T_d);
    k_diss_O2 = k_diss(2 , T_d);
    k_diss_NO(1:3) = 0.41e12/Na*T_d^(-1)*exp(-D(3)/k/T_d);
    k_diss_NO(4:5) = 0.3e12/Na*T_d^(0.5)*exp(-D(3)/k/T_d);
    
    k_rec_N2 = k_diss_N2.*(m(1)/m(4)^2)^(1.5).*h^3.*(2*pi*k*T_d)^(-1.5).*...
                T_d./theta_r(1).*0.5.*exp(-e_i_N2'*ones(1,5)./T + D(1)/k/T_d); 
    k_rec_O2 = k_diss_O2.*(m(2)/m(5)^2)^(1.5).*h^3.*(2*pi*k*T_d)^(-1.5).*...
               T_d./theta_r(2).*0.5.*exp(-e_i_O2'*ones(1,5)./T + D(2)/k/T_d);
  
    %Z_NO_int = sum(exp(-e_i_NO./T))*T_d/theta_r(3);
    k_rec_NO = k_diss_NO*(m(3)/m(4)/m(5))^(1.5)*h^3*(2*pi*k*T_d)^(-1.5)* T_d./theta_r(3)*exp(D(3)/k/T_d);
    
    R_N2_diss = sum(n_c_N2.*(n_N_d^2.*k_rec_N2 - n_N2_d.*k_diss_N2) , 2);
    R_O2_diss = sum(n_c_O2.*(n_O_d^2.*k_rec_O2 - n_O2_d.*k_diss_O2) , 2);
    R_NO_diss = sum([sum(n_N2) sum(n_O2) n_NO n_N n_O].*n_cr.*(n_N_d*n_O_d*k_rec_NO - ...
                n_NO_d*k_diss_NO));
    R_N_diss = -2*sum(R_N2_diss) - R_NO_diss;
    R_O_diss = -2*sum(R_O2_diss) - R_NO_diss;
            
    
    R_N2 = r_cr/n_cr/v_cr.*(R_N2_vibr + R_N2_ex + R_N2_diss);
    R_O2 = r_cr/n_cr/v_cr.*(R_O2_vibr + R_O2_ex + R_O2_diss);
    R_NO = r_cr/n_cr/v_cr.*(R_NO_ex + R_NO_diss);
    R_N = r_cr/n_cr/v_cr.*(R_N_ex + R_N_diss);
    R_O = r_cr/n_cr/v_cr.*(R_O_ex + R_O_diss);
   
    %disp(sum(R_N2) + sum(R_O2) + sum(R_NO) + R_N + R_O);
    
    b(1 : l_N2) = R_N2 - v.*n_N2./S(1).*S(2); 
    b(l_N2 + 1 : l_N2 + l_O2) = R_O2 - v.*n_O2./S(1).*S(2);
    b(l_mol) = R_NO - v*n_NO/S(1)*S(2);
    b(l_mol + 1) = R_N - v*n_N/S(1)*S(2);
    b(l_c) = R_O - v*n_O/S(1)*S(2);
    b(l_v) = 0;
    b(l_T) = -S(2)/S(1)*sum(n.*(T_energy2  + e_i + e_0 + e_f));

    dy = AA\b;
    
    end

end