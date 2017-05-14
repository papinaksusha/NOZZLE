function v_cr = v_critical(n,T)
% n - вектор пятикомпонентный
global m k Na w wx c h sw_o

Molar = [28.0134 31.99880 30.00610 14.0067 15.9994].*1e-3;
M_inverse = sum(m.*n./Molar);
R_bar = k*Na*M_inverse;

Ctr = 1.5*R_bar;
Crot = k*sum(n(1:3));

e = e_i_c;

e_i_N2 = cell2mat(e(1));
e_i_O2 = cell2mat(e(2));
e_i_NO = cell2mat(e(3));

if sw_o == 3
    e_0_N2 = e_i_N2(1);
    e_0_O2 = e_i_O2(1);
    e_0_NO = e_i_NO(1);
else
    e_0_N2 = h*c*(w(1)/2 - wx(1)/4);
    e_0_O2 = h*c*(w(2)/2 - wx(2)/4);
    e_0_NO = h*c*(w(3)/2 - wx(3)/4);
end

Z_vibr_N2 = sum(exp(-e_i_N2./k./T));
Z_vibr_O2 = sum(exp(-e_i_O2./k./T));
Z_vibr_NO = sum(exp(-e_i_NO./k./T));
E_vibr_N2 = sum((e_i_N2 + e_0_N2).*exp(-e_i_N2./k./T))/Z_vibr_N2;
E_vibr_O2 = sum((e_i_O2 + e_0_O2).*exp(-e_i_O2./k./T))/Z_vibr_O2;
E_vibr_NO = sum((e_i_NO + e_0_NO).*exp(-e_i_NO./k./T))/Z_vibr_NO;
C_vibr_N2 = T^(-2)/k*(sum((e_i_N2 + e_0_N2).*e_i_N2.*exp(-e_i_N2./k./T))/Z_vibr_N2 - ...
                 E_vibr_N2*sum(e_i_N2.*exp(-e_i_N2./k./T))/Z_vibr_N2);
C_vibr_O2 = T^(-2)/k*(sum((e_i_O2 + e_0_O2).*e_i_O2.*exp(-e_i_O2./k./T))/Z_vibr_O2 - ...
                 E_vibr_O2*sum(e_i_O2.*exp(-e_i_O2./k./T))/Z_vibr_O2);
C_vibr_NO = T^(-2)/k*(sum((e_i_NO + e_0_NO).*e_i_NO.*exp(-e_i_NO./k./T))/Z_vibr_NO - ...
                 E_vibr_NO*sum(e_i_NO.*exp(-e_i_NO./k./T))/Z_vibr_NO);
     
Cvibr = sum([C_vibr_N2 C_vibr_O2 C_vibr_NO].*n(1:3));
Cv = Ctr + Crot + Cvibr;

Cp = R_bar + Cv;

kappa = Cp/Cv;

v_cr = sqrt(kappa*R_bar/sum(m.*n)*T);


end

