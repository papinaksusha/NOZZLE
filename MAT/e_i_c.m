function e = e_i_c()

global w wx h c sw_o I

i = {0 : I(sw_o,1), 0 : I(sw_o,2), 0 : I(sw_o,3)};
e = cell(3,1);

for sp = 1 : 3
    e(sp) = {h*c.*(w(sp).*cell2mat(i(sp)) - wx(sp).*cell2mat(i(sp)) - ...
            wx(sp).*cell2mat(i(sp)).^2)};   % Vibrational energy of molecules
end

end