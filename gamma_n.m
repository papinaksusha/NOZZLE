function g = gamma_n(sp,T)
% from SSH theory
global h c w wx I sw_o m k

%% TEST GAMMMA_N
%global h c w wx I sw_o m

h_bar = h/2/pi;

i = 0 : I(sw_o,sp) - 1;

mu = m(sp).*m./(m(sp) + m);
r0 = [3.621 3.458 3.47 3.298 2.75].*1e-10; % m
r0 = diag(r0);

for y = 1 : 5
    for j = 1 : 5
        r0(y,j) = 0.5*(r0(y,y) + r0(j,j));
    end
end

alpha = 17.5./r0(sp,:);

e = h*c.*(w(sp).*i - wx(sp).*i - wx(sp).*i.^2); 

g = zeros(I(sw_o,sp),5);

for s = j : length(i) - 1
    g(s,:) = pi*(e(s + 1) - e(s))./h_bar./alpha.*sqrt(mu.*0.5/k/T);
end

end
% â main
% sp = 3;
% 
% T = 500 : 100 : 10000;
% g_n = zeros(I(sw_o,sp),5,length(T));
% 
% for i = 1 : length(T)
%     g_n(:,:,i) = gamma_n(sp,T(i));
% end