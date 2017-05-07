function [fig] = TEST_EXCHANGE(fig,N_T,T)

global sw_o I h k m D theta_r c w wx

k_ex_N2 = zeros(I(sw_o,1) + 1,4,length(T));
k_ex_O2 = zeros(I(sw_o,2) + 1,4,length(T));

k_ex_N = zeros(I(sw_o,1) + 1,4,length(T));
k_ex_O = zeros(I(sw_o,2) + 1,4,length(T));

i = {0 : I(sw_o,1), 0 : I(sw_o,2)};

e_N2 = h*c.*(w(1).*cell2mat(i(1)) - wx(1).*cell2mat(i(1)) - ...
            wx(1).*cell2mat(i(1)).^2);
e_O2 = h*c.*(w(2).*cell2mat(i(2)) - wx(2).*cell2mat(i(2)) - ...
            wx(2).*cell2mat(i(2)).^2);

for j = 1 : length(T);
    k_ex_N2(:,:,j) = k_ex(1,T(j));
    k_ex_O2(:,:,j) = k_ex(2,T(j));
% œ–Œ¬≈–»“‹!
    k_ex_N(:,:,j) = k_ex_N2(:,:,j).*(m(1)*m(5)/m(3)/m(4))^(1.5).*...
                    theta_r(3)/theta_r(1)*0.5.*exp((-e_N2'*ones(1,4) + D(1) - D(3))./k./T(j));
    k_ex_O(:,:,j) = k_ex_O2(:,:,j).*(m(2)*m(4)/m(3)/m(5))^(1.5).*...
                     theta_r(3)/theta_r(2)*0.5.*exp((-e_O2'*ones(1,4) + D(2) - D(3))./k./T(j));
end

%% N2 + O from T
figure(fig)
plot(T, squeeze(log10(k_ex_N2(1,:,:))))
title('i = 0');
ylabel('lg(k_{N2i,NO}^{O,N})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_N(1,:,:))))
title('i = 0');
ylabel('lg(k_{NO,N2i}^{N,O})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_N2(11,:,:))))
title('i = 10');
ylabel('lg(k_{N2i,NO}^{O,N})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_N(11,:,:))))
title('i = 10');
ylabel('lg(k_{NO,N2i}^{N,O})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_N2(21,:,:))))
title('i = 20');
ylabel('lg(k_{N2i,NO}^{O,N})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_N(21,:,:))))
title('i = 20');
ylabel('lg(k_{NO,N2i}^{N,O})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

%% O2 + N from T

figure(fig)
plot(T, squeeze(log10(k_ex_O2(1,:,:))))
title('i = 0');
ylabel('lg(k_{O2i,NO}^{N,O})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_O(1,:,:))))
title('i = 0');
ylabel('lg(k_{NO,O2i}^{O,N})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_O2(11,:,:))))
title('i = 10');
ylabel('lg(k_{O2i,NO}^{N,O})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_O(11,:,:))))
title('i = 10');
ylabel('lg(k_{NO,O2i}^{O,N})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_O2(16,:,:))))
title('i = 15');
ylabel('lg(k_{O2i,NO}^{N,O})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(T, squeeze(log10(k_ex_N(16,:,:))))
title('i = 15');
ylabel('lg(k_{NO,O2i}^{O,N})');
xlabel('T');
legend('1','2','3','4');
fig = fig + 1;


%% N2 + O from i
figure(fig)
plot(cell2mat(i(1)), squeeze(log10(k_ex_N2(:,:,(5000-T(1))/N_T+1))))
title('T = 5000 K');
ylabel('lg(k_{N2i,NO}^{O,N})');
xlabel('i');
xlim([0, I(sw_o,1) - 1]);
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(1)), squeeze(log10(k_ex_N(:,:,(5000-T(1))/N_T+1))))
title('T = 5000 K');
ylabel('lg(k_{NO,N2i}^{N,O})');
xlabel('i');
xlim([0, I(sw_o,1) - 1]);
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(1)), squeeze(log10(k_ex_N2(:,:,(8000-T(1))/N_T+1))))
title('T = 8000 K');
ylabel('lg(k_{N2i,NO}^{O,N})');
xlabel('i');
xlim([0, I(sw_o,1) - 1]);
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(1)), squeeze(log10(k_ex_N(:,:,(8000-T(1))/N_T+1))))
title('T = 8000 K');
ylabel('lg(k_{NO,N2i}^{N,O})');
xlabel('i');
xlim([0, I(sw_o,1) - 1]);
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(1)), squeeze(log10(k_ex_N2(:,:,(14000-T(1))/N_T+1))))
title('T = 14000 K');
ylabel('lg(k_{N2i,NO}^{O,N})');
xlabel('i');
xlim([0, I(sw_o,1) - 1]);
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(1)), squeeze(log10(k_ex_N(:,:,(14000-T(1))/N_T+1))))
title('T = 14000 K');
ylabel('lg(k_{NO,N2i}^{N,O})');
xlabel('i');
xlim([0, I(sw_o,1) - 1]);
legend('1','2','3','4');
fig = fig + 1;

%% O2 + N from T

figure(fig)
plot(cell2mat(i(2)), squeeze(log10(k_ex_O2(:,:,(5000-T(1))/N_T+1))))
title('T = 5000 K');
ylabel('lg(k_{O2i,NO}^{N,O})');
xlabel('i');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(2)), squeeze(log10(k_ex_O(:,:,(5000-T(1))/N_T+1))))
title('T = 5000 K');
ylabel('lg(k_{NO,O2i}^{O,N})');
xlabel('i');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(2)), squeeze(log10(k_ex_O2(:,:,(8000-T(1))/N_T+1))))
title('T = 8000 K');
ylabel('lg(k_{O2i,NO}^{N,O})');
xlabel('i');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(2)), squeeze(log10(k_ex_O(:,:,(8000-T(1))/N_T+1))))
title('T = 8000 K');
ylabel('lg(k_{NO,O2i}^{O,N})');
xlabel('i');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(2)), squeeze(log10(k_ex_O2(:,:,(14000-T(1))/N_T+1))))
title('T = 14000 K');
ylabel('lg(k_{O2i,NO}^{N,O})');
xlabel('i');
legend('1','2','3','4');
fig = fig + 1;

figure(fig)
plot(cell2mat(i(2)), squeeze(log10(k_ex_O(:,:,(14000-T(1))/N_T+1))))
title('T = 14000 K');
ylabel('lg(k_{NO,O2i}^{O,N})');
xlabel('i');
legend('1','2','3','4');
fig = fig + 1;

end

