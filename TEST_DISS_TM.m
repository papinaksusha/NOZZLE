function [fig] = TEST_DISS_TM(fig,N_T,T)

global sw_o I h k m D theta_r c w wx

k_diss_N2 = zeros(I(sw_o,1) + 1,length(T));
k_diss_O2 = zeros(I(sw_o,2) + 1,length(T));

i_N2 = 0 : I(sw_o,1);
i_O2 = 0 : I(sw_o,2);

for j = 1 : length(T)   
    k_diss_N2(:,j) = k_diss(1,T(j)); 
    k_diss_O2(:,j) = k_diss(2,T(j));
end

%% DISS

figure(fig)
hold on
plot(T, log10(k_diss_N2(1,:)))
plot(T, log10(k_diss_N2(16,:)))
plot(T, log10(k_diss_N2(31,:)))
plot(T, log10(k_diss_N2(48,:)))
title('N2');
xlabel('T');
ylabel('lg(k_{N_2i,diss})')
legend('i=0','i=15','i=30','i=47');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(T, log10(k_diss_O2(1,:)))
plot(T, log10(k_diss_O2(16,:)))
plot(T, log10(k_diss_O2(31,:)))
plot(T, log10(k_diss_O2(37,:)))
xlabel('T');
ylabel('lg(k_{O_2i,diss})');
legend('i=0','i=15','i=30','i=36');
title('O2')
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_N2, log10(k_diss_N2(:,(5000-T(1))/N_T+1)))
plot(i_N2, log10(k_diss_N2(:,(8000-T(1))/N_T+1)))
plot(i_N2, log10(k_diss_N2(:,(14000-T(1))/N_T+1)))
title('N2');
xlabel('i');
xlim([0,I(sw_o,1)]);
ylabel('lg(k_{N_2i,diss})')
legend('T = 5000 K','T = 8000 K','T = 14000 K');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_O2, log10(k_diss_O2(:,(5000-T(1))/N_T+1)))
plot(i_O2, log10(k_diss_O2(:,(8000-T(1))/N_T+1)))
plot(i_O2, log10(k_diss_O2(:,(14000-T(1))/N_T+1)))
title('O2');
xlabel('i');
xlim([0,I(sw_o,2)]);
ylabel('lg(k_{O_2i,diss})')
legend('T = 5000 K','T = 8000 K','T = 14000 K');
hold off
fig = fig + 1;

%% REC 

% Detailed balance principle

e_N2 = h*c.*(w(1).*i_N2 - wx(1).*i_N2 - wx(1).*i_N2.^2);
e_O2 = h*c.*(w(2).*i_O2 - wx(2).*i_O2 - wx(2).*i_O2.^2);

k_rec_N2 = zeros(I(sw_o,1) + 1,length(T));
k_rec_O2 = zeros(I(sw_o,2) + 1,length(T));

for j = 1 : length(T)
    k_rec_N2(:,j) = k_diss_N2(:,j).*(m(1)/m(4)^2)^(1.5).*h^3.*(2*pi*k*T(j))^(-1.5).*...
           T(j)./theta_r(1).*0.5.*exp((-e_N2' + D(1))/k/T(j));
    k_rec_O2(:,j) = k_diss_O2(:,j).*(m(2)/m(5)^2)^(1.5).*h^3.*(2*pi*k*T(j))^(-1.5).*...
           T(j)./theta_r(2).*0.5.*exp((-e_O2' + D(2))/k/T(j));
end

figure(fig)
hold on
plot(T, log10(k_rec_N2(1,:)))
plot(T, log10(k_rec_N2(16,:)))
plot(T, log10(k_rec_N2(31,:)))
plot(T, log10(k_rec_N2(48,:)))
title('N2');
xlabel('T');
ylabel('lg(k_{rec, N_2i})')
legend('i=0','i=15','i=30','i=47');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(T, log10(k_rec_O2(1,:)))
plot(T, log10(k_rec_O2(16,:)))
plot(T, log10(k_rec_O2(31,:)))
plot(T, log10(k_rec_O2(37,:)))
xlabel('T');
ylabel('lg(k_{rec, O_2i})');
legend('i=0','i=15','i=30','i=36');
title('O2')
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_N2, log10(k_rec_N2(:,(5000-T(1))/N_T+1)))
plot(i_N2, log10(k_rec_N2(:,(8000-T(1))/N_T+1)))
plot(i_N2, log10(k_rec_N2(:,(14000-T(1))/N_T+1)))
title('N2');
xlabel('i');
xlim([0,I(sw_o,1)]);
ylabel('lg(k_{rec, N_2i})')
legend('T = 5000 K','T = 8000 K','T = 14000 K');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_O2, log10(k_rec_O2(:,(5000-T(1))/N_T+1)))
plot(i_O2, log10(k_rec_O2(:,(8000-T(1))/N_T+1)))
plot(i_O2, log10(k_rec_O2(:,(14000-T(1))/N_T+1)))
title('O2');
xlabel('i');
xlim([0,I(sw_o,2)]);
ylabel('lg(k_{rec, O_2i})')
legend('T = 5000 K','T = 8000 K','T = 14000 K');
hold off
fig = fig + 1;

end

