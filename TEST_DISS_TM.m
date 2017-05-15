%% TEST DISSOSIATION RATE COEFFICIENTS, TRINOR-MARRONE MODEL
% compare with Olya's dissertation 

global sw_o

SW_O = sw_o;

sw_o = 2;

N_T = 100;
TT = 2000 : N_T : 14000;

k_diss_N2 = zeros(I(sw_o,1) + 1,5,length(TT));
k_diss_O2 = zeros(I(sw_o,2) + 1,5,length(TT));

i_N2 = 0 : I(sw_o,1);
i_O2 = 0 : I(sw_o,2);

for j = 1 : length(TT)   
    k_diss_N2(:,:,j) = k_diss(1,TT(j)); 
    k_diss_O2(:,:,j) = k_diss(2,TT(j));
end

%% DISS

figure(fig)
hold on
plot(TT, squeeze(log10(k_diss_N2(1,1,:)))) 
plot(TT, squeeze(log10(k_diss_N2(16,1,:)))) 
plot(TT, squeeze(log10(k_diss_N2(31,1,:)))) 
plot(TT, squeeze(log10(k_diss_N2(48,1,:)))) 

plot(TT,  squeeze(log10(k_diss_N2(1,4,:))), '--')
plot(TT,  squeeze(log10(k_diss_N2(16,4,:))), '--')
plot(TT,  squeeze(log10(k_diss_N2(31,4,:))), '--')
plot(TT,  squeeze(log10(k_diss_N2(48,4,:))), '--')
title('N2');
xlabel('T');
ylabel('lg(k_{N_2i,diss})')
legend('i=0','i=15','i=30','i=47');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(TT, squeeze(log10(k_diss_O2(1,1,:))))
plot(TT, squeeze(log10(k_diss_O2(16,1,:))))
plot(TT, squeeze(log10(k_diss_O2(31,1,:))))
plot(TT, squeeze(log10(k_diss_O2(37,1,:))))
plot(TT, squeeze(log10(k_diss_O2(1,4,:))),'--')
plot(TT, squeeze(log10(k_diss_O2(16,4,:))),'--')
plot(TT, squeeze(log10(k_diss_O2(31,4,:))),'--')
plot(TT, squeeze(log10(k_diss_O2(37,4,:))),'--')
xlabel('T');
ylabel('lg(k_{O_2i,diss})');
legend('i=0','i=15','i=30','i=36');
title('O2')
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_N2, squeeze(log10(k_diss_N2(:,1,(5000-TT(1))/N_T+1))))
plot(i_N2, squeeze(log10(k_diss_N2(:,1,(8000-TT(1))/N_T+1))))
plot(i_N2, squeeze(log10(k_diss_N2(:,1,(14000-TT(1))/N_T+1))))
plot(i_N2, squeeze(log10(k_diss_N2(:,4,(5000-TT(1))/N_T+1))),'--')
plot(i_N2, squeeze(log10(k_diss_N2(:,4,(8000-TT(1))/N_T+1))),'--')
plot(i_N2, squeeze(log10(k_diss_N2(:,4,(14000-TT(1))/N_T+1))),'--')
title('N2');
xlabel('i');
xlim([0,I(sw_o,1)]);
ylabel('lg(k_{N_2i,diss})')
legend('T = 5000 K','T = 8000 K','T = 14000 K');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_O2, squeeze(log10(k_diss_O2(:,1,(5000-TT(1))/N_T+1))))
plot(i_O2, squeeze(log10(k_diss_O2(:,1,(8000-TT(1))/N_T+1))))
plot(i_O2, squeeze(log10(k_diss_O2(:,1,(14000-TT(1))/N_T+1))))
plot(i_O2, squeeze(log10(k_diss_O2(:,4,(5000-TT(1))/N_T+1))),'--')
plot(i_O2, squeeze(log10(k_diss_O2(:,4,(8000-TT(1))/N_T+1))),'--')
plot(i_O2, squeeze(log10(k_diss_O2(:,4,(14000-TT(1))/N_T+1))),'--')
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
 
 k_rec_N2 = zeros(I(sw_o,1) + 1,5,length(TT));
 k_rec_O2 = zeros(I(sw_o,2) + 1,5,length(TT));

for j = 1 : length(TT)
    k_rec_N2(:,:,j) = k_diss_N2(:,:,j).*(m(1)/m(4)^2)^(1.5).*h^3.*(2*pi*k*TT(j))^(-1.5).*...
           TT(j)./theta_r(1).*0.5.*exp((-e_N2'*ones(1,5) + D(1))/k/TT(j));
    k_rec_O2(:,:,j) = k_diss_O2(:,:,j).*(m(2)/m(5)^2)^(1.5).*h^3.*(2*pi*k*TT(j))^(-1.5).*...
           TT(j)./theta_r(2).*0.5.*exp((-e_O2'*ones(1,5) + D(2))/k/TT(j));
end

figure(fig)
hold on
plot(TT, squeeze(log10(k_rec_N2(1,1,:))))
plot(TT, squeeze(log10(k_rec_N2(16,1,:))))
plot(TT, squeeze(log10(k_rec_N2(31,1,:))))
plot(TT, squeeze(log10(k_rec_N2(48,1,:))))
plot(TT, squeeze(log10(k_rec_N2(1,4,:))),'--')
plot(TT, squeeze(log10(k_rec_N2(16,4,:))),'--')
plot(TT, squeeze(log10(k_rec_N2(31,4,:))),'--')
plot(TT, squeeze(log10(k_rec_N2(48,4,:))),'--')
title('N2');
xlabel('T');
ylabel('lg(k_{rec, N_2i})')
legend('i=0','i=15','i=30','i=47');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(TT, squeeze(log10(k_rec_O2(1,1,:))))
plot(TT, squeeze(log10(k_rec_O2(16,1,:))))
plot(TT, squeeze(log10(k_rec_O2(31,1,:))))
plot(TT, squeeze(log10(k_rec_O2(37,1,:))))
plot(TT, squeeze(log10(k_rec_O2(1,4,:))),'--')
plot(TT, squeeze(log10(k_rec_O2(16,4,:))),'--')
plot(TT, squeeze(log10(k_rec_O2(31,4,:))),'--')
plot(TT, squeeze(log10(k_rec_O2(37,4,:))),'--')
xlabel('T');
ylabel('lg(k_{rec, O_2i})');
legend('i=0','i=15','i=30','i=36');
title('O2')
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_N2, squeeze(log10(k_rec_N2(:,1,(5000-TT(1))/N_T+1))))
plot(i_N2, squeeze(log10(k_rec_N2(:,1,(8000-TT(1))/N_T+1))))
plot(i_N2, squeeze(log10(k_rec_N2(:,1,(14000-TT(1))/N_T+1))))
plot(i_N2, squeeze(log10(k_rec_N2(:,4,(5000-TT(1))/N_T+1))),'--')
plot(i_N2, squeeze(log10(k_rec_N2(:,4,(8000-TT(1))/N_T+1))),'--')
plot(i_N2, squeeze(log10(k_rec_N2(:,4,(14000-TT(1))/N_T+1))),'--')
title('N2');
xlabel('i');
xlim([0,I(sw_o,1)]);
ylabel('lg(k_{rec, N_2i})')
legend('T = 5000 K','T = 8000 K','T = 14000 K');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(i_O2, squeeze(log10(k_rec_O2(:,1,(5000-TT(1))/N_T+1))))
plot(i_O2, squeeze(log10(k_rec_O2(:,1,(8000-TT(1))/N_T+1))))
plot(i_O2, squeeze(log10(k_rec_O2(:,1,(14000-TT(1))/N_T+1))))
plot(i_O2, squeeze(log10(k_rec_O2(:,4,(5000-TT(1))/N_T+1))),'--')
plot(i_O2, squeeze(log10(k_rec_O2(:,4,(8000-TT(1))/N_T+1))),'--')
plot(i_O2, squeeze(log10(k_rec_O2(:,4,(14000-TT(1))/N_T+1))),'--')
title('O2');
xlabel('i');
xlim([0,I(sw_o,2)]);
ylabel('lg(k_{rec, O_2i})')
legend('T = 5000 K','T = 8000 K','T = 14000 K');
hold off
fig = fig + 1;

sw_o = SW_O;