%% TEST SSH 
% compare with Alexandrova's dissertation

global sw_o

SW_O = sw_o;
sw_o = 2;

N_T = 100;
TT = 0 : N_T : 9000;


K_ssh_VT_N2  = zeros(I(sw_o,1),5,length(TT));
K_ssh_VV_N2_N2 = zeros(I(sw_o,1),I(sw_o,1),length(TT));
K_ssh_VV_N2_O2 = zeros(I(sw_o,1),I(sw_o,2),length(TT));
k_ssh_VV_N2_NO = zeros(I(sw_o,1),I(sw_o,3),length(TT));

K_ssh_VT_O2  = zeros(I(sw_o,2),5,length(TT));
K_ssh_VV_O2_N2 = zeros(I(sw_o,2),I(sw_o,1),length(TT));
K_ssh_VV_O2_O2 = zeros(I(sw_o,2),I(sw_o,2),length(TT));
k_ssh_VV_O2_NO = zeros(I(sw_o,2),I(sw_o,3),length(TT));


for  i = 1 : length(TT)
   [K_ssh_VT_N2(:,:,i), K_ssh_VV_N2_N2(:,:,i), K_ssh_VV_N2_O2(:,:,i), k_ssh_VV_N2_NO(:,:,i)] = k_ssh(1,TT(i));
   [K_ssh_VT_O2(:,:,i), K_ssh_VV_O2_N2(:,:,i), K_ssh_VV_O2_O2(:,:,i), k_ssh_VV_O2_NO(:,:,i)] = k_ssh(2,TT(i));
end

figure(fig)
semilogy(TT, squeeze(K_ssh_VT_N2(1,1,:))), hold on
semilogy(TT, squeeze(K_ssh_VT_N2(1,4,:)))
legend('mol','at')
xlabel('T')
ylabel('k_{N_2, 1->0}^{M}')
ylim([1e-24,1e-16]);
fig = fig + 1;

figure(fig)
semilogy(TT, squeeze(K_ssh_VV_N2_N2(1,1,:)))
xlabel('T')
ylabel('k_{N_2, 1->0}^{0->1}')
fig = fig + 1;

figure(fig)
semilogy(TT, squeeze(K_ssh_VT_N2(45,1,:))), hold on
semilogy(TT, squeeze(K_ssh_VT_N2(45,4,:)))
legend('mol','at')
xlabel('T')
ylabel('k_{N_2, 45->44}^{M}')
fig = fig + 1;
ylim([1e-20,1e-14]);

figure(fig)
semilogy(TT, squeeze(K_ssh_VV_N2_N2(45,1,:)))
xlabel('T')
ylabel('k_{N_2, 1->0}^{0->1}')
fig = fig + 1;
ylim([1e-20,1e-15]);

figure(fig)
semilogy(TT, squeeze(K_ssh_VT_N2(16,1,:))), hold on
semilogy(TT, squeeze(K_ssh_VT_N2(16,4,:)))
legend('mol','at')
xlabel('T')
ylabel('k_{N_2, 16->15}^{M}')
ylim([1e-21,1e-15]);
fig = fig + 1;

figure(fig)
semilogy(0 : I(sw_o,1) - 1, squeeze(K_ssh_VV_N2_N2(:,1,(4000-TT(1))/N_T+1))), hold on
semilogy(0 : I(sw_o,1) - 1, squeeze(K_ssh_VT_N2(:,1,(4000-TT(1))/N_T+1)))
semilogy(0 : I(sw_o,1) - 1, squeeze(K_ssh_VT_N2(:,4,(4000-TT(1))/N_T+1)))
legend('VV','VT mol','VT at')
xlabel('i');
xlim([0, I(sw_o,1)]);
%ylabel('lg(k_{N_2,i->i-1}^{N_2,0->1})');
hold off
fig = fig + 1;


figure(fig)
semilogy(0 : I(sw_o,1) - 1, squeeze(K_ssh_VV_N2_N2(:,1,(6000-TT(1))/N_T+1))), hold on
semilogy(0 : I(sw_o,1) - 1, squeeze(K_ssh_VT_N2(:,1,(6000-TT(1))/N_T+1)))
semilogy(0 : I(sw_o,1) - 1, squeeze(K_ssh_VT_N2(:,4,(6000-TT(1))/N_T+1)))
legend('VV','VT mol','VT at')
xlabel('i');
xlim([0, I(sw_o,1)]);
%ylabel('lg(k_{N_2,i->i-1}^{N_2,0->1})');
hold off
fig = fig + 1;


figure(fig)
semilogy(TT, squeeze(K_ssh_VV_O2_O2(1,1,:))), hold on
semilogy(TT, squeeze(K_ssh_VT_O2(1,2,:)))
semilogy(TT, squeeze(K_ssh_VT_O2(1,5,:)))
legend('VV','VT mol','VT at')
xlabel('T')
xlim([0, 7000]);
ylabel('k_{O_2, 1->0}^{M}, k_{O_2,1->0}^{O_2,0->1}')
ylim([1e-22,1e-16]);
fig = fig + 1;

figure(fig)
semilogy(TT, squeeze(K_ssh_VV_O2_O2(16,1,:))), hold on
semilogy(TT, squeeze(K_ssh_VT_O2(16,1,:)))
semilogy(TT, squeeze(K_ssh_VT_O2(16,4,:)))
legend('VV','VT mol','VT at')
xlabel('T')
xlim([0, 7000]);
ylabel('k_{O_2, 16->15}^{M}, k_{O_2,16->15}^{O_2,0->1}')
ylim([1e-20,1e-14]);
fig = fig + 1;

figure(fig)
semilogy(0 : I(sw_o,2) - 1, squeeze(K_ssh_VV_O2_O2(:,1,(4000-TT(1))/N_T+1))), hold on
semilogy(0 : I(sw_o,2) - 1, squeeze(K_ssh_VT_O2(:,1,(4000-TT(1))/N_T+1)))
semilogy(0 : I(sw_o,2) - 1, squeeze(K_ssh_VT_O2(:,4,(4000-TT(1))/N_T+1)))
legend('VV','VT mol','VT at')
xlabel('i');
xlim([0, I(sw_o,2)]);
ylabel('k_{O_2,i->i-1}^{M}, k_{O_2,i->i-1}^{O_2,0->1}');
hold off
fig = fig + 1;


figure(fig)
semilogy(0 : I(sw_o,2) - 1, squeeze(K_ssh_VV_O2_O2(:,1,(2000-TT(1))/N_T+1))), hold on
semilogy(0 : I(sw_o,2) - 1, squeeze(K_ssh_VT_O2(:,1,(2000-TT(1))/N_T+1)))
semilogy(0 : I(sw_o,2) - 1, squeeze(K_ssh_VT_O2(:,4,(2000-TT(1))/N_T+1)))
legend('VV','VT mol','VT at')
xlabel('i');
xlim([0, I(sw_o,2)]);
ylabel('k_{O_2,i->i-1}^{M}, k_{O_2,i->i-1}^{O_2,0->1}');
hold off
fig = fig + 1;

figure(fig)
semilogy(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,17,(2000-TT(1))/N_T+1)))),hold on
semilogy(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,17,(4000-TT(1))/N_T+1))))
semilogy(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,17,(6000-TT(1))/N_T+1))))
legend('T = 2000','T = 4000','T = 6000')
xlabel('i');
xlim([0, I(sw_o,2)]);
ylabel('k_{N_2,i->i-1}^{N_2,16->17}');
hold off
fig = fig + 1;

sw_o = SW_O;