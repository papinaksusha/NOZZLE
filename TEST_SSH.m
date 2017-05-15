%% TEST SSH 
% compare with Olya's dissertation

global sw_o

SW_O = sw_o;

sw_o = 2;

N_T = 100;
TT = 2000 : N_T : 14000;

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

%% VT 

f = squeeze(K_ssh_VT_N2(1,:,:));
figure(fig)
plot(TT, log10(f));
title('VT N2 1->0');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{N2,1->0}^{M})')
fig = fig + 1;

figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_O2(1,:,:))));
title('VT O2 1->0');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{O2,1->0}^{M})')
fig = fig + 1;

figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_N2(21,:,:))));
title('VT N2 21->20');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{N2,21->20}^{M})')
fig = fig + 1;

figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_O2(21,:,:))));
title('VT O2 21->20');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{O2,21->20}^{M})')
fig = fig + 1;

figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_N2(47,:,:))));
title('VT N2 47->46');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{N2,47->46}^{M})')
fig = fig + 1;

figure(fig)
plot(TT, log10(squeeze(K_ssh_VT_O2(36,:,:))));
title('VT O2 36->35');
legend('N2','O2','NO','N','O')
xlabel('T');
ylabel('lg(k_{O2,36->35}^{M})')
fig = fig + 1;

figure(fig)
hold on
subplot(2,2,1) 
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VT_N2(:,:,(5000-TT(1))/N_T+1))))
title('VT N2 T = 5000 K');
xlabel('i');
ylabel('lg(k_{N_2,i->i-1}^{M})');
xlim([0,I(sw_o,1) - 1]);
subplot(2,2,2) 
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VT_O2(:,:,(5000-TT(1))/N_T+1))))
title('VT O2 T = 5000 K');
xlabel('i');
ylabel('lg(k_{O_2,i->i-1}^{M})');
xlim([0,I(sw_o,2) - 1]);
subplot(2,2,3) 
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VT_N2(:,:,(8000-TT(1))/N_T+1))))
title('VT N2 T = 8000 K');
xlabel('i');
xlim([0,I(sw_o,1) - 1]);
ylabel('lg(k_{N_2,i->i-1}^{M})');
subplot(2,2,4) 
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VT_O2(:,:,(8000-TT(1))/N_T+1))))
title('VT O2 T = 8000 K');
xlim([0,I(sw_o,2) - 1]);
xlabel('i');
ylabel('lg(k_{O_2,i->i-1}^{M})');
hold off
fig = fig + 1;

%% VV

figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(1,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(11,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(21,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_N2(31,1,:))))
title('VV N2-N2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{N_2,i->i-1}^{N_2,0->1})');
hold off
fig = fig +1;

figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(1,36,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(11,36,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(21,36,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_O2(31,36,:))))
title('VV O2-O2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{O_2,i->i-1}^{O_2,35->36})');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,1,(5000-TT(1))/N_T+1))))
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,1,(8000-TT(1))/N_T+1))))
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_N2(:,1,(14000-TT(1))/N_T+1))))
title('VV N2-N2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(sw_o,1)]);
ylabel('lg(k_{N_2,i->i-1}^{N_2,0->1})');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VV_O2_O2(:,36,(5000-TT(1))/N_T+1))))
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VV_O2_O2(:,36,(8000-TT(1))/N_T+1))))
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VV_O2_O2(:,36,(14000-TT(1))/N_T+1))))
title('VV O2-O2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(sw_o,2)]);
ylabel('lg(k_{O_2,i->i-1}^{O_2,35->36})');
hold off
fig = fig + 1;

% VV'

figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(1,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(11,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(21,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_N2_O2(31,1,:))))
title('VV" N2-O2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{N_2,i->i-1}^{O_2,0->1})');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(1,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(11,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(21,1,:))))
plot(TT, log10(squeeze(K_ssh_VV_O2_N2(31,1,:))))
title('VV" O2-N2');
legend('i = 1','i = 11','i = 21','i = 31')
xlabel('T');
ylabel('lg(k_{O_2,i->i-1}^{N_2,35->36})');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_O2(:,1,(5000-TT(1))/N_T+1))))
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_O2(:,1,(8000-TT(1))/N_T+1))))
plot(0 : I(sw_o,1) - 1, log10(squeeze(K_ssh_VV_N2_O2(:,1,(14000-TT(1))/N_T+1))))
title('VV N2-O2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(sw_o,1)]);
ylabel('lg(k_{N_2,i->i-1}^{O_2,0->1})');
hold off
fig = fig + 1;

figure(fig)
hold on
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VV_O2_N2(:,1,(5000-TT(1))/N_T+1))))
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VV_O2_N2(:,1,(8000-TT(1))/N_T+1))))
plot(0 : I(sw_o,2) - 1, log10(squeeze(K_ssh_VV_O2_N2(:,1,(14000-TT(1))/N_T+1))))
title('VV" O2-N2');
legend('T = 5000','T = 8000','T = 14000')
xlabel('i');
xlim([0, I(sw_o,2)]);
ylabel('lg(k_{O_2,i->i-1}^{N_2,35->36})');
hold off
fig = fig + 1;
%%
figure(fig)

semilogy(TT, squeeze(K_ssh_VT_N2(9,1,:))), hold on
semilogy(TT, squeeze(K_ssh_VT_N2(9,4,:)))
semilogy(TT, squeeze(K_ssh_VV_N2_N2(9,3,:)))
legend('VT whith molecul','VT with atom','VV');
hold off
fig = fig + 1;

sw_o = SW_O;