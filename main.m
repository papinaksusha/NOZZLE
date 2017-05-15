%% 5-components air mixture N2/O2/NO/N/O, state-to-state kinetics in nozzle  

clc
clear 
close all
format long e
addpath('./MAT/')

global Na k h c w wx m theta_r D I sw_o sw_n ex_model e_STELLAR molar

%% CONSTs
% J = kg*m^2/s^2, Pa = kg/m/s^2
% [N2 O2 NO N O]

tic

Na = 6.022141e23;                                                          % mol^-1, Avogadro constant 
k = 1.380648e-23;                                                          % J/K, Boltzmann constant
h = 6.626070e-34;                                                          % J*s, Planck constant
c = 299792458;                                                             % m/s, speed of light 
w = [235857 158019 190420];                                                % m^-1, spectroscopic constant
sw_o = 3;                                                                  % switch on oscillator, 1 -  harmonic oscillator; 2 -  anharmonic oscillator 3 - STELLAR
switch sw_o 
    case 1
        wx = [0 0 0];   
    case 2
        wx = [1432 1198 1407.5];                                           % m^-1, spectroscopic constant
    case 3
        wx = [1432 1198 1407.5];   
end
molar = 1e-3.*[28.0134 31.99880 30.00610 14.0067 15.9994];                 % kg/mol, molar mass         
m = 1.6605402e-24.*molar;                                                  % kg, molecular mass
theta_r = [2.86 2.07 2.42];                                                % K, characteristic rotational emperature
D = h*c.*[7.871e6 4.126e6 5.24e6];                                         % J, dissociation energy
sw_n = 1;                                                                  % switcn on nozzle
ex_model = 3;                                                              % exchange rate coefficients model
                                                       
% Max vibrational levels                                                 
% I(1,:) - harmonic oscillator, I(2,:) - anharmonic oscillator, I(3,:) - STELLAR distributions

I = [33 26 27; 47 36 39; 60, 45, 47]; 

% STELLAR vibrational distributions

e_STELLAR = cell(3,1);

fileID1 = fopen('./MAT/STELLAR/N2-X-v61-Levels-Leroy2006.txt');
fileID2 = fopen('./MAT/STELLAR/O2-X-v46-Levels-RKR.txt');
fileID3 = fopen('./MAT/STELLAR/NO-X-v48-Levels-RKR.txt');

a = textscan(fileID1,'%f');
a = h.*c.*1e2.*cell2mat(a)';
e_STELLAR(1) = {a - a(1)};
a = textscan(fileID2,'%f');
a = h.*c.*1e2.*cell2mat(a)';
e_STELLAR(2) = {a - a(1)};
a = textscan(fileID3,'%f');
a = h.*c.*1e2.*cell2mat(a)';
e_STELLAR(3) = {a - a(1)};
fclose('all');

e_0_N2_STELLAR = h*c*1e-2*1175.78;
e_0_O2_STELLAR = h*c*1e-2*787.38;
e_0_NO_STELLAR = h*c*1e-2*948.646642922263;

%%

fig = 1;

% N_T = 100;
% TT = 1000 : N_T : 9000;

% TEST_SSH;
% TEST_SSH_ALEX
% TEST_DISS_TM;
% TEST_EXCHANGE


% callSTELLAR;
% callSTELLARwithoutNO;
% call5;
% call2;
