%% MVM course
close all; clc; clear; 
warning off

%% ------------------------- Initialization -------------------------

% Simulink model name
simulink_model_name = 'Tyre_Pac96';
% Load tyre data (Pacejka model coefficients)
WheelFile = 'Tyre215_50_19_Comb';      
% Wheel Initialization
pacn = [];
eval(['[Pacejka]=' WheelFile ';'])
pacn = struct2cell(Pacejka);
for ii = 1:size(pacn)    Pace(ii) = pacn{ii}; end
% coefficients of Pacejka tyre model
Pacn = Pace';     


% default coefficient setting
gamma0 = 0;     gamma_slope = 0;
s0 = -1;        s_slope = 0.2;      tau_slip = 0.1; % (tau = Lrel / V)
alfa0 = 0;      alfa_slope = 0.2;   tau_alfa = 0.1;
Fz0 = 1500;     Fz_slope = 0;
mu_slope = 0;   mu0 = 1;
V = 100/3.6; % 100 Km/h
Wheel_Radius = 0.3;

% Plot settings 
F_Size = 14; % FontSize
plotcol = {'--k','-r','-.b',':g','.m','-k'};
