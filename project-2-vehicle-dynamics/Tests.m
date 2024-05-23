close all; clc; clear; 
warning off

%% ------------------------- Initialization ------------------------- %%
simulink_model_name = 'Tyre_Pac96';
Time_sim = 10; % Simulation time

% Load tyre data (Pacejka model coefficients)
WheelFile = 'Tyre215_50_19_Comb';
eval(['[Pacejka]=' WheelFile ';'])
pacn = struct2cell(Pacejka);
for ii = 1:size(pacn) Pace(ii)=pacn{ii}; end
Pacn = Pace';


% ------ vehicle parameters ------ %
m = 1812; % kerb weight [Kg]
g = 9.81; %gravità
wheelbase = 2.77; %  [m]
cg_height = 0.55; % center of gravity height [m]
% TODO mancano valori dimensione car
% ...
Fzf_ratio= 0.5;
Fzr_ratio = 0.5;

% wheel
wheel_radius = 0.348; % [m] 
Ir = 1.46; % wheel inertia [Kg*m^2] % 0.5*(11+13)*(0.348)^2

% motor
peak_power = 150000; % [W]
max_torque = 310; % [Nm]
max_motor_speed = 16000 * 2*pi/60; % [rpm] -> [rad/sec]
min_motor_speed = peak_power / max_torque; % [rad/sec] basespeed
gear_ratio = 10.5;
motor_eff = 0.9;
transm_eff = 0.95;
torsional_stiffness = 9000; % OPTIONAL!!

% motor delay
motor_delay = 0.02; % [s]
motor_risetime = 0.05; % [s]

% battery
battery_cap = 58000; % [Wh]
battery_volt = 800; % [V]

% inverter
inverter_eff = 0.9;

% brakes
brake_delay = 0.02; % [s]
brake_risetime = 0.025; % [s]
brake_front = 0.25; % brake distribution 75:25 
brake_rear = 0.75;

% rolling resistance
f0 = 0.009; 
f2 = 6.5e-6; % [s^2/m^2]

% air drag
rho = 1.225; % [kg/m³] air density 
Af = 2.36; % [m^2] frontal area
Cx = 0.27; % drag coefficient


% ------ signals coefficients setting ------ %
% camber angle
gamma0 = 0;     
gamma_slope = 0;

% slip
s0 = 0.1;         
s_slope = 0;      
tau_slip = 0.1; % (tau = Lrel / V)

% slip angle
alfa0 = 0;      
alfa_slope = 0;   
tau_alfa = 0.1;

% vertical load
Fz0 = m*9.81; % car weigth [N], 1812*9.81 = 17 kN
Fz_slope = 0;

% friction
mu0 = 1;        
mu_slope = 0;   

% velocity
V = 100/3.6; % 100 [Km/h] = 27.7 [m/s]
V0 = 0.1;
w0 = 0.1; 

% acceleration pedal
Tm0 = 2100;

% brakes pedal
Tb0 = 500;

K = 10;

% testname = 'test 2';
% risetime_motor = 5;
% t_stop_acc = 7;
% risetime_brake = 10;
% t_braking = 3;

% save("test.mat", "testname", "risetime_motor", "t_stop_acc", "risetime_brake", "t_braking")

% ------ Plot settings ------ %
F_Size = 14; % FontSize
plotcol = {'k','r','b','g','m','k'};

%% ------------------------- Test 1 ------------------------- %%
% TIP IN - TIP OFF
sim(simulink_model_name);


% subplot(2,1,1)
% plot(sr.Time, Fxr.Data);
% grid on
% subplot(2,1,2)
% plot(sr.Time, Fxr.Data);

%% ------------------------- Test 2 ------------------------- %%
% Acceleration and braking
% testname = 'test 2';
% risetime_motor = 5;
% t_stop_acc = 7;
% risetime_brake = 10;
% t_braking = 7;

%% ------------------------- Test 3 ------------------------- %%
% Emergency braking
% testname = 'test 3';
% risetime_motor = 5;
% t_stop_acc = 7;
% risetime_brake = 10;
% t_braking = 7;










% Time_sim = 1;
% risetime_motor_signal_vec = [1, 5, 10];
% figure
% for i = 1:3 
% % Run simulation
% risetime_motor_signal = risetime_motor_signal_vec(i);
% sim(simulink_model_name);
% plot(sf.Time, Fxr.Data, plotcol{i},'DisplayName', strcat('risetime ', num2str (risetime_motor_signal)));
% hold on
% end
% grid on
% legend('show', 'Location', 'northeast', 'FontSize', 12, 'Box', 'off');



