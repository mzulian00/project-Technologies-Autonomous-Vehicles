close all; clc; clear; 
warning off

%% ------------------------- Initialization ------------------------- %%
simulink_model_name = 'Tyre_Pac96';
Time_sim = 12; % Simulation time

% Load tyre data (Pacejka model coefficients)
WheelFile = 'Tyre215_50_19_Comb';
eval(['[Pacejka]=' WheelFile ';'])
pacn = struct2cell(Pacejka);
for ii = 1:size(pacn) Pace(ii)=pacn{ii}; end
Pacn = Pace';     


% ------ vehicle parameters ------ %
m = 1812; % kerb weight [Kg]
wheelbase = 2.77; %  [m]
cg_height = 0.55; % center of gravity height [m]
% TODO mancano valori dimensione car
% ...
Fzf_ratio= 0.5;
Fzr_ratio = 0.5;

% wheel
wheel_radius = 0.3; % [m] % TODO find correct value 
Ir = 0.675; % wheel inertia TODO find value [0.675 è di chatgpt]

% motor
peak_power = 150; % [kW]
max_torque = 310; % [Nm]
max_speed = 16000; % [rpm]
gear_ratio = 10.5;
motor_eff = 0.9;
transm_eff = 0.95;
torsional_stiffness = 9000; % OPTIONAL!!

% motor delay
motor_delay = 0.02; % [s]
motor_risetime = 0.05; % [s]

% battery
battery_cap = 58; % [KWh]
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
V = 100/3.6; % 100 [Km/h]
V0 = 0.1;
w0 = 0.1; 

% acceleration pedal
Tm0 = 2100;
accel_time = 0.001;

% brakes pedal
Tb0 = 500;
brake_time = 2;






% ------ Plot settings ------ %
% F_Size = 14; % FontSize
% plotcol = {'--k','-r','-.b',':g','.m','-k'};



%% ------------------------- Test 1 ------------------------- %%

% Run simulation
sim(simulink_model_name);







