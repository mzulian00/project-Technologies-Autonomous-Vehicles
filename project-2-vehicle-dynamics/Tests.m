close all; clc; clear; 
warning off

%% ------------------------- Initialization ------------------------- %%
simulink_model_name = 'Tyre_Pac96';
Time_sim = 12; % Simulation time
early_stop = 1; % Stop if Vel vehicle == 0

% Load tyre data (Pacejka model coefficients)
WheelFile = 'Tyre215_50_19_Comb';
eval(['[Pacejka]=' WheelFile ';'])
pacn = struct2cell(Pacejka);
for ii = 1:size(pacn) Pace(ii)=pacn{ii}; end
Pacn = Pace';

% ------ vehicle parameters ------ %
m = 1812; % kerb weight [Kg], car weigth = 1812*9.81 = 17 [kN]
g = 9.81; % gravity [m/s^2]
wheelbase = 2.77; %  [m]
cg_height = 0.55; % center of gravity height [m]

Fzf_ratio= 0.5; % equal weight distribution when V=0
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
inverter_eff = 0.9;

% brakes
brake_delay = 0.02; % [s]
brake_risetime = 0.025; % [s]
brake_front = 0.75; % brake distribution 75:25 
brake_rear = 0.25;

% rolling resistance
f0 = 0.009; 
f2 = 6.5e-6; % [s^2/m^2]

% air drag
rho = 1.225; % [kg/m³] air density 
Af = 2.36; % [m^2] frontal area
Cx = 0.27; % drag coefficient

% ------ signals coefficients setting ------ %
% camber angle
gamma0 = 0; gamma_slope = 0;

% slip
s0 = 0.1;         
s_slope = 0;      
tau_slip = 0.1; % (tau = Lrel / V)

% slip angle
alfa0 = 0; alfa_slope = 0; tau_alfa = 0.1;

% vertical load
Fz0 = m*9.81; % 
Fz_slope = 0;

% friction
mu = 1;        

% velocity
V = 100/3.6; % 100 [Km/h] = 27.7 [m/s]
V0 = 0.1;
w0 = 0.1; 

% acceleration pedal

Tm_max = max_torque*gear_ratio*motor_eff;

% brakes pedal
Tb_max = 4000;

% PID ABS
Kp_ABS = 100;
Ki_ABS = 100;
Kd_ABS = 10;
slip_reference_ABS = -0.15;
max_Torque_ABS = Tb_max*brake_rear*0.5*0.8; 
% ABS can act on 80% of the rear brakes torque

% PID Torque Control System
Kp_TCS = 200;
Ki_TCS= 100;
Kd_TCS = 100;
slip_reference_TCS = +0.15;
max_Torque_TCS = Tm_max*0.5*0.5; 
% TCS can act on 50% of the motor torque

% ------ Plot settings ------ %
F_Size = 14; % FontSize
plotcol = {'r','b','g','m','k'};



%% ------------------------- Test 1 ------------------------- %%
% TIP IN - TIP OFF

test_number = 1;
Time_sim = 20;
sim(simulink_model_name);





%% ------------------------- Test 2 ------------------------- %%
% Acceleration and braking

test_number = 2;

%% ------------------------- Test 3 ------------------------- %%
% Emergency braking

test_number = 3;



%% ------------------------- Tuning Torque control  ------------------------- %%
% test_number = 2;
% Time_sim = 2;
% 
% Kp_v = [0, 200];
% Ki_v = [0, 100, 50, 100];
% Kd_v = [0, 100, 50, 100];
% L = length(Kp_v);
% s = [];
% w = [];
% Tmv = [];
% Vv = [];
% plotlegend = cell(L,1);
% for i = 1:L
%     Kp_TCS = Kp_v(i);
%     Ki_TCS = Ki_v(i);
%     Kd_TCS = Kd_v(i);
%     sim(simulink_model_name);
%     s(i,:) = sr;
%     w(i,:) = wr;
%     Tmv(i, :) = Tm;
%     Vv(i, :) = V;
%     plotlegend{i} = sprintf('Kp=%d-Ki=%d,Kd=%d', Kp_TCS, Ki_TCS, Kd_TCS);
% end
% 
% PLOT(t, s, L, plotlegend, 'Torque Control System s')
% PLOT(t, w, L, plotlegend, 'Torque Control System w')
% PLOT(t, Tmv, L, plotlegend, 'Torque Control System Tm')
% PLOT(t, Vv, L, plotlegend, 'Torque Control System V')
% 
% legend(plotlegend)
%% ------------------------- Tuning ABS ------------------------- %%
test_number = 2;

% Tb_max = 4000;
% Kp_v = [0, 100];
% Ki_v = [0, 100, 100, 100];
% Kd_v = [0, 10, 10, 100];
% L = length(Kp_v);
% s = [];
% w = [];
% plotlegend = cell(L,1);
% for i = 1:L
%     Kp = Kp_v(i);
%     Ki = Ki_v(i);
%     Kd = Kd_v(i);
%     sim(simulink_model_name);
%     s(i,:) = sr;
%     w(i,:) = wr;
%     plotlegend{i} = sprintf('Kp=%d-Ki=%d,Kd=%d', Kp_ABS, Ki_ABS, Kd_ABS);
% end
% 
% PLOT(t, s, L, plotlegend, 'ABS s')
% PLOT(t, w, L, plotlegend, 'ABS w')


%% ------------------------------------------------------------------


function PLOT(t, x, L, plotlegend, nome_fig)
F_Size = 14; % FontSize
plotcol = {'r','b','g','m','k'};
figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
for i = 1:L
    plot(t, x(i, :), plotcol{i});
    hold on
end
grid on
legend(plotlegend)

end

