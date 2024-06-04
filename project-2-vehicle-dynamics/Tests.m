close all; clc; clear; 
warning off

%% ------------------------- Initialization ------------------------- %%
simulink_model_name = 'Tyre_Pac96';
Time_sim = 20; % Simulation time
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

Fzf_ratio = 0.5; % equal weight distribution when V=0
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
torsional_stiffness = 9000; % OPTIONAL!! (NOT USED)
% motor delay
motor_delay = 0.02; % [s]
motor_risetime = 0.05; % [s]

% battery
battery_cap = 58000; % [Wh]
battery_volt = 800; % [V] (NOT USED)
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
% camber angle (NOT USED)
gamma0 = 0; gamma_slope = 0;
% slip angle (NOT USED)
alfa0 = 0; alfa_slope = 0; tau_alfa = 0.1;

% slip
s0 = 0.1;          
Lrel = 0.3; % 0.12<->0.45 [m]
% (tau = Lrel / V)

% friction
mu = 0.85; % asphalt dry        

% velocity
% V = 100/3.6; % 100 [Km/h] = 27.7 [m/s]
V0 = 0.1;
w0 = 0.1; 

% acceleration pedal
Tm_max = max_torque*gear_ratio*motor_eff;
t_stop_acc = 12;

% brakes pedal
Tb_max = Tm_max*4;
t_braking = 12;

% PID ABS
Kp_ABS = 1e5;
Ki_ABS = 100;
Kd_ABS = 10;
slip_reference_ABS = -0.15;
max_Torque_ABS = Tb_max*brake_front*0.5*0.5; 
% ABS can act only on 50% of the rear brakes torque

% PID Torque Control System
Kp_TCS = 200;
Ki_TCS= 100;
Kd_TCS = 100;
slip_reference_TCS = +0.15;
max_Torque_TCS = Tm_max*0.5 *0.5;
% TCS can act only on 50% of the motor torque

test_number = 2;
k_brake_antireverse = 10;

return;




%% ------------------------- Test 1 ------------------------- %%
disp('----------- Test 1 - Longitudinal Acceleration -----------')
test_number = 1;
t_stop_acc = Time_sim;
Time_sim = 5;
t_stop_acc = Time_sim;
mu = 0.9;
V0=0.1;
w0=0.1;
av = [];
plotlegend = cell(2,1);

% PID Torque Control System OFF
Kp_TCS = 0;
Ki_TCS= 0;
Kd_TCS = 0;

sim(simulink_model_name);
av(:, 1) = ax;

% PID Torque Control System ON
Kp_TCS = 200;
Ki_TCS= 100;
Kd_TCS = 100;

sim(simulink_model_name);
av(:, 2) = ax;

close all
PLOT(t, av, 2, {'TCS OFF', 'TCS ON'}, 'a [m/s^2]','TEST 1 - Longitudinal acceleration')

mu=0.85;V0=0.1;w0=0.1;t_stop_acc=12;



%% ------------------------- Test 2 ------------------------- %%
disp('----------- Test 2 - Acceleration Times -----------')
test_number = 2;
Time_sim = 50;
t_stop_acc = Time_sim;
sim(simulink_model_name);


close all
PLOT(t, V*3.6, 1, {}, 'V [Km/h]','TEST 2 - Vehicle velocity')
hold on
plot([t(1) t(end) ], [100 100]) , hold on
plot([t(1) t(end) ], [150 150]) , hold on
plot([t(1) t(end) ], [200 200])
% t_100 = find(V*3.6<=100.1 && V*3.6 >= 99.9)
t_150 = find(V*3.6==150.0)
t_200 = find(V*3.6==200)

%% ------------------------- Test 3 ------------------------- %%
disp('----------- Test 3 - Power Loss -----------')
test_number = 3;
Time_sim = 50;
t_stop_acc = Time_sim;

sim(simulink_model_name);

close all
PLOT(t, V, 1, {}, 'V [Km/h]','TEST 3 - Vehicle velocity')
%rolling resistance, aerodynamic drag, electric powertrain,transmission and longitudinal tyre slip power losses



%% ------------------------- Test 4 ------------------------- %%
disp('---- Test 4 - Energy consumption and achievable range ----')
test_number = 4;
Time_sim = 150;
t_stop_acc = Time_sim;

sim(simulink_model_name);


fprintf('Peak velocity = %3d\n', round(max(V)*3.6))
fprintf('Stopping distance = %.2f [m]\n', delta_x(end))


close all
PLOT(t, V, 1, {}, 'V [Km/h]','TEST 4 - Vehicle velocity')
PLOT(t, E, 1, {}, 'E [kWh]','TEST 4 - Energy Consumption'), hold on
plot([t(1) t(end)], [battery_cap battery_cap]), hold on
plot([t(1) t(end)], [0 0])

%% ------------------------- Test 5 ------------------------- %%
disp('----------- Test 5 - TIP IN - TIP OFF -----------')

test_number = 5;
Time_sim = 20;
V0 = 50/3.6;
w0 = V0/wheel_radius;
t_stop_acc = 1.1;
t_braking = 5;

sim(simulink_model_name);

close all
PLOT(t, ax, 1, {}, 'a [m/s^2]', 'TEST 5 - Longitudinal acceleration')
PLOT(t, V*3.6, 1, {}, 'V [Km/h]','TEST 5 - Vehicle velocity')
PLOT(t, Trr, 1, {}, 'Tr','TEST 5 - Rolling resistance')
PLOT(t, air_drag, 1, {}, 'Fa','TEST 5 - Air drag')
PLOT(t, [sf, sr], 2, {'s front', 's rear'}, 's','TEST 5 - Slip')

V0 = 0.1;w0 = 0.1;

%% ------------------------- Test 6 ------------------------- %%
disp('------ Test 6 - Recuperated Energy  ------')

test_number = 6;
Time_sim = 100;
t_stop_acc = 15;
V0 = 0.1; w0 = 0.1;
sim(simulink_model_name);

close all
PLOT(t, ax, 1, {'ax'}, 'a','TEST 6 - Longitudinal acceleration')
PLOT(t, V*3.6, 1, {'V [Km/h]'}, 'V [Km/h]','TEST 6 - Vehicle velocity')
PLOT(t, E, 1, {'E'},'E', 'TEST 6 - Recuperated Energy')
fprintf('Stopping distance = %.2f [m]\n', delta_x(end)-delta_x(find(t==t_braking)))
fprintf('Recuperated Energy = %.2f\n', E(end)-min(E))
fprintf('Used Energy = %.2f\n', battery_cap-min(E))
fprintf('Total Wasted Energy = %.2f\n', battery_cap-E(end))


%% ------------------------- Test 7 ------------------------- %%
disp('--------- Test 7 - Emergency braking ---------')
disp('mu = 0.90 - Cobblestone')
disp('mu = 0.85 - Asphalt dry')
disp('mu = 0.75 - Asphalt wet')

test_number = 7;
Time_sim = 100;
V0 = 100/3.6;
w0 = V0/wheel_radius-20;
t_stop_acc = 0.2;
t_braking = 0.2;
mu_v = [0.9, 0.85, 0.75];
L = length(mu_v);
s = [];
w = [];
Tbv = [];
Vv = [];
av = [];
tv = [];
plotlegend = cell(L,1);
for i = 1:L
    mu = mu_v(i);
    sim(simulink_model_name);
    if length(tv)<length(t), tv=t; end
    s = padding(s, sf);
    w = padding(w, wf);
    Tbv = padding(Tbv, Tb_PID);
    Vv = padding(Vv, V);
    av = padding(av, ax);
    plotlegend{i} = sprintf('mu = %.2f', mu);
    fprintf('mu = %.2f, Stopping distance = %.2f [m]\n', mu, delta_x(end)-delta_x(find(t==t_braking)))
end

close all
PLOT(t, s, L, plotlegend,'s', 'TEST 7 - Slip Front')
PLOT(t, w, L, plotlegend, 'w [rad/s]','TEST 7 - Front Wheel Speed')
PLOT(t, Tbv, L, plotlegend, 'Tb [Nm]','TEST 7 - Front Brakes')
PLOT(t, Vv*3.6, L, plotlegend, 'V [Km/h]','TEST 7 - Vehicle Speed')
PLOT(t, av, L, plotlegend, 'a [m/s^2]','TEST 7 - Longitudinal acceleration')

mu=0.85;V0=0.1;w0=0.1;t_stop_acc=12;t_braking=12;


%% -------------------- Test 3 (COMPARE ABS ON/OFF) -------------------- %%
disp('--------- Test 3 - COMPARE ABS ON/OFF ---------')
test_number = 3;
Time_sim = 20;

% ABS OFF
Kp_ABS = 0;Ki_ABS = 0;Kd_ABS = 0;
sim(simulink_model_name);
V_OFF = V;
sf_OFF = sf; % sr_OFF = sr;
wf_OFF = wf; % wr_OFF = wr;
t_OFF = t;
fprintf('Stopping distance ABS OFF = %.2f\n', delta_x(end)-delta_x(find(t==t_braking)))

% ABS ON
Kp_ABS = 1e5;Ki_ABS = 100;Kd_ABS = 10;
sim(simulink_model_name);

V = padding(V_OFF, V);
s = padding(sf_OFF, sf);
w = padding(wf_OFF, wf);

fprintf('Stopping distance ABS ON  = %.2f\n', delta_x(end)-delta_x(find(t==t_braking)))

close all
PLOT(t, [Tbf, Tb_PID], 2, {'Tbf ABS OFF', 'Tbf ABS ON'},'Tb [Nm]', 'TEST 3 - Front brakes')
if length(t_OFF)>length(t),t = t_OFF;end
PLOT(t, V*3.6, 2, {'V ABS OFF [Km/h]', 'V [Km/h]','V ABS ON [Km/h]'}, 'TEST 3 - Vehicle velocity')
PLOT(t, s, 2, {'s ABS OFF', 's ABS ON'},'s', 'TEST 3 - Front Slip')
PLOT(t, w, 2, {'w ABS OFF', 'w ABS ON'},'w [rad/s]', 'TEST 3 - Front Wheels angular vel')


%% ------------------------- Tuning ABS ------------------------- %%
disp('-------------------- Tuning ABS --------------------')

test_number = 3;
Time_sim = 20;
k_brake_antireverse = 10;
Kp_v = [0, 1e5, 1e5, 1e5];
Ki_v = [0, 1e2, 1e2, 1e3, 1e4];
Kd_v = [0, 1, 1e2, 1e3, 1];
L = length(Kp_v);
s = [];
w = [];
Tbv = [];
Vv = [];
tv = [];
plotlegend = cell(L,1);
for i = 1:L
    Kp_ABS = Kp_v(i);
    Ki_ABS = Ki_v(i);
    Kd_ABS = Kd_v(i);
    sim(simulink_model_name);
    if length(tv)<length(t), tv=t; end
    s = padding(s, sf);
    w = padding(w, wf);
    Tbv = padding(Tbv, Tb_PID);
    Vv = padding(Vv, V);
    plotlegend{i} = sprintf('Kp=%d-Ki=%d,Kd=%d', Kp_ABS, Ki_ABS, Kd_ABS);
end

close all
PLOT(t, s, L, plotlegend, 's','ABS s')
PLOT(t, w, L, plotlegend, 'w [rad/s]','ABS w')
PLOT(t, Tbv, L, plotlegend, 'Tb [Nm]','ABS Tb')
PLOT(t, Vv*3.6, L, plotlegend,'V [Km/h]', 'ABS V')


%% ------------------------- Tuning Torque control  ------------------------- %%
disp('-------------------- Tuning TCS --------------------')
test_number = 2;
Time_sim = 2;

Kp_v = [0, 1000, 500, 20];
Ki_v = [0, 100, 500, 1000];
Kd_v = [0, 100, 100, 100];
L = length(Kp_v);
s = [];
w = [];
Tmv = [];
Vv = [];
plotlegend = cell(L,1);
for i = 1:L
    Kp_TCS = Kp_v(i);
    Ki_TCS = Ki_v(i);
    Kd_TCS = Kd_v(i);
    sim(simulink_model_name);
    s = padding(s, sr);
    w = padding(w, wr);
    Tmv = padding(Tmv, Tm);
    Vv = padding(Vv, V);
    plotlegend{i} = sprintf('Kp=%d-Ki=%d,Kd=%d', Kp_TCS, Ki_TCS, Kd_TCS);
end

close all
PLOT(t, s, L, plotlegend, 's', 'Torque Control System s')
PLOT(t, w, L, plotlegend, 'w [rad/s]', 'Torque Control System w')
PLOT(t, Tmv, L, plotlegend, 'Tm [Nm]', 'Torque Control System Tm')
PLOT(t, Vv, L, plotlegend, 'V [Km/h]','Torque Control System V')

legend(plotlegend)

%%
% close all, clc
% PLOT(t, V, 1, {'s'}, 's', 'Torque Control System s')

%% -----------------------------------------------------------------------%%

function PLOT(t, x, L, plotlegend, yaxis, nome_fig)
    linewidth = 1;
    plotcol = {'r','b','g','m','k'};
    figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
    if L == 1
        plot(t, x, plotcol{1},'LineWidth', linewidth);
    else
        for i = 1:L
            plot(t, x(:, i), 'Color', plotcol{i},'LineWidth', linewidth);
            hold on
        end
    end
    grid on
    xlabel('time [seconds]');
    ylabel(yaxis);
    title(nome_fig)

    if isempty(plotlegend) == 0
        legend(plotlegend)
    end
    saveas(gcf, strcat('imgs\', nome_fig, '.png'));  % Save as PNG file
end
    
function V = padding(v1, v2)
    
    [L1, N] = size(v1);
    L2 = length(v2);
    
    if L1 == 0
        V1 = v2;
        V2 = [];
        V = [V1, V2];
        return
    end
    
    if L1 > L2
        V1=v1;
        V2=[v2;zeros(L1-L2,1)];
    else
        V1 = zeros(L2, N);
        V1(1:L1,:)=v1;
        V2=v2;
    end
    V = [V1, V2];
    
end