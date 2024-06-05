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
battery_cap = 58*3.6*1e3; % [KWh] -> [K Joule]
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
Lrel = 0.12; % range 0.12<->0.45 [m]
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
Kp_TCS = 500;
Ki_TCS= 100;
Kd_TCS = 100;
slip_reference_TCS = +0.15;
max_Torque_TCS = Tm_max*0.5 *0.5;
% TCS can act only on 50% of the motor torque
Kp_TCS = 0;Ki_TCS= 0;Kd_TCS = 0; % OFF

test_number = 2;
k_brake_antireverse = 10;

return;




%% ------------------------- Test 1 ------------------------- %%
disp('----------- Test 1 - Longitudinal Acceleration -----------')
test_number = 1;

Time_sim = 3;

mu = 0.9; % high friction asphalt
V0=0.1;
w0=0.1;
av = [];
srear = [];
Vv = [];
plotlegend = cell(2,1);

% PID Torque Control System OFF
Kp_TCS = 0;
Ki_TCS= 0;
Kd_TCS = 0;

sim(simulink_model_name);
av(:, 1) = ax;
srear(:,1) = sf;
Vv(:,1) = V;

% PID Torque Control System ON
Kp_TCS = 1000;
Ki_TCS= 100;
Kd_TCS = 100;

sim(simulink_model_name);
av(:, 2) = ax;
srear(:,2) = sf;
Vv(:,2) = V;
close all
PLOT(t, av, 2, {'TCS OFF','TCS ON'}, 'a [m/s^2]','TEST 1 - Longitudinal acceleration')
PLOT(t, srear, 2, {'TCS OFF','TCS ON'}, 's','TEST 1 - Rear Slip')
PLOT(t, Vv*3.6, 2, {'TCS OFF','TCS ON'}, 'V [Km/h]','TEST 1 - Vehicle velocity')
PLOT(t, [Tm, Tm_TCS], 2, {'TCS OFF','TCS ON'}, 'Tm [Nm]','TEST 1 - Motor Torque')

mu=0.85;V0=0.1;w0=0.1;
Kp_TCS = 0;Ki_TCS= 0;Kd_TCS = 0;% TCS OFF


%% ------------------------- Test 2 ------------------------- %%
disp('----------- Test 2 - Acceleration Times -----------')
test_number = 2;
Time_sim = 250;
t_stop_acc = Time_sim;

sim(simulink_model_name);

close all
PLOT(t, V*3.6, 1, {}, 'V [Km/h]','TEST 2 - Vehicle velocity')
find_time_given_vel = @(v,vv,tt) round(tt(find(vv*3.6>=v,1)),1) 
t_100 = find_time_given_vel(100,V,t);
t_150 = find_time_given_vel(150,V,t);
t_200 = find_time_given_vel(200,V,t);
peak_vel = round(max(V)*3.6);
t_peak = find_time_given_vel(max(V)*3.6,V,t);
hold on, plot(t_100, 100, '*b','Markersize', 10) 
hold on, plot(t_150, 150, '*b','Markersize', 10) 
hold on, plot(t_200, 200, '*b','Markersize', 10) 
hold on, plot(t_peak, peak_vel, '*b','Markersize', 10) 

saveas(gcf, strcat('imgs\', 'TEST 2 - Vehicle velocity', '.png'));  % Save as PNG file
fprintf('100 [Km/h], acceleration time = %.2f [sec]\n', t_100)
fprintf('150 [Km/h], acceleration time = %.2f [sec]\n', t_150)
fprintf('200 [Km/h], acceleration time = %.2f [sec]\n', t_200)
fprintf('Peak velocity = %3d [Km/h], at time %.2f [sec]\n', peak_vel , t_peak)

%% ------------------------- Test 3 ------------------------- %%
disp('----------- Test 3 - Power Loss -----------')
test_number = 3;
Time_sim = 20;
t_stop_acc = Time_sim/2;

sim(simulink_model_name);

close all
PLOT(t, V, 1, {}, 'V [Km/h]','TEST 3 - Vehicle velocity')
PLOT(t, [Pot_Trr, Pot_Trf]*1e-3, 2, {'P Tr rear', 'P Tr front'}, 'P [kW]','TEST 3 - Power lost Rolling resistance')
PLOT(t, Pot_air*1e-3, 1, {}, 'P [kW]','TEST 3 - Power lost air drag')
PLOT(t, Pot_mec*1e-3, 1, {}, 'P [kW]','TEST 3 - Power lost electric powertrain')
PLOT(t, Pot_transm*1e-3, 1, {}, 'P [kW]','TEST 3 - Power lost transmission')
PLOT(t, [Pot_Fxr, Pot_Fxf]*1e-3, 2, {'P Fx rear', 'P Fx front'}, 'P [kW]','TEST 3 - Power lost longitudinal slip')
plotlegend = {'Transmission','motor','air drag','Tr rear', 'Tr front'};
PLOT(t, [Pot_transm, Pot_mec, Pot_air, Pot_Trr, Pot_Trf]*1e-3, 5, plotlegend , 'P [kW]','TEST 3 - Power lost')

%% ------------------------- Test 4 ------------------------- %%
disp('---- Test 4 - Energy consumption and achievable range ----')
test_number = 4;
Time_sim = 2000;

t_stop_acc_v = [0.5, 0.7, 1];
L = length(t_stop_acc_v);
Ev = [];
Vv = [];
tv = [];
plotlegend = cell(L,1);
for i = 1:L
    t_stop_acc = t_stop_acc_v(i);
    sim(simulink_model_name);
    if length(tv)<length(t), tv=t; end
    Ev = padding(Ev, E);
    Vv = padding(Vv, V);
    plotlegend{i} = sprintf('maxV = %.2f', round(max(V)*3.6));
    fprintf('Peak velocity = %3d [Km/h]\n', round(max(V)*3.6))
    fprintf('Stopping distance = %.2f [m]\n', delta_x(end))
end

close all
PLOT(tv, Vv, L, plotlegend, 'V [Km/h]','TEST 4 - Vehicle velocity')
PLOT(tv, Ev, L, plotlegend, 'E [KJ]','TEST 4 - Energy Consumption'), hold on
plot([tv(1) tv(end)], [battery_cap battery_cap], '--b'), hold on
plot([tv(1) tv(end)], [0 0], '--b')
plotlegend{L+1} = 'full battery';
plotlegend{L+2} = 'empty battery';
legend(plotlegend)


%% ------------------------- Test 5 ------------------------- %%
disp('----------- Test 5 - TIP IN - TIP OFF -----------')

test_number = 5;
Time_sim = 7;
V0 = 50/3.6;
w0 = V0/wheel_radius;
t_stop_acc = 2;
t_braking = 5;

sim(simulink_model_name);

close all
PLOT(t, ax, 1, {}, 'a [m/s^2]', 'TEST 5 - Longitudinal acceleration')
PLOT(t, V*3.6, 1, {}, 'V [Km/h]','TEST 5 - Vehicle velocity')
PLOT(t, [sf, sr], 2, {'s front', 's rear'}, 's','TEST 5 - Slip')

V0 = 0.1; w0 = 0.1;

%% ------------------------- Test 6 ------------------------- %%
disp('----------- Test 6 - Recuperated Energy  -----------')

test_number = 6;
Time_sim = 100;
t_stop_acc = 15;
V0 = 0.1; w0 = 0.1;
sim(simulink_model_name);

close all
PLOT(t, ax, 1, {'ax'}, 'a','TEST 6 - Longitudinal acceleration')
PLOT(t, V*3.6, 1, {'V'}, 'V [Km/h]','TEST 6 - Vehicle velocity')
PLOT(t, E, 1, {'E'},'E [KJ]', 'TEST 6 - Recuperated Energy'), hold on
plot([t(1) t(end)], [battery_cap battery_cap], '--b')
legend({'Energy','full battery'})
fprintf('Stopping distance = %.2f [m]\n', delta_x(end)-delta_x(find(t==t_braking)))
fprintf('Recuperated Energy = %.2f [KJ]\n', E(end)-min(E))
fprintf('Used Energy = %.2f [KJ]\n', battery_cap-min(E))
fprintf('Total Wasted Energy = %.2f [KJ]\n', battery_cap-E(end))


%% ------------------------- Test 7 ------------------------- %%
disp('-------------- Test 7 - Emergency braking --------------')
disp('mu = 0.85 - Asphalt dry')
disp('mu = 0.75 - Asphalt wet')

test_number = 7;
Time_sim = 100;
V0 = 100/3.6;
w0 = V0/wheel_radius-20;
t_stop_acc = 0.2;
t_braking = 0.2;
mu_v = [0.85, 0.75];
L = length(mu_v);
s = [];
w = [];
Tbv = [];
Vv = [];
av = [];
tv = [];
plotlegend = cell(L,1);
% ABS ON
Kp_ABS = 1e5;Ki_ABS = 100;Kd_ABS = 10;
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
index_short = 1:find(tv>=0.6,1)
t_short = tv(index_short)
PLOT(t_short, s(index_short, :), L, plotlegend,'s', 'TEST 7.1 - Slip Front')
PLOT(t_short, w(index_short, :), L, plotlegend, 'w [rad/s]','TEST 7.1 - Front Wheel Speed')
PLOT(t_short, Tbv(index_short, :), L, plotlegend, 'Tb [Nm]','TEST 7.1 - Front Brakes')
PLOT(tv, Vv*3.6, L, plotlegend, 'V [Km/h]','TEST 7.1 - Vehicle Speed')
PLOT(t_short, av(index_short, :), L, plotlegend, 'a [m/s^2]','TEST 7.1 - Longitudinal acceleration')

mu=0.85;V0=0.1;w0=0.1;


%% -------------------- Test 7 (COMPARE ABS ON/OFF) -------------------- %%
disp('--------- Test 7 - COMPARE ABS ON/OFF ---------')
test_number = 7;
Time_sim = 5;
V0 = 100/3.6;
w0 = V0/wheel_radius;
t_stop_acc = 0.2;
t_braking = 0.2;
mu_v = [0.85, 0.75];
L = length(mu_v);
s = [];
w = [];
Tbv = [];
Vv = [];
av = [];
tv = [];

% ABS OFF
Kp_ABS = 0;Ki_ABS = 0;Kd_ABS = 0;
for i = 1:L
    mu = mu_v(i);
    sim(simulink_model_name);
    if length(tv)<length(t), tv=t; end
    s = padding(s, sf);
    w = padding(w, wf);
    Tbv = padding(Tbv, Tb_PID);
    Vv = padding(Vv, V);
    av = padding(av, ax);
    plotlegend{i} = sprintf('ABS OFF mu = %.2f', mu);
    fprintf('ABS OFF, mu = %.2f, Stopping distance = %.2f [m]\n', mu, delta_x(end)-delta_x(find(t==t_braking)))
end

% ABS ON
Kp_ABS = 1e5;Ki_ABS = 100;Kd_ABS = 10;
for i = 1:L
    mu = mu_v(i);
    sim(simulink_model_name);
    if length(tv)<length(t), tv=t; end
    s = padding(s, sf);
    w = padding(w, wf);
    Tbv = padding(Tbv, Tb_PID);
    Vv = padding(Vv, V);
    av = padding(av, ax);
    plotlegend{i+2} = sprintf('ABS ON  mu = %.2f', mu);
    fprintf('ABS ON,  mu = %.2f, Stopping distance = %.2f [m]\n', mu, delta_x(end)-delta_x(find(t==t_braking)))
end

close all, clc
L=L*2;
index_short = 1:find(tv>=0.6,1)
t_short = tv(index_short)
PLOT(t_short, s(index_short, :), L, plotlegend,'s', 'TEST 7.2 - Slip Front')
PLOT(t_short, w(index_short, :), L, plotlegend, 'w [rad/s]','TEST 7.2 - Front Wheel Speed')
PLOT(t_short, Tbv(index_short, :), L, plotlegend, 'Tb [Nm]','TEST 7.2 - Front Brakes')
PLOT(tv, Vv*3.6, L, plotlegend, 'V [Km/h]','TEST 7.2 - Vehicle Speed')
PLOT(t_short, av(index_short, :), L, plotlegend, 'a [m/s^2]','TEST 7.2 - Longitudinal acceleration')

mu=0.85;V0=0.1;w0=0.1;t_stop_acc=12;t_braking=12;

%% ------------------------- Tuning ABS ------------------------- %%
disp('------------------- Tuning ABS ------------------')

test_number = 7;
Time_sim = 2;
V0 = 100/3.6;
w0 = V0/wheel_radius;
t_stop_acc = 0.2;
t_braking = 0.2;
k_brake_antireverse = 10;
Kp_v = [0, 1e5, 1e10, 1e5];
Ki_v = [0, 1, 1, 100, 1];
Kd_v = [0, 1, 1, 1, 1];
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
PLOT(tv, s, L, plotlegend, 's','ABS s')
PLOT(tv, w, L, plotlegend, 'w [rad/s]','ABS w')
PLOT(tv, Tbv, L, plotlegend, 'Tb [Nm]','ABS Tb')
PLOT(tv, Vv*3.6, L, plotlegend,'V [Km/h]', 'ABS V')


%% ------------------------- Tuning Torque control  ------------------------- %%
disp('------------------- Tuning TCS ------------------')
test_number = 1;
Time_sim = 2;

Kp_v = [0, 1000, 500, 2000];
Ki_v = [0, 100, 500, 1000];
Kd_v = [0, 100, 100, 100];
L = length(Kp_v);
s = [];
w = [];
Tmv = [];
Vv = [];
tv = [];
plotlegend = cell(L,1);
for i = 1:L
    Kp_TCS = Kp_v(i);
    Ki_TCS = Ki_v(i);
    Kd_TCS = Kd_v(i);
    sim(simulink_model_name);
    if length(tv)<length(t), tv=t; end
    s = padding(s, sr);
    w = padding(w, wr);
    Tmv = padding(Tmv, Tm_TCS);
    Vv = padding(Vv, V);
    plotlegend{i} = sprintf('Kp=%d-Ki=%d,Kd=%d', Kp_TCS, Ki_TCS, Kd_TCS);
end

close all
PLOT(tv, s, L, plotlegend, 's', 'Torque Control System s')
PLOT(tv, w, L, plotlegend, 'w [rad/s]', 'Torque Control System w')
PLOT(tv, Tmv, L, plotlegend, 'Tm [Nm]', 'Torque Control System Tm')
PLOT(tv, Vv, L, plotlegend, 'V [Km/h]','Torque Control System V')



%%
close all, clc
PLOT(t, [V,V,V,V,V,V], 5, {}, 'V', 'Prova')

%% -----------------------------------------------------------------------%%

function PLOT(t, x, L, plotlegend, yaxis, nome_fig)
    linewidth = 1.4;
    plotcol = {'r','b','g','m','k'};
    figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
    if L == 1
        plot(t, x, plotcol{1},'LineWidth', linewidth);
    else
        for i = 1:L
            ppp = plot(t, x(:, i), 'Color', plotcol{i},'LineWidth', linewidth);
            ppp.Color(4) = 0.5;

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