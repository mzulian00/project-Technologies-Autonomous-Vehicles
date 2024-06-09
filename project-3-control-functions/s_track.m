%% Technologies for Autonomous Vehicles
% SINGLE TRACK MODEL ========== FULL Version ============
% close all; 
clear; clc
UIControl_FontSize_bak = get(0, 'DefaultUIControlFontSize');
% Select simulation
set(0, 'DefaultUIControlFontSize', 18);
caso = menu('Choose CoG position','a = L/2-10cm','a = L/2','a = L/2+10cm','Drift: a=b; Ca=Ca_nom/100');       
set(0, 'DefaultUIControlFontSize', UIControl_FontSize_bak);
if isempty(caso)
    caso=1;
end
% 4 cases depending on COG position: 
% (1) a = L/2-10cm 
% (2) a = L/2
% (3) a = L/2+10cm 
% (4) a=b; Cr=Cr_nom/100

%% --------- Vehicle and tyre DATA ----------
m = 1997.6*1.2;       %[kg] mass
L = 2.85;       % [m] wheelbase
a_vet =[L/2-0.1 L/2 L/2+0.1 L/2];  %[m] front wheelbase
a = a_vet(caso);  
b = L-a ;       %[m] rear wheelbase        
Tf=1.54;        %[m] front track
Tr=1.54;        %[m] rear track
Jz=3728;        %[kg m^2] mass moment of inertia 
g = 9.81;       %[m/s^2]
tau_s = 15;     % steering ratio

%----- static loads
FzF = m*g*b/L; 
FzR = m*g*a/L;
Perc_F = FzF /(m*g)*100;
Perc_R = 100-Perc_F;
disp(['static load distribution (%F - %R): ',num2str(round(Perc_F)),'-',num2str(round(Perc_R))])

%------ Cornering stiffness: no load transfer
eval(['load Dati_txt' filesep 'CornStiff_Vs_Fz'])
% interpolate wheel cornering stiffness versus vertical load
% CF_w = interp1(Fz_vet,C_alpha_vet,FzF/2);
% CR_w = interp1(Fz_vet,C_alpha_vet,FzR/2);
DFz_tot = 0*4000;
p_F = 0.2; 
DFz_F = p_F*DFz_tot;
DFz_R = (1-p_F)*DFz_tot;

FzFL = FzF/2 - DFz_F;
FzFR = FzF/2 + DFz_F;
FzRL = FzR/2 - DFz_R;
FzRR = FzR/2 + DFz_R;

CF_wL = interp1(Fz_vet,C_alpha_vet,FzFL);
CF_wR = interp1(Fz_vet,C_alpha_vet,FzFR);
CF = CF_wL + CF_wR;

CR_wL = interp1(Fz_vet,C_alpha_vet,FzRL);
CR_wR = interp1(Fz_vet,C_alpha_vet,FzRR);
CR = CR_wL + CR_wR;

% 
% 
% % axle stiffness
% CF = 2*CF_w;    % front
if caso==4
    CR_wL = CR_wL/100 % /100 for drift simulation
    CR = 2*CR_wL;
    C_alpha_vet_drift = C_alpha_vet/100;
else
    CR = 2*CR_wL;    % rear
end

figure; hold on
plot(FzFL,CF_wL/1000,'rs','markersize',15,'linewidth',2,'DisplayName','$C_{FL}$'); 
plot(FzRL,CR_wL/1000,'ko','markersize',15,'linewidth',2,'DisplayName','$C_{RL}$'); 
plot(Fz_vet,C_alpha_vet/1000,'linewidth',2,'DisplayName','$C_{\alpha,FL}$')

if caso ==4
plot(Fz_vet,C_alpha_vet_drift/1000,'linewidth',2,'DisplayName','$C_{\alpha,RL}$')
end

grid on
set(gca,'FontName','Times New Roman','FontSize',14)
legend('interpreter','latex','FontSize',20,'location','best')
xlabel('F_Z [N]'); ylabel('C_{\alpha} [kN/rad]')
set(gcf,'Position',[150 38 985 667]);

disp('axle cornering stiffness')
disp(['CF = ',num2str(CF),' N/rad'])
disp(['CR = ',num2str(CR),' N/rad'])
disp(' ')
% check under/over steering
if CR*b-CF*a>0
    disp('============ understeering ============')
elseif CF*a-CR*b==0
    disp('============ neutral vehicle ============')
else
    disp('============ oversteering ============')
    V_cr = sqrt(CF*CR*L^2/(m*(a*CF-b*CR)))*3.6;
    disp(['critical speed: ',num2str(round(V_cr*10)/10),' km/h'])
end
% understeering and slip angle gradients
mf = m*b/L; mr = m-mf;
K_us_an = (mf/CF-mr/CR);
K_beta_an = -mr/CR;
% K_beta_an_Guiggiani = -m/L^2*(CF*a^2+CR*b^2)/(CF*CR);

% tangent speed (beta=0)
V_beta0 = sqrt(b*L*CR/a/m)*3.6;
disp(' ')
disp(['understeering gradient: K_US = ',num2str(K_us_an),' rad/(m/s^2)'])
disp(['slip angle gradient: K_beta = ',num2str(K_beta_an),' rad/(m/s^2)'])
disp(['tangent speed: V_beta = ',num2str(V_beta0),' km/h'])


pause

%% SINGLE TRACK model analysis:
% Poles, Damping Factors, Natural Frequencies versus velocity
% vehicle speed
vel=[[1:0.1:20],[20:0.5:50],[50:5:450]]/3.6;
if caso==4
    vel=[1:0.1:30]/3.6;
end
delta_f = 1*pi/180;
delta_f = 1;
delta_r = 0;
% matrix initialization 
POLES=[zeros(length(vel),2)];       % Poles
ZETA=[zeros(length(vel),2)];        % damping fators
FREQ=[zeros(length(vel),2)];        % natural frequencies
DET_A=[zeros(length(vel),1)];       % determinant of A
TR_A=[zeros(length(vel),1)];        % trace of A
X_r = [zeros(length(vel),2)];
Y_r = [zeros(length(vel),6)];
ay_r = [zeros(length(vel),1)];
u_r = [delta_f;0];

% Definition of state space variables
StateNames ={'\beta','r'};
InputNames={'\delta_F','\delta_R'};
OutputNames={'\beta','r','\rho','\alpha_F','\alpha_R','a_y'};

% for cycle (vehicle speed)
for k=1:length(vel)
    Vv=vel(k);      % vehicle speed
    % state space matrices: A,B,C,D
    A=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
        (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
    B=[CF/(m*Vv) CR/(m*Vv);
        (CF*a/Jz) -(CR*b/Jz)];
    C= [1,0
        0,1
        (-CR-CF)/(m*Vv^2),(-CF*a+CR*b)/(m*Vv^3)
        -1, -a/Vv
        -1, b/Vv
        (-CR-CF)/(m),(-CF*a+CR*b)/(m*Vv)];
    D = [0 0;
        0 0;
        CF/(m*Vv^2) CR/(m*Vv^2)
        1 0
        0 1
        CF/m CR/m];
    % state space system
    G=ss(A,B,C,D);
    % Wn, Z e P: fnatural frequencies, damping factors and poles
    [Wn,Z,P]=damp(G);
    % storing values in matrices (k-th column)
    POLES(k,:)=P;
    ZETA(k,:)=Z;
    FREQ(k,:)=Wn;
    % Determinant and trace of matrix A
    DET_A(k)=det(A);
    TR_A(k)=trace(A);
    % steady state response to delta_F
    X_r(k,:) = -A^-1*B*u_r;
    Y_r(k,:) = C*X_r(k,:)'+D*u_r;
    ay_r(k,:) = X_r(k,2)'*Vv;
%     [V,D] = eig(A);
end

% poles
P1=POLES(:,1);
P2=POLES(:,2);
% real and imaginary part (first pole)
P1_Real=real(P1);
P1_Im=imag(P1);
% real and imaginary part (second pole)
P2_Real=real(P2);
P2_Im=imag(P2);

fig_1=figure('Name','Poles, Damping, Natural Frequencies','NumberTitle','off','PaperType','A4');
figure(fig_1)
subplot(2,2,[3 4])
% Plot eigenvalues vs speed
scatter(P1_Real,P1_Im,[],[1:1:length(vel)])
hold on
grid on
scatter(P2_Real,P2_Im,[],[1:1:length(vel)])
xlabel('Real'),ylabel('Im'),title('Poles'),
% xlim([-200 20]),ylim([-10 10])
plot([0,0],ylim,'--k')
hold on
colorbar %axis equal
set(gca,'FontName','Times New Roman','FontSize',14)

subplot(2,2,1),plot(vel*3.6,FREQ,'LineWidth',2)
title('Natural Frequencies'),xlabel('vel [km/h]'),ylabel('Wn [rad/s]')
grid on
set(gca,'FontName','Times New Roman','FontSize',14)

subplot(2,2,2),plot(vel*3.6,ZETA,'square')
xlabel('vel [km/h]'),ylabel('Damping Ratio'),title('Damping Ratio');
grid on
%plot(POLES,'x'),hold on,
set(gca,'FontName','Times New Roman','FontSize',14)

% Plot determinant versus speed
fig_11=figure('Name','Determinant of [A]','NumberTitle','off','PaperType','A4');
figure(fig_11)
hold all
plot(vel*3.6,DET_A,'linewidth',2),
plot(vel*3.6,zeros(length(vel),1),'--k')
xlabel('vel [km/h]'),ylabel('det(A)'),title('$\lambda_1 \lambda_2$','interpreter','latex','Fontsize',18)
ylim([-10 100])
set(gca,'FontName','Times New Roman','FontSize',14)

figure('Name','steady state response vs. velocity')
subplot(2,2,1)
plot(vel*3.6,Y_r(:,1),'o','linewidth',2,'Displayname','$\beta$'); %xlabel('vel [km/h]'),%ylabel('Y_r'),
title('$\beta /\delta_F$','interpreter','latex','Fontsize',18); hold on; 
plot(xlim,[0 0],'--k'); grid on
% legend('\beta [rad]','r [rad/s]','\rho [1/m]','\alpha_f [rad]','\alpha_r [rad]')

subplot(2,2,2)
plot(vel*3.6,Y_r(:,2),'o','linewidth',2)
title('$r /\delta_F$','interpreter','latex','Fontsize',18); grid on

subplot(2,2,3)
plot(vel*3.6,Y_r(:,3),'o','linewidth',2); title('$\rho /\delta_F$','interpreter','latex','Fontsize',18); xlabel('vel [km/h]'); grid on

subplot(2,2,4)
plot(vel*3.6,Y_r(:,4:5),'o','linewidth',1);legend('$\alpha_f /\delta_F$','$\alpha_r /\delta_F$','Location','best','interpreter','latex','Fontsize',18)
xlabel('vel [km/h]')
grid on

figure('name','lateral acceleration')
plot(vel*3.6,Y_r(:,6),'o','linewidth',2); 
title('$a_y /\delta_F$','interpreter','latex','Fontsize',18)
xlabel('vel [km/h]'); grid on; xlabel('vel [km/h]')
delta_meno_delta_0 = u_r(1)-L*Y_r(:,3);
beta_meno_beta_0 = Y_r(:,1)-b*Y_r(:,3);

ay_r = Y_r(:,6);
figure('name','understeering diagram')
subplot(1,2,1)
plot(ay_r,delta_meno_delta_0*180/pi*tau_s,'linewidth',2)
xlabel('a_y [m/s^2]'); grid on; 
ylabel('$\delta-\delta_0$ [deg]','interpreter','latex','Fontsize',18); xlim([0 5])
subplot(1,2,2)
plot(ay_r,beta_meno_beta_0*180/pi,'linewidth',2); 
xlabel('a_y [m/s^2]'); grid on; 
ylabel('$\beta-\beta_0$ [deg]','interpreter','latex','Fontsize',18); xlim([0 5])

K_us_num = (delta_meno_delta_0(end)- delta_meno_delta_0(1))/(ay_r(end)-ay_r(1))
K_beta_num = (beta_meno_beta_0(end)- beta_meno_beta_0(1))/(ay_r(end)-ay_r(1))

pause

%% BODE diagram
% speed selection for Bode diagrams
clear A B C D
if caso == 4
    V=input('Select speed for Bode diagrams (default = 5 km/h): ');
else
    V=input('Select speed for Bode diagrams (default = 40 km/h): ');
end
V=V/3.6;
if isempty(V)
    V=40/3.6;
    if caso == 4
        V=5/3.6;
    end
end

% state space
Vv = V;
% state space matricesi: A,B,C,D
    A=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
        (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
    B=[CF/(m*Vv) CR/(m*Vv);
        (CF*a/Jz) -(CR*b/Jz)];
    C = [1,0
        0,1
        (-CR-CF)/(m*Vv^2),(-CF*a+CR*b)/(m*Vv^3)
        -1, -a/Vv
        -1, b/Vv
        (-CR-CF)/(m),(-CF*a+CR*b)/(m*Vv)];
    D = [0 0;
        0 0;
        CF/(m*Vv^2) CR/(m*Vv^2)
        1 0
        0 1
        CF/m CR/m];
%
% choice of Bode outputs
disp(' -------------- ')
disp('choose indices of the desired outputs ')
disp(['1=',(OutputNames{1}),'; 2=',OutputNames{2},'; 3=',OutputNames{3},'; 4=',OutputNames{4},'; 5=',OutputNames{5},'; 6=',OutputNames{6}])
out_index = input('Default [1:2] ');
if isempty(out_index)==1
    out_index = 1:2;
end
% out_index = 1:2;
sys=ss(A,B,C(out_index,:),D(out_index,:),...
    'StateName',StateNames,'InputName',InputNames,'OutputName',OutputNames(out_index));

% figure; step(sys); figure; ltiview(sys)
% G1=ss2tf(sys(1,1)),
% G2=ss2tf(sys(2,2)),
fig_8=figure('Name','Bode Plots','NumberTitle','off','PaperType','A4');
figure(fig_8)
h = bodeplot(sys);
% Change units to Hz and make phase plot invisible
setoptions(h,'FreqUnits','Hz','PhaseVisible','on','Grid','On','Xlim',[0.1 10],'MagScale','linear','MagUnits','abs'); %,'PhaseWrapping','on'
set(findall(gcf,'-property','FontSize'),'FontSize',13)

pause

%% run Simulink model
% default parameters
t_end_sim = 20;
delta_vol_max_ramp = 0;
dvol_max = 0;
t_drift_cs_1 = 0; t_drift_cs_2 = 0;
dvol_drift1 = 0; dvol_drift2 = 0; dvol_drift3 = 0;
% eval(['load Dati_txt' filesep 'sweep.txt'])

set(0, 'DefaultUIControlFontSize', 18);
sel_man = menu('select manoeuvre',' Step steer','Ramp steer 15 deg/s','Sine Sweep','Drift (with reduced rear cornerng stiffness)');
set(0, 'DefaultUIControlFontSize', UIControl_FontSize_bak);
if isempty(sel_man)
    sel_man=1;
end
switch sel_man
    case 1
        %         (Step steer)
        t_end_sim = 10*1;      % [s]
        Vd = 100;                % [km/h]
        dvol_max = input('max steering angle (deg) [default = 100 deg]:');
    case 2
        %        (ramp steer)
        t_end_sim = 36*0+20;
        Vd = 40;                % [km/h]
        delta_vol_max_ramp = input('max steering angle (deg) [default = 200 deg]:');;
    case 3
        %        (Sine Sweep)
        t_end_sim=5;
        Vd = 100;               % [km/h]
        dvol_max = 20;              % [deg] steering wheel angle
    case 4
        Vd = 40;                % [km/h]
        t_drift_cs_1 = 5*1 +3*0;
        t_drift_cs_2 = 14.38*1+ 7.136*0;
        dvol_drift= [30 -85*1+(-85*0.6*0) 0];
        dvol_drift1 = dvol_drift(1);
        dvol_drift2 = -dvol_drift(2)+dvol_drift(1);
        dvol_drift3 = dvol_drift(3)-dvol_drift(2);
end
if isempty(dvol_max)==1
    dvol_max=10; % max steering angle
end
if isempty(delta_vol_max_ramp)==1
    delta_vol_max_ramp = 200; % max steering angle
end
%
V = input(['choose speed for Simulink manoeuvre (default = ' num2str(Vd) ' km/h): ']);
if isempty(V)
    V=Vd;       % km/h
end
V=V/3.6; 
Vv=V;
% matrix initialization
% state space matricesi: A,B,C,D
    A_sim=[(-CF-CR)/(m*Vv),(-CF*a+CR*b-m*Vv^2)/(m*Vv^2);
        (-CF*a+CR*b)/Jz,(-CF*a^2-CR*b^2)/(Jz*Vv)];
    B_sim=[CF/(m*Vv) CR/(m*Vv);
        (CF*a/Jz) -(CR*b/Jz)];
    C_sim = [1,0
        0,1
        (-CR-CF)/(m*Vv^2),(-CF*a+CR*b)/(m*Vv^3)
        -1, -a/Vv
        -1, b/Vv
        (-CR-CF)/(m),(-CF*a+CR*b)/(m*Vv)];
    D_sim = [0 0;
        0 0;
        CF/(m*Vv^2) CR/(m*Vv^2)
        1 0
        0 1
        CF/m CR/m];
%

% Run simualtion
dt_sim = 1e-3;          % time step
ay_max = 10;            % [m/s^2] limit acceleration: stop simulation
beta_max = 80;          % [deg] limit slip angle: stop simulation
sim('s_track_model')    % run Simulink model

%% POST PROCESSING
%--------- Plot Steering angle
% if exist('fig_6')==0
%     fig_6 = figure('Name','Steering Angle','NumberTitle','off','PaperType','A4');
% else
%     figure(fig_6)
% end
figure('Name','steering angle')
hold all; grid on
plot(delta_steer,'LineWidth',2),xlabel('time [s]'),
hold on; plot(L*ro*180/pi*tau_s,'--k'); 
legend('\delta','\delta_0','Fontsize',18,'location','best')
title('Steering Angle \delta_s')

%--------- Plot beta and psi_dot
figure('Name','States')
hold all; grid on
subplot(2,1,1),plot(beta,'LineWidth',2),xlabel('time [s]')
hold on
plot(b*ro*180/pi,'--k'); 
legend('\beta','\beta_0','Fontsize',16,'location','best')
title('slip angle \beta [deg]','Fontsize',16)
grid on
subplot(2,1,2)
plot(r,'LineWidth',2)
hold on
plot(delta_steer/180*pi,'LineWidth',2),xlabel('time [s]'),
title('r [deg/s]','Fontsize',16), xlabel('time [s]')
ylabel('')
legend('r [deg/s]','\delta [rad]','Fontsize',16,'location','best')
%ylim([-150 150])
grid on

%--------- Plot ay vs t
figure('Name','a_y(t)')
hold all; grid on
plot(ay,'LineWidth',2)
title('Lateral Acceleration a_y [m/s^2]','Fontsize',16)
xlabel('time [s]')
ylabel('')
% ---- plot beta beta_dot
figure('Name','\beta, \beta_dot')
hold all; grid on
plot(beta.Data,beta_dot(:,2),'LineWidth',2) %title('Lateral Acceleration a_y [m/s^2]'),
xlabel('\beta [deg]','Fontsize',18)
ylabel('\beta_{dot} [rad/s]','Fontsize',18)

if sel_man == 2  % plot vs ay only for ramp steer 
%--------- Plot beta vs ay
figure('Name','beta vs ay')
% plot(a_y(:,2),Beta(:,2))
scatter(a_y(:,2),Beta(:,2),[],a_y(:,1)); colorbar
xlabel('a_y [m/s^2]')
ylabel('\beta [deg]'); grid on
text(11,2.2,['time[s]'])
set(gca,'FontName','Times New Roman','FontSize',16)
%--------- Plot delta vs ay
figure('Name','delta vs ay')
% plot(a_y(:,2),Beta(:,2))
scatter(a_y(:,2),delta_rad(:,2)*180/pi,[],a_y(:,1)); colorbar
xlabel('a_y [m/s^2]'); ylabel('\delta_{vol} [deg]'); grid on
text(11,22,['time[s]'])
set(gca,'FontName','Times New Roman','FontSize',16)
%--------- Plot delta-delta0 vs ay
figure('Name','delta-delta_0 vs ay')
delta0 =L*ro.Data;
% plot(a_y(:,2),Beta(:,2))
plot(a_y(:,2),(delta_rad(:,2)-delta0*tau_s)*180/pi,'linewidth',2); 
xlabel('a_y [m/s^2]'); ylabel('\delta_{vol}-\delta_0 [deg]'); grid on
set(gca,'FontName','Times New Roman','FontSize',16)
end

%--------- Plot alpha_F e alpha_R 
figure('Name','alphaF e R'); hold all
plot(alfaF*180/pi,'LineWidth',2); plot(alfaR*180/pi,'LineWidth',2); 
xlabel('time [s]'); ylabel('\alpha [deg]'); grid on; 
legend('\alpha_F','\alpha_R','Fontsize',16,'location','best')
set(gca,'FontName','Times New Roman','FontSize',14)
legend({},'FontSize',16)

%--------- Plot curvature
figure('Name','rho'); hold all
plot(ro,'LineWidth',2); 
xlabel('time [s]'); ylabel('\rho [1/m]'); grid on; 
set(gca,'FontName','Times New Roman','FontSize',14)
%--------- Plot trajectory
spost_x=Var_trajectory(:,1);
spost_y=Var_trajectory(:,2);

figure
hold all; grid on
% plot(spost_x,spost_y,'LineWidth',2)
scatter(spost_x,spost_y,[],a_y(:,1))
title('trajectory'),axis equal,xlabel('X [m]'),ylabel('Y[m]');colorbar
text(49,12,['time[s]'])

pause(3)

%% Vehicle Trajectory and Animation
F_Size=13;

figure('Name','Vehicle CG location','NumberTitle','off','PaperType','A4');
hold all; grid on
plot(Var_trajectory(:,1),Var_trajectory(:,2),'Linewidth',2);
title('CG Trajectory'); axis equal
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('X [m]'); ylabel('Y [m]')
X_G = Var_trajectory(:,1);
Y_G = Var_trajectory(:,2);

% figure
%hold on
axis equal
dt_frame = 0.1; % [s] vehicle frame refresh (time) interval 
decim_frame = dt_frame/dt_sim;

% cycle to show vehicle motion
for cont1=1:decim_frame:length(Psi)
    X = X_G(cont1)+[-b a a -b]';
    Y = Y_G(cont1)+[-Tf/2 -Tf/2 Tf/2 Tf/2]';
    vert = [X,Y];
    fac = [1 2 3 4];
    hVeicolo = patch('Faces',fac,'Vertices',vert,'FaceColor','red','FaceAlpha',.5);
    direction = [0 0 1];
    xlim([X(1)-5 X(2)+3])
    ylim([Y(1)-5 Y(3)+3])
    
    x0 = X_G(cont1);
    y0 = Y_G(cont1);
    z0 = 0;
    ORIGIN = [x0,y0,z0];
    rotate(hVeicolo,direction,Psi(cont1,2),ORIGIN);
    pause(0.1)
end

plot(Var_trajectory(:,1),Var_trajectory(:,2),'Linewidth',2);
axis auto

%% TF estimation 
if sel_man == 3
    fs=1/dt_sim; % Hz
    t_wind = 5; % s
    nfft=t_wind*fs;	         	% number of points (fft 512)
    nolap=0.9*nfft;				% overlapping percentage (Hanning window) 90%
    [T_ay_d,freq1] = tfestimate(delta_steer.Data,ay.Data,nfft,nolap,[],fs);
    figure
    plot(freq1,abs(T_ay_d))
    grid on
    xlim([0 5])
    xlabel('frequency [Hz]')
    ylabel('a_y/\delta [g/deg]')
end
