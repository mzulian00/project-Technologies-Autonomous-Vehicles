%%
clc, clear all, close all, warning off

%% ---------------- Initialization ----------------
t_end_sim = 300; 
dt_sim = 1e-3;
model_name = 'path_tracking';

m = 1575;       %[kg] mass
a = 1.3;  
b = 1.5;        %[m] rear wheelbase       
L = a+b;        %[m] wheelbase
Jz=2875;        %[kg m^2] mass moment of inertia 
g = 9.81;       %[m/s^2]
Cf = 2*60000;   %Front axle cornering stiffness
Cr = 2*57000;   %Rear  axle cornering stiffness
max_delta = 25*2*pi/360;
V = 80/3.6;

A = [ 0     1               0                   0;
      0 -(Cf+Cr)/(m*V)      (Cf+Cr)/m      (Cr*b-Cf*a)/(m*V);
      0     0               0                   1;
      0 (Cr*b-Cf*a)/(Jz*V) (Cf*a-Cr*b)/Jz -(Cr*b^2+Cf*a^2)/(Jz*V)
    ];

B1 = [0;  Cf/m;                 0;  Cf*a/Jz];
B2 = [0; (Cr*b-Cf*a)/(m*V)-V;   0;  -(Cr*b^2+Cf*a^2)/(Jz*V)]; 
B = [B1, B2];
C = eye(4);
D = zeros(4,2);

% FEEDBACK POLE PLACEMENT
K_feedback_pole = place(A, B1, [-1, -0.5, -1.1, -0.6]*10);
% [eig(A), eig(A-B1*K_feedback_pole)];

% FEEDBACK LQR
Q = diag([1, 1, 1, 1]);
R = 1;
K_feedback_LQR = lqr(A, B1, Q, R);

% FEEDFORWARD
FeedForward = @(K) m*V^2/L * (b/Cf - a/Cr + a/Cr*K(3)) + L - b*K(3);
K_feedforward = FeedForward(K_feedback_LQR);
% INTEGRATIVE
K_integrative = [1 1];

% K_fb = K_feedback_LQR;
% K_ff = K_feedforward;
% K_i = K_integrative;

K_fb = zeros(1,4);
K_ff = 0;
K_i = [0 0];

x0 = 0;
ff_delay = 0;
trajectory_type = 1;
return;

%% ---------------- TEST 0 ----------------
disp('---- TEST 0.1 - SYSTEM BEHAVIOUR ----')
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;


trajectory_type = 1;
K_fb = zeros(1,4);
K_ff = 0;
K_i = [0 0];

x0 = 0;
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);
x0 = 10;
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

close all
plotlegend = {'x0 = 0','x0 = 10'};
PLOT_TRAJ(t, V, Kl, traj, N, plotlegend, 'TEST 0.1 TRAJ')
PLOT_U(t, u1v, u2v, N, plotlegend, 'TEST 0.1 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 0.1 E')

%%
disp('---- TEST 0.2 - STEP RESPONSE ----')
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;

trajectory_type = 5;
V_prec = V; V=0; % psi = 0
max_delta_prec = max_delta;
max_delta = 1;
K_ff = 1; % delta=1
x0 = 0; 
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);
close all
plotlegend = {};
PLOT_U(t, u1v*pi/180, u2v, N, plotlegend, 'TEST 0.2 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 0.2 E')
V = V_prec;
max_delta = max_delta_prec;

%%
disp('---- TEST 0.3 - LRQ vary R and Q ----')
K_LQR = [];

% FEEDBACK LQR
Q = diag([1, 1, 1, 1]);
R = 1;
K_LQR(:,1) = lqr(A, B1, Q, R);

Q = diag([2, 2, 1, 1]);
R = 1;
K_LQR(:,2) = lqr(A, B1, Q, R);

Q = diag([10, 10, 1, 1]);
R = 1;
K_LQR(:,3) = lqr(A, B1, Q, R);

Q = diag([0.1, 0.1, 1, 1]);
R = 1;
K_LQR(:,4) = lqr(A, B1, Q, R);

K_LQR

plotlegend = {};
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;
for i = 1:4
    K_fb = K_LQR(:,i)';
    K_ff = FeedForward(K_fb);
    sim(model_name)
    [e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);
    plotlegend{i} = sprintf('K%d', i);
end
close all
PLOT_U(t, u1v, u2v, N, plotlegend, 'TEST 0.3 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 0.3 E')
K_fb = K_feedback_LQR;
K_ff = FeedForward(K_fb);

%% ---------------- TEST 1 ----------------
disp('---- TEST 1 - FB OFF VS ON ----')
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;

trajectory_type = 1;
x0 = 0;

% open loop
K_fb = zeros(1,4);
K_ff=FeedForward(K_fb);
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

% lqr
K_fb = K_feedback_LQR;
K_ff=FeedForward(K_fb);
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

% pole (small)
K_fb = place(A, B1, [-1, -0.5, -1.1, -0.6]);
K_ff=FeedForward(K_fb);
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

% pole (correct)
K_fb = K_feedback_pole;
K_ff=FeedForward(K_fb);
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);


close all
plotlegend = {'fb OFF','fb LQR', 'fb POLE 1', 'fb POLE 10'};
PLOT_TRAJ(t, V, Kl, traj, N, plotlegend, 'TEST 1 TRAJ')
PLOT_U(t, u1v, u2v, N, plotlegend, 'TEST 1 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 1 E')


%% ---------------- TEST 2 ----------------
disp('---- TEST 2 - FF and FB contributions ----')

trajectory_type = 1;
x0 = 0;

K_fb=K_feedback_LQR;
K_ff=FeedForward(K_fb);
sim(model_name)

close all
plotlegend = {'delta', 'FF','FB'};                           
PLOT(t, [u(:,1), FF, FB], 3, plotlegend, '', 'FF and FB')
saveas(gcf, strcat('imgs\', 'TEST 2 FF-FB', '.png'));  % Save as PNG file

%% ---------------- TEST 3 ----------------
disp('---- TEST 3 - Integrative term ----')
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;

trajectory_type = 1;

% K_fb = K_feedback_pole;
% K_fb = place(A, B1, [-1, -0.5, -1.1, -0.6]);
K_fb = zeros(1,4);
K_ff=FeedForward(K_fb);
K_integrative = [1 -1]*-1;


K_i = [0 0];
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

K_i = K_integrative*1e-7;
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

K_i = K_integrative*1e-8;
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

close all
plotlegend = {'integ off','integ 1','integ 2','integ 3'};
PLOT_U(t, u1v, u2v, N, plotlegend, 'TEST 3 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 3 E')
PLOT_TRAJ(t, V, Kl, traj, N, plotlegend, 'TEST 3 TRAJ')



%% ---------------- TEST 4 ----------------
disp('---- TEST 4 - vary Vx ----')
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;
V_prec=V;A_prec=A;B_prec=B;
trajectory_type = 1;


plotlegend = {};
Vv = [80, 160, 220]/3.6;
L = length(Vv);
for i = 1:L
    V = Vv(i);
    plotlegend{i} = sprintf('V = %d', V*3.6);
    A = [ 0     1               0                   0;
      0 -(Cf+Cr)/(m*V)      (Cf+Cr)/m      (Cr*b-Cf*a)/(m*V);
      0     0               0                   1;
      0 (Cr*b-Cf*a)/(Jz*V) (Cf*a-Cr*b)/Jz -(Cr*b^2+Cf*a^2)/(Jz*V)
    ];
    B1 = [0;  Cf/m;                 0;  Cf*a/Jz];
    B2 = [0; (Cr*b-Cf*a)/(m*V)-V;   0;  -(Cr*b^2+Cf*a^2)/(Jz*V)]; 
    B = [B1, B2];
    
    K_fb = lqr(A, B1, Q, R);
    K_ff = FeedForward(K_fb);

    sim(model_name)
    [e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);
    
    fprintf('V = %3d, K_ff = %3.1f, K_fb = [ %3.2f, %3.2f, %3.2f, %3.2f ]\n', V*3.6,K_ff,K_fb(1),K_fb(2),K_fb(3),K_fb(4))
end

close all
PLOT_U(t, u1v, u2v, N, plotlegend, 'TEST 4 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 4 E')

V = V_prec; A = A_prec; B = B_prec;

%% ---------------- TEST 5 ----------------
disp('---- TEST 5 - Trajectories ----')


K_fb = K_feedback_LQR;
K_ff = FeedForward(K_fb);
plotlegend = {'curve', 'circle', 'skid pad', 'obst avoid'};
close all
for i = 1:4
    trajectory_type = i;
    e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;
    sim(model_name)
    [e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

    PLOT_TRAJ(t, V, Kl, traj, 1, {plotlegend{i}}, sprintf('TEST 5.%d TRAJ', i))

end

%% ---------------- TEST 6 ----------------
disp('---- TEST 6 - Cf and Cr vary ----')
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;
A_prec=A;B_prec=B;Cf_prec=Cf;Cr_prec=Cr;a_prec=a;b_prec=b;
trajectory_type = 1;

plotlegend = {};
Cf_v = [Cf_prec, Cr_prec, (Cf_prec+Cr_prec)/2];
Cr_v = [Cr_prec, Cf_prec, (Cf_prec+Cr_prec)/2];
a_v = [a_prec, b_prec, (a_prec+b_prec)/2];
b_v = [b_prec, a_prec, (a_prec+b_prec)/2];
L = length(Cf_v);
for i = 1:L
    Cf = Cf_v(i);
    Cr = Cr_v(i);
    a = a_v(i);
    b = b_v(i);
    
    A = [ 0     1               0                   0;
      0 -(Cf+Cr)/(m*V)      (Cf+Cr)/m      (Cr*b-Cf*a)/(m*V);
      0     0               0                   1;
      0 (Cr*b-Cf*a)/(Jz*V) (Cf*a-Cr*b)/Jz -(Cr*b^2+Cf*a^2)/(Jz*V)
    ];
    B1 = [0;  Cf/m;                 0;  Cf*a/Jz];
    B2 = [0; (Cr*b-Cf*a)/(m*V)-V;   0;  -(Cr*b^2+Cf*a^2)/(Jz*V)]; 
    B = [B1, B2];
    
    K_fb = lqr(A, B1, Q, R);
    K_ff = FeedForward(K_fb);

    sim(model_name)
    [e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);
    if Cr*b-Cf*a>0
        disp('============ understeering ============')
        plotlegend{i} = 'understeering';
    elseif Cf*a-Cr*b==0
        disp('============ neutral vehicle ============')
        plotlegend{i} = 'neutral';
    else
        disp('============ oversteering ============')
        plotlegend{i} = 'oversteering';
    end
    fprintf('Cf = %d, Cr = %d, a = %.1f, b = %.1f\n', Cf, Cr, a, b)
    % fprintf('K_ff = %3.1f, K_fb = [ %3.2f, %3.2f, %3.2f, %3.2f ]\n',K_ff,K_fb(1),K_fb(2),K_fb(3),K_fb(4))
end

close all
PLOT_U(t, u1v, u2v, N, plotlegend, 'TEST 6 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 6 E')

A = A_prec; B = B_prec;
Cf = Cf_prec; Cr = Cr_prec;
a = a_prec; b = b_prec;

%% ---------------- TEST 7 ----------------
disp('---- TEST 7 - FF delay ----')
e1v=[];e2v=[];e3v=[];e4v=[];u1v=[];u2v=[];traj=[];N=0;

K_fb = K_feedback_LQR;
K_ff = FeedForward(K_fb);

ff_delay = 0;
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

ff_delay = 0.01;
sim(model_name)
[e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N);

close all
plotlegend = {'FF delay OFF', 'FF delay ON'};
PLOT_U(t, u1v, u2v, N, plotlegend, 'TEST 7 U')
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'TEST 7 E')
% PLOT_TRAJ(t, V, Kl, traj, N, plotlegend, 'TEST 7 TRAJ')

ff_delay = 0; 
%% prova PLOT
close all
PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, 'prova')

%% ---------------------------- FUNCTIONS ---------------------------------
function [e1v, e2v, e3v, e4v, u1v, u2v, traj, N] = take_res(e, u, Var_trajectory, e1v, e2v, e3v, e4v, u1v, u2v, traj, N)
    N = N+1;
    e1v(:,N) = e(:,1);
    e2v(:,N) = e(:,2);
    e3v(:,N) = e(:,3);
    e4v(:,N) = e(:,4);
    u1v(:,N) = u(:,1);
    u2v(:,N) = u(:,2);
    traj(N, :, :) = Var_trajectory;
end

function PLOT(t, x, L, plotlegend, yaxis, nome_fig)
    linewidth = 1.4;
    plotcol = {'b','r','g','m','k'};
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
    xlabel('time [s]');
    ylabel(yaxis);
    xlim([t(1), t(end)]);
    title(nome_fig, Interpreter="latex", FontSize=18)
    if isempty(plotlegend) == 0
        legend(plotlegend{1:L})
    end
    end

function PLOT_E(t, e1v, e2v, e3v, e4v, N, plotlegend, nome_fig )
    figure('Name','State space','NumberTitle','off','PaperType','A4')
    set(gcf, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    subplot(2,2,1)
    PLOT(t, e1v, N, plotlegend,'[m]','$\textbf{e1 = y - $y_{ref}$}$')
    subplot(2,2,2)
    PLOT(t, e2v, N, plotlegend,'[m/s]','$\textbf{e2 = $\dot{e1}$}$')
    subplot(2,2,3)
    PLOT(t, e3v, N, plotlegend,'[rad]','$\textbf{e3 = $\psi$ - $\psi_{des}$}$')
    subplot(2,2,4)
    PLOT(t, e4v, N, plotlegend,'[rad/s]','$\textbf{e4 = $\dot{e3}$}$')
    saveas(gcf, strcat('imgs\', nome_fig, '.png'));  % Save as PNG file

end

function PLOT_U(t, u1v, u2v, N, plotlegend, nome_fig)
    figure('Name','Input delta','NumberTitle','off','PaperType','A4')
    subplot(2,1,1)
    PLOT(t, u1v*360/2/pi, N, plotlegend,'[deg]','$\textbf{u1 $\delta$}$')
    subplot(2,1,2)
    % figure('Name','Input psi','NumberTitle','off','PaperType','A4')
    PLOT(t, u2v, N, plotlegend,'[rad/s]','$\textbf{u2 $\psi$}$')
    saveas(gcf, strcat('imgs\', nome_fig, '.png'));  % Save as PNG file

end

function PLOT_TRAJ(t, V, Kl, traj, N, plotlegend, nome_fig)
    linewidth = 1.4;
    plotcol = {'b','r','g','m','k'};
    dt_sim = 1e-3;
    th = V*dt_sim*Kl;
    cum_th = cumsum(th);
    ds = V*dt_sim;
    x = ds * cos(cum_th);
    y = ds * sin(cum_th);
    cum_x = cumsum(x);
    cum_y = cumsum(y);
    
    figure('Name','Curvature','NumberTitle','off','PaperType','A4')
    plot(t, Kl, 'b', LineWidth=linewidth)
    ylabel('[1/m]')
    xlabel('time [s]')
    xlim([t(1), t(end)])
    title('Curvature Kl');
    saveas(gcf, strcat('imgs\', nome_fig, ' Kl.png'));  % Save as PNG file
    
    figure('Name','Trajectory','NumberTitle','off','PaperType','A4')
    ppp = plot(cum_x, cum_y, '--k', LineWidth=2.5);
    ppp.Color(4) = 0.5;
    for i=1:N
        hold on
        ppp = plot(traj(i,:,1), traj(i,:,2), plotcol{i}, LineWidth=linewidth);
        ppp.Color(4) = 0.5;
    end
    x_max = max(max(abs(traj(:,:,1))));
    y_max = max(max(abs(traj(:,:,2))));
    % xlim([min(min(traj(:,:,1))) - x_max*0.1, max(max(traj(:,:,1))) + x_max*0.1 ])
    % ylim([min(min(traj(:,:,2))) - y_max*0.1, max(max(traj(:,:,2))) + y_max*0.1 ])
    xlabel('x [m]')
    ylabel('y [m]')
    plot(0,0,'*r','Markersize', 15)
    plotlegend = {'Traj desired', plotlegend{1:N}};
    legend(plotlegend)
    title('Trajectories');
    saveas(gcf, strcat('imgs\', nome_fig, '.png'));  % Save as PNG file

end