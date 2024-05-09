%% MVM course
close all; clc; clear; 
warning off

%% Load tyre data (Pacejka model coefficients)
% File with the tyre data
WheelFile = 'Tyre225_50_17_Comb';      
% Wheel Initialization
pacn = [];
eval(['[R_in,largh,AR,Kp,Cp,Jns,lon_r,lat_r,sigma_lon_r,sigma_lat_r,eps_f,V_low,kvl0,camb_0F,camb_0R,toe_0F,toe_0R,KF,Pacejka,f0_rot,K_rot,Tresh_slip]=' WheelFile ';'])
pacn = struct2cell(Pacejka);
for ii = 1:size(pacn)
    Pace(ii) = pacn{ii};
end
Pacn = Pace';     % array containing the coefficients of Pacejka tyre model

% Simulink model name
simulink_model_name = 'Tyre_Pac96';

% friction coefficient setting
mu_slope = 0;
mu0 = 1;
gamma0 = 0;
gamma_slope = 0;

% -------------------------Plot settings --------------------------
F_Size = 14; % FontSize
plotcol = {'--k','-r','-.b',':g','.m','-k'};

%% ------------------------- Test 1 --------------------------
% Plot settings
plotleg = {'F_Z = 1500 N','F_Z = 3000 N','F_Z = 4500 N','F_Z = 6000 N'};
% Plot (Fx - s) 
nome_fig = 'Fx-s';
figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
hold all; grid on

% Setting data for  Test  1
alfa0 = 0;
alfa_slope = 0;
Fz_vet = [1500 3000 4500 6000];
Fz_slope = 0;
s0 = -1;
s_slope = 0.2;

disp('Test (Fx-s) ...')
% Ciclo for per l'esecuzione del numero di simulazioni richieste
for cont1 = 1: length(Fz_vet)
    Fz0 = Fz_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
    plot(s(:,2),Fx(:,2),plotcol{cont1},'Linewidth',2); 
    cont1;
end
% Graph settings and saving as .fig and .png
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('s[-]'); ylabel('F_X [N]')
% print -dpng Immagini_Pac96\Fx_s.png
% saveas(gcf, 'Immagini_Pac96\Fx_s', 'fig')
% set(gcf,'Position',[150 38 985 667]);

%% ------------------------- Test  2 --------------------------
% Plot settings
% Plot (Fy - alfa) 
plotleg = {'F_Z = 1500 N','F_Z = 3000 N','F_Z = 4500 N','F_Z = 6000 N'};
nome_fig = 'Fy-alfa';
h_Fy_ = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz Fy';
h_Mz_Fy = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz alpha';
h_Mz_alpha = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on


% Setting data for  Test  2
alfa0 = -25;
alfa_slope = 5;
Fz_vet = [1500 3000 4500 6000];
s0 = 0;
s_slope = 0;

disp('Test (Fy - alfa) ...')
for cont1 = 1: length(Fz_vet)
    Fz0 = Fz_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
    %     sim PneumaticoCorretto_mod
    figure(h_Fy_)
    plot(alfa(:,2),Fy(:,2),plotcol{cont1},'Linewidth',2);
    figure(h_Mz_Fy)
    plot(Mz(:,2),Fy(:,2),'Linewidth',2);
    figure(h_Mz_alpha)
    plot(alfa(:,2),Mz(:,2),'Linewidth',2);
    cont1;
end
figure(h_Fy_)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [deg]'); ylabel('F_Y [N]')
xlim([-1 1]*25)

figure(h_Mz_Fy)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('Mz [Nm]'); ylabel('F_Y [N]')

figure(h_Mz_alpha)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [°]'); ylabel('M_Z [Nm]')
xlim([-1 1]*25)

figure(h_Fy_)
% print -dpng Immagini_Pac96\Fy_alpha.png
% saveas(gcf, 'Immagini_Pac96\Fy_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
xlim([-1 1]*25)

figure(h_Mz_Fy)
% print -dpng Immagini_Pac96\Mz_Fy.png
% saveas(gcf, 'Immagini_Pac96\Mz_Fy', 'fig')
% set(gcf,'Position',[150 38 985 667]);

figure(h_Mz_alpha)
% print -dpng Immagini_Pac96\Mz_alpha.png
% saveas(gcf, 'Immagini_Pac96\Mz_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
xlim([-1 1]*25)

%% ------------------------- Test  1 mu --------------------------
% Plot settings
% Plot (Fy - alfa) par mu 
% plotleg = {'F_Z = 1500 N','F_Z = 3000 N','F_Z = 4500 N','F_Z = 6000 N'};
nome_fig = 'Fx-slip';
h_Fx_ = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

% Setting data for  Test  2
Fz0 = 4000;
Fz_slope = 0;
mu_slope = 0;
mu0 = 1;
alfa0 = 0;
alfa_slope = 0;
mu_vet = [0.1 0.3 0.6 0.8 1];
s0 = -1;
s_slope = 0.2;
clear plotleg_auto
disp('Test (Fx-slip mu) ...')
for cont1 = 1: length(mu_vet)
    mu0 = mu_vet(cont1);
    plotleg_auto{cont1} = ['\mu = ',num2str(mu0)];
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
figure(h_Fx_)
plot(s(:,2),Fx(:,2),plotcol{cont1},'Linewidth',2); 
     cont1;
end
figure(h_Fx_)
legend(plotleg_auto,'location','southeast')
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\kappa [-]'); ylabel('F_X [N]')
% xlim([-0.2 0.2])

% print -dpng Immagini_Pac96\Fx_mu.png
% saveas(gcf, 'Immagini_Pac96\Fx_mu', 'fig')
% set(gcf,'Position',[150 38 985 667]);

% error('exit')

%% ------------------------- Test  2 mu --------------------------
% Plot settings
% Plot (Fy - alfa) 
% plotleg = {'F_Z = 1500 N','F_Z = 3000 N','F_Z = 4500 N','F_Z = 6000 N'};
nome_fig = 'Fy-alfa';
h_Fy_ = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz Fy';
h_Mz_Fy = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz alpha';
h_Mz_alpha = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on


% Setting data for  Test  2
alfa0 = -25;
alfa_slope = 5;
Fz0 = 4000;
Fz_slope = 0;
s0 = 0;
s_slope = 0;
mu_vet = [0.1 0.3 0.6 0.8 1];
disp('Test (Fy-alfa Mz-Fy Mz-alpha mi) ...')
for cont1 = 1: length(mu_vet)
    mu0 = mu_vet(cont1);
    plotleg_auto{cont1} = ['\mu = ',num2str(mu0)];
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
    
    figure(h_Fy_)
    plot(alfa(:,2),Fy(:,2),plotcol{cont1},'Linewidth',2);
    figure(h_Mz_Fy)
    plot(Mz(:,2),Fy(:,2),'Linewidth',2);
    figure(h_Mz_alpha)
    plot(alfa(:,2),Mz(:,2),'Linewidth',2);
    cont1;
end
figure(h_Fy_)
legend(plotleg_auto,'location','southeast')
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [deg]'); ylabel('F_Y [N]')
xlim([-1 1]*25)

figure(h_Mz_Fy)
legend(plotleg_auto,'location','southeast')
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('Mz [Nm]'); ylabel('F_Y [N]')

figure(h_Mz_alpha)
legend(plotleg_auto,'location','southeast')
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [°]'); ylabel('M_Z [Nm]')
xlim([-1 1]*25)

figure(h_Fy_)
% print -dpng Immagini_Pac96\Fy_alpha_mu.png
% saveas(gcf, 'Immagini_Pac96\Fy_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
xlim([-1 1]*25)

figure(h_Mz_Fy)
% print -dpng Immagini_Pac96\Mz_Fy_mu.png
% saveas(gcf, 'Immagini_Pac96\Mz_Fy', 'fig')
% set(gcf,'Position',[150 38 985 667]);

figure(h_Mz_alpha)
% print -dpng Immagini_Pac96\Mz_alpha_mu.png
% saveas(gcf, 'Immagini_Pac96\Mz_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
xlim([-1 1]*25)

%% ------------------------- Test  3 --------------------------
% Plot settings
% plotleg = {'\alpha = 0.1°','\alpha = 5°','\alpha = 10°'};
% Plot (Fx - s) 
nome_fig = 'Cornering Stiffness';
figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
hold all; grid on
% Setting data for  Test 
alfa_vet = [0.1];
alfa_slope = 0;
Fz0 = 0;
Fz_slope = 1000;
s0 = 0;
s_slope = 0;

disp('Test (Cornering Stiffness) ...')
for cont1 = 1: length(alfa_vet)
    alfa0 = alfa_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
%     plot(Fz(:,2),Fy(:,2)./alfa(:,2),plotcol{cont1},'Linewidth',2); 
    plot(Fz(:,2),-C_F(:,3),'Linewidth',2); 
    cont1;
end
% legend(plotleg)

set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('F_Z [N]'); ylabel('K_{\alpha} [N/rad]')
% print -dpng Immagini_Pac96\kalfa_Fz.png
% saveas(gcf, 'Immagini_Pac96\kalfa_Fz', 'fig')
% set(gcf,'Position',[150 38 985 667]);

%% ------------------------- Test  4 --------------------------
% Plot settings
plotleg = {'\alpha = 2°','\alpha = 4°','\alpha = 6°','\alpha = 8°'};
% Plot (Fx - s) 
nome_fig = 'Ellipse of adherence';
figure1 = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on
% Setting data for  Test  2
alfa_vet = 2:2:8;
alfa_slope = 0;
Fz0 = 5000;
Fz_slope = 0;
s0 = -0.3;
s_slope = 0.06;

disp('Test (Ellipse of adherence) ...')
for cont1 = 1: length(alfa_vet)
    alfa0 = alfa_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
    plot(-Fx(:,2),Fy(:,2),plotcol{cont1},'Linewidth',2); 
    cont1;
end
legend(plotleg,'Location','NorthEast')

% Setting data for  Test  4b
alfa0 = 0;
alfa_slope = 0.8;
Fz0 = 5000;
Fz_slope = 0;
s_vet = [-0.1 -0.05 -0.025 0.025 0.05 0.1];
s_slope = 0;

for cont1 = 1: length(s_vet)
    s0 = s_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
    plot(-Fx(:,2),Fy(:,2),'--k','Linewidth',1); 
    cont1;
end

set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('F_X [N]'); ylabel('F_Y [N]')
ylim([0 5500]);

% Create textbox
annotation1 = annotation(...
  figure1,'textbox',...
  'Position',[0.6294 0.1349 0.07614 0.04048],...
  'BackgroundColor',[1 1 1],...
  'LineStyle','none',...
  'FitHeightToText','off',...
  'FontName','Times New Roman',...
  'FontSize',14,...
  'String',{'s=0.025'});
 
% Create textbox
annotation2 = annotation(...
  figure1,'textbox',...
  'Position',[0.7439 0.1344 0.06497 0.04048],...
  'BackgroundColor',[1 1 1],...
    'LineStyle','none',...
  'FitHeightToText','off',...
  'FontName','Times New Roman',...
  'FontSize',14,...
  'String',{'s=0.05'});
 
% Create textbox
annotation3 = annotation(...
  figure1,'textbox',...
  'Position',[0.8319 0.1339 0.05239 0.04048],...
  'BackgroundColor',[1 1 1],...
    'LineStyle','none',...
  'FitHeightToText','off',...
  'FontName','Times New Roman',...
  'FontSize',14,...
  'String',{'s=0.1'});
 
% Create textbox
annotation4 = annotation(...
  figure1,'textbox',...
  'Position',[0.3127 0.1364 0.07888 0.04048],...
  'BackgroundColor',[1 1 1],...
    'LineStyle','none',...
  'FitHeightToText','off',...
  'FontName','Times New Roman',...
  'FontSize',14,...
  'String',{'s=-0.025'});
 
% Create textbox
annotation5 = annotation(...
  figure1,'textbox',...
  'Position',[0.2153 0.1999 0.07249 0.04048],...
  'BackgroundColor',[1 1 1],...
    'LineStyle','none',...
  'FitHeightToText','off',...
  'FontName','Times New Roman',...
  'FontSize',14,...
  'String',{'s=-0.05'});
 
% Create textbox
annotation6 = annotation(...
  figure1,'textbox',...
  'Position',[0.1421 0.1344 0.06467 0.04048],...
  'BackgroundColor',[1 1 1],...
    'LineStyle','none',...
  'FitHeightToText','off',...
  'FontName','Times New Roman',...
  'FontSize',14,...
  'String',{'s=-0.1'});
 
% print -dpng Immagini_Pac96\Ellisse.png
% saveas(gcf, 'Immagini_Pac96\Ellisse', 'fig')
% set(gcf,'Position',[150 38 985 667]);

%% ------------------------- Test  5 --------------------------
% Plot settings
plotleg = {'Fx','Fy(\alpha=5°)'};
% Plot (Fx - s) 
nome_fig = 'Combined load';
figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
hold all; grid on
% Setting data for  Test 
alfa_vet = [5];
alfa_slope = 0;
Fz0 = 4500;
Fz_slope = 0;
s0 = -1;
s_slope = 0.2;

disp('Test (Combined load) ...')
for cont1 = 1: length(alfa_vet)
    alfa0 = alfa_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
    plot(s(:,2),Fx(:,2),plotcol{cont1},'Linewidth',2); 
    plot(s(:,2),Fy(:,2),plotcol{cont1+2},'Linewidth',2); 
    cont1;
end
legend(plotleg)

set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('s [-]'); ylabel('F [N]')
% print -dpng Immagini_Pac96\Fx_Fy_Comb.png
% saveas(gcf, 'Immagini_Pac96\Fx_Fy_Comb', 'fig')
% set(gcf,'Position',[150 38 985 667]);


%% ------------------------- Test  6 --------------------------
% Plot settings
plotleg = {'s=0.025','s=0.05','s=0.075','s=0.1'};
% Plot (Fx - alpha) 
nome_fig = 'Combined load: Fx (alpha)';
figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
hold all; grid on
% Setting data for  Test 
alfa0 = -25;
alfa_slope = 5;
Fz0 = 4500;
Fz_slope = 0;
s_vet = 0.025:0.025:0.1;
s_slope = 0;

disp('Test (Combined load: Fx-alpha)...')
for cont1 = 1: length(s_vet)
    s0 = s_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
    plot(alfa(:,2),Fx(:,2),plotcol{cont1},'Linewidth',2); 
    cont1;
end
legend(plotleg)

set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [deg]'); ylabel('F_X [N]')
% print -dpng Immagini_Pac96\Fx_Comb.png
% saveas(gcf, 'Immagini_Pac96\Fx_Comb', 'fig')
% set(gcf,'Position',[150 38 985 667]);

%% ------------------------- Test  7 --------------------------
% Plot settings
plotleg = {'\alpha=2°','\alpha=4°','\alpha=6°','\alpha=8°'};
% Plot (Fy - s) 
nome_fig = 'Combined load: Fy (s)';
figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
hold all; grid on
% Setting data for  Test 
alfa_vet = 2:2:8;
alfa_slope = 0;
Fz0 = 4500;
Fz_slope = 0;
s0 = -1;
s_slope = 0.2;

disp('Test (Combined load: Fy-s)...')
for cont1 = 1: length(alfa_vet)
    alfa0 = alfa_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
    plot(s(:,2),Fy(:,2),plotcol{cont1},'Linewidth',2); 
    cont1;
end
legend(plotleg)

set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('s [-]'); ylabel('F_Y [N]')
% print -dpng Immagini_Pac96\Fy_Comb.png
% saveas(gcf, 'Immagini_Pac96\Fy_Comb', 'fig')
% set(gcf,'Position',[150 38 985 667]);

%% ------------------------- Test  8 - Combined load: Fx (alpha) , Fy (alpha) --------------------------
% Plot settings
plotleg = {'s=0','s=0.05','s=0.01','s=0.2','s=1'};
% Plot (Fx - alpha) 
nome_fig = 'Combined load: Fx (alpha)';
h_Fx_alpha = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Combined load: Fy (alpha) for different slip';
h_Fy_alpha = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Combined load: Mz (alpha) for different slip';
h_Mz_alpha = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

% Setting data for  Test 
alfa0 = -25;
alfa_slope = 5;
Fz0 = 4500;
Fz_slope = 0;
% s_vet = (2.5:2.5:10)/100;
s_vet = [0 5 10 20 100]./100;
s_slope = 0;

disp('Test (Combined load: Fx-alpha , Fy-alpha)...')
for cont1 = 1: length(s_vet)
    s0 = s_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
    plotleg_auto{cont1} = num2str(s0*100);
%     sim PneumaticoCorretto_mod
    figure(h_Fx_alpha)
    plot(alfa(:,2),Fx(:,2),'Linewidth',2); 
    figure(h_Fy_alpha)
    plot(alfa(:,2),Fy(:,2),'Linewidth',2); 
    figure(h_Mz_alpha)
    plot(alfa(:,2),Mz(:,2),'Linewidth',2); 
    cont1;
end
figure(h_Fx_alpha)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [deg]'); ylabel('F_X [N]')

figure(h_Fy_alpha)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [deg]'); ylabel('F_Y [N]')

figure(h_Mz_alpha)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [deg]'); ylabel('M_Z [Nm]')


% print -dpng Immagini_Pac89\Fx_Comb.png
% saveas(gcf, 'Immagini_Pac89\Fx_Comb', 'fig')
% set(gcf,'Position',[150 38 985 667]);

clear plotleg_auto
error('')
%% ------------------------- Test  9 'Combined load: Fx (s) , Fy (s)'--------------------------
% Plot settings
plotleg = {'\alpha=2°','\alpha=4°','\alpha=6°','\alpha=8°'};
% Plot (Fy - s) 

nome_fig = 'Combined load: Fx (s)';
h_fig_Fxs = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Combined load: Fy (s)';
h_fig_Fys = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Combined load: Mz (s)';
h_Mz_s = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Combined load: Mz (s)';
h_Mz_Fx = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

% Setting data for  Test 
alfa_vet = 0:2:6;
% alfa_vet = 0:5:30;
% alfa_vet_pac = [2:2:6 12 30];
% alfa_vet = alfa_vet_pac;

alfa_slope = 0;
Fz0 = 4500;
Fz_slope = 0;
s0 = -1;
s_slope = 0.2;

disp('Test (Combined load: Fx-s , Fy-s)...')
for cont1 = 1: length(alfa_vet)
    alfa0 = alfa_vet(cont1);
    plotleg_auto{cont1} = ['\alpha = ',num2str(alfa0),'°'];
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
figure(h_fig_Fxs)
plot(s(:,2),Fx(:,2),'Linewidth',2);
figure(h_fig_Fys)
plot(s(:,2),Fy(:,2),'Linewidth',2);
cont1;
figure(h_Mz_Fx)
plot(Fx(:,2),Mz(:,2),'Linewidth',2);
figure(h_Mz_s)
plot(s(:,2),Mz(:,2),'Linewidth',2);
end
figure(h_fig_Fxs)
legend(plotleg_auto)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('s [-]'); ylabel('F_X [N]')

figure(h_fig_Fys)
legend(plotleg_auto)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('s [-]'); ylabel('F_Y [N]')

figure(h_Mz_s)
legend(plotleg_auto)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('slip [-]'); ylabel('M_Z [Nm]')

figure(h_Mz_Fx)
legend(plotleg_auto)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('F_X [N]'); ylabel('M_Z [Nm]')

%% ------------------------- Test  Fy-Camber par Fz --------------------------
% Effect of Camber
% Plot settings
plotleg = {'F_Z = 1500 N','F_Z = 3000 N','F_Z = 4500 N','F_Z = 6000 N'};
% Plot (Fx - s) 
nome_fig = 'Fy-camber';
figure('Name',nome_fig,'NumberTitle','off','PaperType','A4')
hold all; grid on

% Setting data for  Test  1
alfa0 = 0;
alfa_slope = 0;
gamma0 = -30;
gamma_slope = -gamma0*2/10;
Fz_vet = [1500 3000 4500 6000];
Fz_slope = 0;
s0 = 0;
s_slope = 0;

disp('Test (Fy-camber)...')
% Ciclo for per l'esecuzione del numero di simulazioni richieste
for cont1 = 1: length(Fz_vet)
    Fz0 = Fz_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
%     sim PneumaticoCorretto_mod
    plot(gamma(:,2),Fy(:,2),plotcol{cont1},'Linewidth',2); 
    cont1;
end
% Graph settings and saving as .fig and .png
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\gamma[-]'); ylabel('F_Y [N]')
% print -dpng Immagini_Pac96\Fx_s.png
% saveas(gcf, 'Immagini_Pac96\Fx_s', 'fig')
% set(gcf,'Position',[150 38 985 667]);

%% ------------------------- Test Fy-alfa par gamma --------------------------
% Plot settings
% Plot (Fy - alfa) 
plotleg = {'\gamma = -5°','\gamma = 0°','\gamma = 5°','\gamma = 10°'};
nome_fig = 'Fy-alfa';
h_Fy_ = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz Fy';
h_Mz_Fy = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz alpha';
h_Mz_alpha = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

% Setting data for  Test  2
alfa0 = -25;
alfa_slope = 5;
gamma_vet = [-5 0 5];
gamma_slope = 0;
s0 = 0;
s_slope = 0;
Fz0 = 4000;
Fz_slope = 0;

disp('Test (Fy-alfa, Fy-Mz, Mz-alfa par gamma)...')
for cont1 = 1: length(gamma_vet)
    gamma0 = gamma_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
    %     sim PneumaticoCorretto_mod
    figure(h_Fy_)
    plot(alfa(:,2),Fy(:,2),plotcol{cont1},'Linewidth',2);
    figure(h_Mz_Fy)
    plot(Mz(:,2),Fy(:,2),'Linewidth',2);
    figure(h_Mz_alpha)
    plot(alfa(:,2),Mz(:,2),'Linewidth',2);
    cont1;
end
warning off
figure(h_Fy_)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [deg]'); ylabel('F_Y [N]')
xlim([-1 1]*25)

figure(h_Mz_Fy)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('Mz [Nm]'); ylabel('F_Y [N]')

figure(h_Mz_alpha)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\alpha [°]'); ylabel('M_Z [Nm]')
xlim([-1 1]*25)

figure(h_Fy_)
% print -dpng Immagini_Pac96\Fy_alpha.png
% saveas(gcf, 'Immagini_Pac96\Fy_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
xlim([-1 1]*25)

figure(h_Mz_Fy)
% print -dpng Immagini_Pac96\Mz_Fy.png
% saveas(gcf, 'Immagini_Pac96\Mz_Fy', 'fig')
% set(gcf,'Position',[150 38 985 667]);

figure(h_Mz_alpha)
% print -dpng Immagini_Pac96\Mz_alpha.png
% saveas(gcf, 'Immagini_Pac96\Mz_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
xlim([-1 1]*25)
%% ------------------------- Test  300 --------------------------
% Plot settings
% Plot (Fy - alfa) 
plotleg = {'F_Z = 1500 N','F_Z = 3000 N','F_Z = 4500 N','F_Z = 6000 N'};
nome_fig = 'Fy-gamma';
h_Fy_ = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz Fy';
h_Mz_Fy = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

nome_fig = 'Mz gamma';
h_Mz_alpha = figure('Name',nome_fig,'NumberTitle','off','PaperType','A4');
hold all; grid on

% Setting data for  Test  2
alfa0 = 0;
alfa_slope = 0;
gamma0 = -20;
gamma_slope = 4;
Fz_vet = [1500 3000 4500 6000];
Fz_slope = 0;
s0 = 0;
s_slope = 0;

disp('Test (Fy-gamma, Fy-Mz, Mz-gamma par Fz) ...')
for cont1 = 1: length(Fz_vet)
    Fz0 = Fz_vet(cont1);
    % Run Simulink model
    string_sim = strcat({'sim '},{simulink_model_name});
    eval(string_sim{1})
    %     sim PneumaticoCorretto_mod
    figure(h_Fy_)
    plot(gamma(:,2),Fy(:,2),plotcol{cont1},'Linewidth',2);
    figure(h_Mz_Fy)
    plot(Mz(:,2),Fy(:,2),'Linewidth',2);
    figure(h_Mz_alpha)
    plot(gamma(:,2),Mz(:,2),'Linewidth',2);
    cont1;
end
figure(h_Fy_)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\gamma [deg]'); ylabel('F_Y [N]')
% xlim([-1 1]*25)

figure(h_Mz_Fy)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('Mz [Nm]'); ylabel('F_Y [N]')

figure(h_Mz_alpha)
legend(plotleg)
set(gca,'FontName','Times New Roman','FontSize',F_Size)
xlabel('\gamma [°]'); ylabel('M_Z [Nm]')
% xlim([-1 1]*25)

figure(h_Fy_)
% print -dpng Immagini_Pac96\Fy_alpha.png
% saveas(gcf, 'Immagini_Pac96\Fy_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
% xlim([-1 1]*25)

figure(h_Mz_Fy)
% print -dpng Immagini_Pac96\Mz_Fy.png
% saveas(gcf, 'Immagini_Pac96\Mz_Fy', 'fig')
% set(gcf,'Position',[150 38 985 667]);

figure(h_Mz_alpha)
% print -dpng Immagini_Pac96\Mz_alpha.png
% saveas(gcf, 'Immagini_Pac96\Mz_alpha', 'fig')
% set(gcf,'Position',[150 38 985 667]);
% xlim([-1 1]*25)
warning on

