%% TAV module
close all; clc; clear; 
warning off

% Load tyre data (Pacejka model coefficients)
% File with the tyre data
WheelFile = 'Tyre215_50_19_Comb';      
% Wheel Initialization
pacn = [];
eval(['[Pacejka]=' WheelFile ';'])
pacn = struct2cell(Pacejka);
for ii = 1:size(pacn)
    Pace(ii) = pacn{ii};
end
Pacn = Pace';     % array containing the coefficients of Pacejka tyre model


