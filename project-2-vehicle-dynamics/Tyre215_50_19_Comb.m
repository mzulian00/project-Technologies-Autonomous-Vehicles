function [Pacejka]=Tyre
%%%%
% This script initializes the coefficients for Pacejka'96 tyre model (the
% rolling resistance contribution is not included in this formulation)
% Tyre code: 225/50 19
%%%%

% Tyre model
% Pacejka96 coefficients
%nominal load Fnom
Pacejka.FZ0_  = 4000;
% scaling factors
% pure slip
Pacejka.LFZ0 = 1;
Pacejka.LmuV = 1;
Pacejka.LKXK = 1;
Pacejka.LKYA = 1;
Pacejka.LCX_ = 1;
Pacejka.LCY_ = 1;
Pacejka.LEX_ = 1;
Pacejka.LEY_ = 1;
Pacejka.LHX_ = 0;   % = 1 horizontal shift (X)
Pacejka.LHY_ = 0;   % = 1 horizontal shift (Y)
Pacejka.LVX_ = 0;   % = 1 vertical shift (X)
Pacejka.LVY_ = 0;   % = 1 vertical shift (Y)
Pacejka.LKYG = 1;
Pacejka.LKZG = 1;
Pacejka.LT__ = 1;
Pacejka.MR__ = 1;
% combined slip
Pacejka.LXA_ = 1;
Pacejka.LYK_ = 1;
Pacejka.LVYK = 1;
Pacejka.LS__ = 1;
%other
Pacejka.LCZ_ = 1;
Pacejka.LMX_ = 1;
Pacejka.LMY_ = 1;
%longitudinal
Pacejka.PCX1 = 1.6184;
Pacejka.PDX1 = 1.1787;
Pacejka.PDX2 = -0.11958;
Pacejka.PEX1 = 0.041417;
Pacejka.PEX2 = 0.095349;
Pacejka.PEX3 = 0.17843;
Pacejka.PEX4 = -0.12208;
Pacejka.PKX1 = 21.906;
Pacejka.PKX2 = 0.32899;
Pacejka.PKX3 = 0.21657;
Pacejka.PHX1 = -0.0004653;
Pacejka.PHX2 = 0.00009455;
Pacejka.PVX1 = -0.0079168;
Pacejka.PVX2 = 0.0048096;
%----Parameters from Pacejka book --------
Pacejka.RBX1 = 5.0;
Pacejka.RBX2 = 8.0;
Pacejka.RCX1 = 1.0;
Pacejka.RHX1 = 0.0;
Pacejka.RTX1 = 7.1035;
Pacejka.RTX2 = 0.0;
Pacejka.RTX3 = 0.0;
%---------------------------------------
%overturning coefficients
Pacejka.QSX1 = -0.0095549;
Pacejka.QSX2 = 0.87619;
Pacejka.QSX3 = 0.008298;
%lateral
Pacejka.PCY1 = 1.4828;
Pacejka.PDY1 = 1.0864;
Pacejka.PDY2 = -0.080154;
Pacejka.PDY3 = 3.7815;
Pacejka.PEY1 = 0.055063;
Pacejka.PEY2 = -0.80469;
Pacejka.PEY3 = -0.039428;
Pacejka.PEY4 = -6.5532;
Pacejka.PKY1 = -42.977;
Pacejka.PKY2 = 4.0656;
Pacejka.PKY3 = 0.12675;
Pacejka.PHY1 = 0.0012807;
Pacejka.PHY2 = 0.0007145;
Pacejka.PHY3 = 0.026934;
Pacejka.PVY1 = -0.011461;
Pacejka.PVY2 = 0.00673;
Pacejka.PVY3 = -0.56448;
Pacejka.PVY4 = -0.14748;
%----Parameters from Pacejka book --------
Pacejka.RBY1 = 7.0;
Pacejka.RBY2 = 2.5;
Pacejka.RBY3 = 0.0;
Pacejka.RCY1 = 1.0;
Pacejka.RHY1 = 0.02;
Pacejka.RVY1 = 0.0;
Pacejka.RVY2 = 0.0;
Pacejka.RVY3 = -0.2;
Pacejka.RVY4 = 14.0;
Pacejka.RVY5 = 1.9;
Pacejka.RVY6 = 10.0;
Pacejka.PTY1 = 0.0;
Pacejka.PTY2 = 0.0;
%---------------------------------------
% rolling resistance
Pacejka.QSY1 = 0.0;
Pacejka.QSY2 = 0.0;
% Aligning torque
Pacejka.QBZ1 = 11.312;
Pacejka.QBZ2 = -0.6980;
Pacejka.QBZ3 = 0.2136;
Pacejka.QBZ4 = -0.01457;
Pacejka.QBZ5 = 0.0881;
Pacejka.QBZ6 = 0;
Pacejka.QBZ9 = 19.097;
Pacejka.QBZZ = 0;% corresponds to QBZ10
Pacejka.QCZ1 = 1.1432;
Pacejka.QDZ1 = 0.05457;
Pacejka.QDZ2 = 0.00306;
Pacejka.QDZ3 = -0.0306;
Pacejka.QDZ4 = 1.2509;
Pacejka.QDZ6 = -0.004950;
Pacejka.QDZ7 = 0.000906;
Pacejka.QDZ8 = -0.42417;
Pacejka.QDZ9 = 0.02017;
Pacejka.QEZ1 = -4.4005;
Pacejka.QEZ2 = 0.1613;
Pacejka.QEZ3 = 0.0;
Pacejka.QEZ4 = 0.2326;
Pacejka.QEZ5 = -0.8557;
Pacejka.QHZ1 = 0.0020128;
Pacejka.QHZ2 = -0.000272;
Pacejka.QHZ3 = 0.0366;
Pacejka.QHZ4 = 0.03202;
% From Pac. book, Tab 4.2
Pacejka.SSZ1 = 0.0;
Pacejka.SSZ2 = -0.1; % from Pac Tab4.2
Pacejka.SSZ3 = -1.0; % from Pac Tab4.2
Pacejka.SSZ4 = 0;
Pacejka.QTZ1 = 0.0;