%-------------------------------------------------------------------------------
% Beginning of file Prop_CWODE_7w.m
%-------------------------------------------------------------------------------
%
%                P1        P2
%                |         |
%                |         |    
%                |         |
%                |   S&I   |    
%                |    |    |    
%                |    |    |    
% ------------------------------------------> wavelength(nm)
%                A1   A3   A2       
%                w1   w3   w2 
%               (A1) (A0) (A2)  
%
%                P1        P2
%                |         |
%      P3        |         |         P4
%      |         |   DSI   |         |
%      |    S1   |    |    |    S2   |
%      |    |    |    |    |    |    |
%      |    |    |    |    |    |    |
% ------------------------------------------> wavelength(nm)
%      A4   A6   A1   A3   A2   A7   A5    
%      d    g    a    c    b    h    f
%     (A5) (A3) (A1) (A0) (A2) (A4) (A6) 
%
% Project : Calculation of nonlinear propagation in fiber for FOPA based PSA 
%           using 7-wave model
% FileName: Prop_CWODE_7w.m
% Function: Numerical calculate the nonlinear propagation of FOPA based PIA/PSA 
%           in fiber using coupled wave ordinary differential equation of 
%           7-wave model, including:
%           13 Non-Degenerate Four-wave-mixing (NDFWM) processes
%            9     Degenerate Four-wave-mixing ( DFWM) processes
% Version : 
%   v0.0 @2015.03.20 created  by W.L. Xie 
%        Using the ordinary differential equation (ODE) of the evolution of 
%        the power and phase resolved from the field evolution of the CNLSE.
%   v0.1 @2015.03.31 released by W.L. Xie
%        Directly using the 7 complex ODE of the field evolution.
%   v0.2 @2017.01.01 released by W.L. Xie
%        Reconstruct the architecture of the code
% Describe: 
%   Simulatie signal gain and delta_beta vs. pump-signal separation
%   Simulatie signal gain and delta_beta vs. pump-signal relative phase
%-------------------------------------------------------------------------------

%% Reset and formation 
close all;
clear;
format long e;
tic;
clc;

fprintf(1,'----------------------------------------------------------------------------------\n');
fprintf(1,'Calculation of nonlinear propagation based on coupled wave ODE using 7-wave model.\n');
%-------------------------------------------------------------------------------

%% Global constant 

global c;           % 
c      = 299792458; % m/s, Speed of light in vacuum
%-------------------------------------------------------------------------------

%% Simulation mode parameters 
fprintf(1,'Configurating simulation mode and parameters... ');
                   
global SimMod;

SimMod = [];

SimMod.Wave_mode = '3W2P'; 
                 % '3W1P'
                 % '3W2P'
                 % '4W2P'
SimMod.Wave_num  =  7; % Simulation in types of wave case
                 %  3, 3-wave model 
                 %  5, 5-wave model, S&I P1 P2 P3 P4
                 %  7, 7-wave model, S&I P1 P2 P3 P4 S1 S2 (default)
                 % 11,11-wave model 
SimMod.Wave_set  =  'EqualSpaced'; % Wave allocation
                 %  'EqualSpaced', equal spaced
                 %  'Distributed', 5-wave model, 
                 %  'RandomLoacated', 7-wave model, 
SimMod.FLD_nPWR  =  1; % Simulation using Power or Field evolution
                 %  1: Field evolution (default)
                 %  0: Power/Phase evolution
SimMod.PSA_nPIA  =  1; % Simulation in PS or PI condition
                 %  1: PS, phase-sensitive manner (default)
                 %  0: PI, phase-insensitive manner
SimMod.Scan_mode = 'Signal_Power';
                 % 'Cen_Wavelength': Scan center Wvlen of excitation
                 % 'PP_Separation' : Scan pump-pump wavelength separation
                 % 'Pump_Power'    : Scan pump power
                 % 'Signal_Power'  : Scan signal power
                 % 'Relative_Phase': Scan relative phase
                 % 'Propagation_Z' : Propagate along the fiber
                 % 'Constellation' : Constellation for noise 
                 % 'Invalid Condition' 
SimMod.Scan_cond = 'Normal'; 
                 % Cen_Wavelength: Normal
                 % PP_Separation : Normal 
                 % Pump_Power    : Normal Propagate Len_Z_var  
                 % Signal_Power  : Normal 
                 % Relative_Phase: Normal 
                 % Propagation_Z : Normal Propagate Len_Z_var 
                 % Constellation : None   BPSK      QPSK 
                 % 'Invalid Condition'
SimMod.Wlofs_num =    3+1; % number of center wavelength offset, default: 200+1
SimMod.WlSep_num =    3+1; % number of pump pump wavelength separation, default: 210+1
SimMod.Power_num =    3+1; % number of pump/signal power, 250+1
SimMod.Phase_num =  420+1; % number of signal phase, default: 820+1 for 2pi, 420+1 for pi
SimMod.PropZ_num =    0+1; % number of propagate dz, default: 1011*4+1
SimMod.Len_Mul   =    0+1; % number of length multiplication 
SimMod.Mod_type  = 'None'; % 'None' 'BPSK' 'QPSK' 
SimMod.Noise_num =      0+SimMod.Phase_num; % number of additive white gaussian noise, default: 100+wPhase_num
SimMod.Freq_rang =  105*125e9; % frequency range of the scan, default: 105*125e9
SimMod.Freq_res  = SimMod.Freq_rang/(SimMod.WlSep_num-1); % Frequency resolution

fprintf(1,'[COMPLETED]\n'); 
%-------------------------------------------------------------------------------

%% Paramter of nonlinear medium 
fprintf(1,'Fabricating the NLM... ');

global NLM;     % Nonlinear medium parameters

NLM    = [];    % Class of nonlinear medium parameters

NLM.SN = 00;    % Series number of NLM

% Load NLM parameters
f_NLM( NLM.SN );   % NLM select function
                   % 00: simulation trial
                   % 04: 01+02+03, NLM standard @ ZDW=1547.50nm
                   % 07: 05+06,    NLM Spine @ ZDW=1566.00nm
                   % 11: 2014.ECOC Maxime, 6-wave
                   % 12: 2009.CLEO M. Gao, 9-wave 
                   % 13: 2012.OFC  M. Gao, 7-wave
                   % 14: 2012.OL   M. Gao, 7-wave
                   % 15: 2008.OE   Bismuth-Oxide-based NLM
                   % 16: 2011.OL   C. Lun, 3-w analytical
                   % 17: 2012.OE   C. Lun, 3-w numerical
                   % 18: 2014,JLT  D. Liu, 7-w numerical
                   % 19: Photonic Crystal Waveguide from Alfredo

fprintf(1,'[COMPLETED]\n');
%-------------------------------------------------------------------------------

%% Initial Excitation parameters 
fprintf(1,'Setting parameter of excitations ... ');

global ExcSet;
ExcSet = [];

switch SimMod.Wave_num,  
  case  3,
    ExcSet.HP = 0; % High-order pumps
    ExcSet.HS = 0; % High-order signals
  case  5,
    ExcSet.HP = 1; % High-order pumps
    ExcSet.HS = 0; % High-order signals
  case  7,
    ExcSet.HP = 1; % High-order pumps
    ExcSet.HS = 1; % High-order signals
end;
%-------------------------------------------------------------------------------

% Assign the wavelength offset of the degenerate signal and idler
if SimMod.Wlofs_num == 1,

% ExcSet.Wvlen_ofs = -10.0e-9; % m, wavelength offset 
% ExcSet.Wvlen_ofs = +00.0e-9; % m, wavelength offset
% ExcSet.Wvlen_ofs = +10.0e-9; % m, wavelength offset 
  ExcSet.Wvlen_ofs = + 3.0e-9; % m, wavelength offset 
  
% ExcSet.Lambda_s  = NLM.ZDW + ExcSet.Wvlen_ofs; % m, wavelength of signal
  ExcSet.Lambda_c  = zeros(1,1);
  
 
  
else
  %% 
  ExcSet.Wvlen_ofs = linspace(-20.0e-9,+20.0e-9,SimMod.Wlofs_num); % m, wavelength offset
% ExcSet.Wvlen_ofs = [-20.0e-9, +20.0e-9, +1.0e-9]; % m, wavelength offset

% ExcSet.Lambda_s  = 1480.0e-9 : 1.0e-9 : 1620.0e-9; % m, wavelength of signal
% ExcSet.Lambda_s  = 1527.5e-9 : 0.1e-9 : 1567.5e-9; % m, wavelength of signal
  ExcSet.Lambda_c  = zeros(1,1);
  
end;

% Assign the frequency detuning between signal and Pump1/Pump2  
if SimMod.WlSep_num == 1,

% ExcSet.Sc_PPS_fq =   4.8*125e9; % Hz, frequency detuning between S0 and P1/2
% ExcSet.Sc_PPS_fq =   9.5*125e9; % 30.0*125e9  
  ExcSet.Sc_PPS_fq =   5.625e11;  % 
  ExcSet.Dl_PPS_fq = zeros(SimMod.WlSep_num,1); % Hz, delta Lambda pump-pump separation in frequency 

  ExcSet.Sc_PPS_wl = zeros(SimMod.WlSep_num,1); 
  ExcSet.Dl_PPS_wl = zeros(SimMod.WlSep_num,1);
  
else %
  %% 
  ExcSet.Sc_PPS_fq = linspace(0*125e9,SimMod.Freq_rang,SimMod.WlSep_num); % Hz, frequency detuning between SI0 and Pump1/2
% ExcSet.Sc_PPS_fq = [ 1.0*125e9, 4.5*125e9, 7.5*125e9,10.0*125e9,...
%                     15.1*125e9,20.0*125e9,30.0*125e9,48.0*125e9];
  ExcSet.Dl_PPS_fq = zeros(SimMod.WlSep_num,1);
  
% ExcSet.Sc_PPS_wl = 0 : 0.2e-9 : 50e-9;
  ExcSet.Sc_PPS_wl = zeros(SimMod.WlSep_num,1);
  ExcSet.Dl_PPS_wl = zeros(SimMod.WlSep_num,1);
  
end;

% Assign the propagation nonlinear medium length
switch SimMod.Scan_mode, 
  case 'Cen_Wavelength', 
    ExcSet.PropZ_Len = NLM.Len; 
  case 'PP_Separation',  
    ExcSet.PropZ_Len = NLM.Len;
  case 'Pump_Power',     
    switch SimMod.Scan_cond, 
      case 'Normal',    
        ExcSet.PropZ_Len = NLM.Len*SimMod.Len_Mul;
        ExcSet.PropZ_Res = (NLM.Len*SimMod.Len_Mul-0)/(SimMod.PropZ_num-1);
      case 'Propagate', 
        ExcSet.PropZ_Len = linspace(0,NLM.Len*SimMod.Len_Mul,SimMod.PropZ_num);
        ExcSet.PropZ_Res = (NLM.Len*SimMod.Len_Mul-0)/(SimMod.PropZ_num-1);
      case 'Len_Z_var', % Not proper 
        ExcSet.PropZ_Len = linspace(0,NLM.Len*SimMod.Len_Mul,SimMod.PropZ_num);
        ExcSet.PropZ_Res = (NLM.Len*SimMod.Len_Mul-0)/(SimMod.PropZ_num-1);
      otherwise, 
    end;
  case 'Signal_Power',   
    ExcSet.PropZ_Len = NLM.Len*SimMod.Len_Mul;
  case 'Relative_Phase', 
    ExcSet.PropZ_Len = NLM.Len;
  case 'Propagation_Z',  
    switch SimMod.Scan_cond, 
      case 'Propagate',  
        ExcSet.PropZ_Len = linspace(0,NLM.Len*SimMod.Len_Mul,SimMod.PropZ_num);
        ExcSet.PropZ_Res = (NLM.Len*SimMod.Len_Mul-0)/(SimMod.PropZ_num-1);
      case 'Len_Z_var',  % Not proper
        ExcSet.PropZ_Len = linspace(0,NLM.Len*SimMod.Len_Mul,SimMod.PropZ_num);
        ExcSet.PropZ_Res = (NLM.Len*SimMod.Len_Mul-0)/(SimMod.PropZ_num-1);
      case 'Normal',    
        ExcSet.PropZ_Len = NLM.Len*SimMod.Len_Mul;
        ExcSet.PropZ_Res = (NLM.Len*SimMod.Len_Mul-0)/(SimMod.PropZ_num-1);
      otherwise, 
    end;
  case 'Constellation',  
    ExcSet.PropZ_Len = linspace(0,NLM.Len*SimMod.Len_Mul,SimMod.PropZ_num);
    ExcSet.PropZ_Res = (NLM.Len*SimMod.Len_Mul-0)/(SimMod.PropZ_num-1);
  otherwise,
end;
  
fprintf(1,'[COMPLETED]\n');
%-------------------------------------------------------------------------------

%% Initial field of all the waves 
fprintf(1,'Initializing excitation field of waves ... ');

global Exc;

% Exc.Pwr; % W, initial power of all the waves (@Z=0)
Exc.Pwr      = zeros(SimMod.Power_num, SimMod.Wave_num); % Initial the initial power of waves
switch SimMod.Wave_mode 
  case '3W1P'  
    if SimMod.Power_num == 1,
      Exc.Pwr(:,1) = 1.0;  % W, Initial power of Pump1  (P1), non-degenerate
      Exc.Pwr(:,2) = 1.0;  % W, Initial power of Pump2  (P2), non-degenerate
      Exc.Pwr(:,3) = 1.0e-6;  % W, Initial power of Signal0(S0), degenerate signal & Idler 
      Exc.Pwr(:,4) = 0.0e-3;  % W, Initial power of Pump3  (P3), (@Z=0)
      Exc.Pwr(:,5) = 0.0e-3;  % W, Initial power of Pump4  (P4), (@Z=0)
      Exc.Pwr(:,6) = 0.0e-9;  % W, Initial power of Signal1(S1), (@Z=0)
      Exc.Pwr(:,7) = 0.0e-9;  % W, Initial power of Signal2(S2), (@Z=0)
    else
      if strcmp(SimMod.Scan_mode,'Pump_Power'),
        Exc.Pwr(:,1) = 1.0e-6; % W, Initial power of Signal1(S1), non-deplete
        Exc.Pwr(:,2) = 1.0e-6; % W, Initial power of Signal1(S2), non-deplete
      % Exc.Pwr(:,3) = linspace(1.0e-6, 10e-3,SimMod.Power_num); % W, Initial power of Pump1  (P1), degenerate pump
        Exc.Pwr(:,3) = 10.^(linspace(  5, 30,SimMod.Power_num)/10)*1e-3; 
        Exc.Pwr(:,4) = 0.0e-3; % W, Initial power of Pump3  (P3), (@Z=0)
        Exc.Pwr(:,5) = 0.0e-3; % W, Initial power of Pump4  (P4), (@Z=0)
        Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
        Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
      elseif strcmp(SimMod.Scan_mode,'Signal_Power'),
      % Exc.Pwr(:,1) = linspace( 10e-3,500e-3,SimMod.Power_num); % W, Initial power of Signal1(S1), non-deplete
        Exc.Pwr(:,1) = 10.^(linspace(-10, 40,SimMod.Power_num)/10)*1e-3; 
      % Exc.Pwr(:,2) = linspace( 10e-3,500e-3,SimMod.Power_num); % W, Initial power of Signal1(S2), non-deplete
        Exc.Pwr(:,2) = 10.^(linspace(-10, 40,SimMod.Power_num)/10)*1e-3; 
        Exc.Pwr(:,3) = 100e-3; % W, Initial power of Pump1  (P1), degenerate pump
        Exc.Pwr(:,4) = 0.0e-3; % W, Initial power of Pump3  (P3), (@Z=0)
        Exc.Pwr(:,5) = 0.0e-3; % W, Initial power of Pump4  (P4), (@Z=0)
        Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
        Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
      else
        Exc.Pwr(:,1) = 100e-3; % W, Initial power of Signal1(S1), non-deplete
        Exc.Pwr(:,2) = 100e-3; % W, Initial power of Signal1(S2), non-deplete
        Exc.Pwr(:,3) = 1.0e-6; % W, Initial power of Pump1  (P1), degenerate pump
        Exc.Pwr(:,4) = 0.0e-3; % W, Initial power of Pump3  (P3), (@Z=0)
        Exc.Pwr(:,5) = 0.0e-3; % W, Initial power of Pump4  (P4), (@Z=0)
        Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
        Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
      end; % end of if-else
    end; % end of if-else
  % end of case 1 of switch SimMod.Wave_mode
  case '3W2P'  
    if SimMod.Power_num == 1 
      Exc.Pwr(:,1) = 100e-3;  % W, Initial power of Pump1  (P1), non-degenerate
      Exc.Pwr(:,2) = 100e-3;  % W, Initial power of Pump2  (P2), non-degenerate
      Exc.Pwr(:,3) = 1.0e-6;  % W, Initial power of Signal0(S0), degenerate signal & Idler 
    % Exc.Pwr(:,3) = 10.^((7.5)/10)*1e-3; 
      Exc.Pwr(:,4) = 0.0e-6;  % W, Initial power of Pump3  (P3), (@Z=0)
      Exc.Pwr(:,5) = 0.0e-6;  % W, Initial power of Pump4  (P4), (@Z=0)
      Exc.Pwr(:,6) = 0.0e-9;  % W, Initial power of Signal1(S1), (@Z=0)
      Exc.Pwr(:,7) = 0.0e-9;  % W, Initial power of Signal2(S2), (@Z=0)
    else % SimMod.Power_num != 1
      switch SimMod.Scan_mode, 
        case 'Cen_Wavelength', % 
          Exc.Pwr(:,1) = 100e-3; % W, Initial power of Signal1(S1), non-deplete
          Exc.Pwr(:,2) = 100e-3; % W, Initial power of Signal1(S2), non-deplete
          Exc.Pwr(:,3) = 1.0e-6; % W, Initial power of Pump1  (P1), degenerate pump
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
        case 'PP_Separation',  % 
          Exc.Pwr(:,1) = 100e-3; % W, Initial power of Signal1(S1), non-deplete
          Exc.Pwr(:,2) = 100e-3; % W, Initial power of Signal1(S2), non-deplete
          Exc.Pwr(:,3) = 1.0e-6; % W, Initial power of Pump1  (P1), degenerate pump
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
        case 'Pump_Power',     % Scan pump power 
          Exc.Pwr(:,1) = linspace(  1e-3, 400e-3,SimMod.Power_num); % W, Initial power of Pump1(P1)
        % Exc.Pwr(:,1) = 10.^(linspace(  0, 23,SimMod.Power_num)/10)*1e-3; 
          Exc.Pwr(:,2) = linspace(  1e-3, 400e-3,SimMod.Power_num); % W, Initial power of Pump1(P2)
        % Exc.Pwr(:,2) = 10.^(linspace(  0, 23,SimMod.Power_num)/10)*1e-3; 
          Exc.Pwr(:,3) = 1.0e-3; % W, Initial power of degenerate signal&idler(SI)
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-6; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-6; % W, Initial power of Signal2(S2), (@Z=0)
        case 'Signal_Power',   % Scan signal power 
        % Exc.Pwr(:,1) = 100e-3; % W, Initial power of Pump1(P1), 
        % Exc.Pwr(:,2) = 100e-3; % W, Initial power of Pump2(P2), 
        % Exc.Pwr(:,3) = linspace(Exc.Pwr(1,1)*1e-3,Exc.Pwr(1,1),SimMod.Power_num); % W, 
        % Exc.Pwr(:,3) = 10.^(linspace(-40,+10,SimMod.Power_num)/10); % W, 
        
          Exc.Pwr(:,1) = 3.940375e-02; % W, Initial power of Pump1(P1), 
          Exc.Pwr(:,2) = 3.940375e-02; % W, Initial power of Pump2(P2),
          Exc.Pwr(:,3) = 10.^(linspace(10*log10(Exc.Pwr(1,1)+Exc.Pwr(1,2))-30,10*log10(Exc.Pwr(1,1)+Exc.Pwr(1,2)),SimMod.Power_num)/10); 
          
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
        case 'Relative_Phase', % 
          Exc.Pwr(:,1) = 100e-3; % W, Initial power of Signal1(S1), non-deplete
          Exc.Pwr(:,2) = 100e-3; % W, Initial power of Signal1(S2), non-deplete
          Exc.Pwr(:,3) = 1.0e-6; % W, Initial power of Pump1  (P1), degenerate pump
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
        case 'Propagation_Z',  % 
          Exc.Pwr(:,1) = 100e-3; % W, Initial power of Signal1(S1), non-deplete
          Exc.Pwr(:,2) = 100e-3; % W, Initial power of Signal1(S2), non-deplete
          Exc.Pwr(:,3) = 1.0e-6; % W, Initial power of Pump1  (P1), degenerate pump
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
        case 'Constellation',  % Scan signal power 
          
          Exc.Pwr(:,1) = 100e-3; % W, Initial power of Pump1(P1)
          Exc.Pwr(:,2) = 100e-3; % W, Initial power of Pump2(P2)
          Exc.Pwr(:,3) = 10.^(linspace(-0, 15,SimMod.Power_num)/10)*1e-3; 
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
        otherwise, 
          Exc.Pwr(:,1) = 100e-3; % W, Initial power of Signal1(S1), non-deplete
          Exc.Pwr(:,2) = 100e-3; % W, Initial power of Signal1(S2), non-deplete
          Exc.Pwr(:,3) = 1.0e-6; % W, Initial power of Pump1  (P1), degenerate pump
          Exc.Pwr(:,4) = 0.0e-6; % W, Initial power of Pump3  (P3), (@Z=0)
          Exc.Pwr(:,5) = 0.0e-6; % W, Initial power of Pump4  (P4), (@Z=0)
          Exc.Pwr(:,6) = 0.0e-9; % W, Initial power of Signal1(S1), (@Z=0)
          Exc.Pwr(:,7) = 0.0e-9; % W, Initial power of Signal2(S2), (@Z=0)
      end; % end of switch 
    end; % end of if-else 
  % end of case 1 of switch SimMod.Wave_mode
  case '4W2P'  
    if SimMod.Power_num == 1,
      %
    else
      %
    end; % end of if-else
  % end of case '4W2P' of switch SimMod.Wave_mode
  otherwise    
    %
end;

% Exc.Phi; % rad, initial phase of all the waves (@Z=0)
if strcmp(SimMod.Mod_type,'None'), 
  Exc.Phi = zeros(SimMod.Phase_num, SimMod.Wave_num); % Initial the initial phase of 7 waves
else  
  Exc.Phi = zeros(SimMod.Noise_num, SimMod.Wave_num); % Initial the initial phase of 7 waves
end;
if SimMod.Phase_num == 1,
  Exc.Phi(:,1) = 0; % rad, initial phase of Pump1  (P1), non-degenerate
  Exc.Phi(:,2) = 0; % rad, initial phase of Pump2  (P2), non-degenerate
  Exc.Phi(:,3) = 0; % rad, initial phase of Signal0(S0), degenerate signal & Idler 
  Exc.Phi(:,4) = 0; % rad, initial phase of Pump4, (P4), A(4)
  Exc.Phi(:,5) = 0; % rad, initial phase of Pump3, (P3), A(5)
  Exc.Phi(:,6) = 0; % rad, initi al phase of Signal1(S1), A(6)
  Exc.Phi(:,7) = 0; % rad, initial phase of Signal2(S2), A(7)
else % SimMod.Phase_num != 1
  switch SimMod.Scan_mode, 
    case 'Cen_Wavelength', 
      Exc.Phi(:,1) = 0.0; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,2) = 0.0; % rad, initial phase of Pump2  (P2), non-degenerate
      Exc.Phi(:,3) = linspace(-0.025*pi,+1.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
      Exc.Phi(:,4) = 0.0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0.0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0.0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0.0; % rad, initial phase of Signal2(S2), A(7)
    case 'PP_Separation',  
      Exc.Phi(:,1) = 0.0; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,2) = 0.0; % rad, initial phase of Pump2  (P2), non-degenerate
      Exc.Phi(:,3) = linspace(-0.025*pi,+1.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
      Exc.Phi(:,4) = 0.0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0.0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0.0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0.0; % rad, initial phase of Signal2(S2), A(7)
    case 'Pump_Power',     % Scan pump power 
      Exc.Phi(:,1) = 0.0; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,2) = 0.0; % rad, initial phase of Pump2  (P2), non-degenerate
      Exc.Phi(:,3) = linspace(-0.025*pi,+1.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
      Exc.Phi(:,4) = 0.0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0.0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0.0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0.0; % rad, initial phase of Signal2(S2), A(7)
    case 'Signal_Power',   % Scan signal power 
    % Exc.Phi(:,1) = 0.0; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,1) = linspace(-0.025*pi,+2.025*pi,SimMod.Phase_num); 
      Exc.Phi(:,2) = 0.0; % rad, initial phase of Pump2  (P2), non-degenerate
      Exc.Phi(:,3) = 0.0; ...linspace(-0.025*pi,+2.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
      Exc.Phi(:,4) = 0.0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0.0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0.0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0.0; % rad, initial phase of Signal2(S2), A(7)
    case 'Relative_Phase', % 
      Exc.Phi(:,1) = 0.0; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,2) = 0.0; % rad, initial phase of Pump2  (P2), non-degenerate
      Exc.Phi(:,3) = linspace(-0.025*pi,+1.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
      Exc.Phi(:,4) = 0.0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0.0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0.0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0.0; % rad, initial phase of Signal2(S2), A(7)
    case 'Propagation_Z',  % 
      Exc.Phi(:,1) = 0.0; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,2) = 0.0; % rad, initial phase of Pump2  (P2), non-degenerate
      Exc.Phi(:,3) = linspace(-0.025*pi,+1.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
      Exc.Phi(:,4) = 0.0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0.0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0.0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0.0; % rad, initial phase of Signal2(S2), A(7)
    case 'Constellation',  % Scan signal power 
      Exc.Phi(:,1) =+8.35*pi/8; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,2) = 0.00; % rad, initial phase of Pump2  (P2), non-degenerate
    % Exc.Phi(:,3) = linspace(-0.025*pi,+1.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
    % Exc.Phi(1:(SimMod.Phase_num-2)/2,      3) = 0.0*pi + Nis.pha_nis_wgn(1:(SimMod.Phase_num-2)/2,      1);
    % Exc.Phi(1+(SimMod.Phase_num-2)/2:end-2,3) = 1.0*pi + Nis.pha_nis_wgn(1+(SimMod.Phase_num-2)/2:end-2,1);
      Exc.Phi(1                         :  (SimMod.Noise_num-4)/4,3) = 0.0*pi + pi/4; 
      Exc.Phi(1+  (SimMod.Noise_num-4)/4:2*(SimMod.Noise_num-4)/4,3) = 0.5*pi + pi/4; 
      Exc.Phi(1+2*(SimMod.Noise_num-4)/4:3*(SimMod.Noise_num-4)/4,3) = 1.0*pi + pi/4; 
      Exc.Phi(1+3*(SimMod.Noise_num-4)/4:                   end-4,3) = 1.5*pi + pi/4; 
      Exc.Phi(end-3,3) =  0.0*pi + pi/4; 
      Exc.Phi(end-2,3) =  0.5*pi + pi/4; 
      Exc.Phi(end-1,3) =  1.0*pi + pi/4; 
      Exc.Phi(end  ,3) =  1.5*pi + pi/4; 
      Exc.Phi(:,4) = 0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0; % rad, initial phase of Signal2(S2), A(7)
    otherwise, 
      Exc.Phi(:,1) = 0.0; % rad, initial phase of Pump1  (P1), non-degenerate
      Exc.Phi(:,2) = 0.0; % rad, initial phase of Pump2  (P2), non-degenerate
      Exc.Phi(:,3) = linspace(-0.025*pi,+1.025*pi,SimMod.Phase_num); % rad, initial phase of Signal0(S0), degenerate signal & Idler 
      Exc.Phi(:,4) = 0.0; % rad, initial phase of Pump4, (P4), A(4)
      Exc.Phi(:,5) = 0.0; % rad, initial phase of Pump3, (P3), A(5)
      Exc.Phi(:,6) = 0.0; % rad, initial phase of Signal1(S1), A(6)
      Exc.Phi(:,7) = 0.0; % rad, initial phase of Signal2(S2), A(7)
  end; % end of switch 
end; % end of if-else SimMod.Power_num == 1

% Exc.Fld; % V, initial field of all the waves (@Z=0)
Exc.Fld      = zeros(1,                SimMod.Wave_num);

global Nis;
Nis = []; 

Nis.Num  = SimMod.Noise_num; 
Nis.wgn1 = randn(Nis.Num,1); 
Nis.wgn2 = randn(Nis.Num,1); 

fprintf(1,'[COMPLETED]\n');
%-------------------------------------------------------------------------------
  
%% Simulation results of all the waves 
fprintf(1,'Initializing memory and storage ... ');

global beta_matrix;
global d_kappa_Lin;

beta_matrix = zeros( 5,SimMod.WlSep_num); 
d_kappa_Lin = zeros(22,SimMod.WlSep_num); 

global theta;       % rad, relative phase of 13 NDFWM and 9 DFWM
global delta_theta; % rad, phase
global delta_beta;  % delta_Beta 01 ~ 22

theta       = zeros(22,1);
delta_theta = zeros(22,1);
delta_beta  = zeros(22,1);

global Var; 
global Raw; 
global Max; 
global PCS; 

f_Init_Mem();

fprintf(1,'[COMPLETED]\n');
%-------------------------------------------------------------------------------

%% Calculate the propagation 
fprintf(1,'Propagating in nonlinear medium ... ');

switch SimMod.Scan_mode,
  case 'Cen_Wavelength', 
    fprintf(1,'with varying Center_Wavelength ... \n');
    
    for j1 = 1 : SimMod.Wlofs_num,
      
      % Loop j1 for Cen_Wavelength
      
      ExcSet.Lambda_s = NLM.ZDW + ExcSet.Wvlen_ofs(j1);
      
      %-------------------------------------------------------------------------
      for j2 = 1 : SimMod.WlSep_num,

        % Loop j2 for PP_Separation

        % Calculate beta_vector and delta_beta(22,1) by Freq detuning
        [beta_matrix(:,j2),ExcSet.PP_Sp_Wvl(j2)] = f_delta_beta(ExcSet.Lambda_s,ExcSet.Sc_PPS_fq(j2)); 
        ExcSet.PP_Sp_Wvl(j2) = ExcSet.PP_Sp_Wvl(j2)*2e9; % Transfer SI to nm
        d_kappa_Lin(:,j2)    = delta_beta; 

        %-----------------------------------------------------------------------
        for j3 = 1 : SimMod.Power_num,

          % Loop j3 for initial Pump/Signal_Power
          
          %---------------------------------------------------------------------
          for j4 = 1 : SimMod.Phase_num,

            % Loop j4 for initial Relative_Phase
            
            %-------------------------------------------------------------------
            for j5 = 1 : 1,
              
              % Create fields
              f_Create_Fld(j1,j2,j3,j4,j5);
            
              % Calculate ODEs of 7-wave
              if SimMod.FLD_nPWR == 0,
                %% Calculating the propagation by solving ODEs of Power/Phase evolution
                if SimMod.PSA_nPIA == 0,
                  % Calculate the propagation in PIA mode
                  disp(['Calculate PIA @',num2str(ExcSet.Lambda_s*1e9,'%10.3f'),...
                        'nm with PP=',num2str(ExcSet.PP_Sp_Wvl(j2),'%10.3f'),...
                        'nm, P_{p}=',num2str(10*log10(Exc.Pwr(j3,1))+30,'%10.2f'), ...
                        'dBm, P_{s}=',num2str(10*log10(Exc.Pwr(j3,3))+30,'%10.2f'), ...
                        'dBm, using signal phase ', num2str(Exc.Phi(j4,3)/pi,'%10.3f'),...
                        'rad by 7w power']);         

                  options = odeset('RelTol', 1e-9, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                   Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                   Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                   Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                   Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                   options);  % Initial Conditions    
                  %-------------------------------------------------------------
                else % SimMod.PSA_nPIA == 1,
                  % Calculate the propagation in PSA mode
                  disp(['Calculate PSA @',num2str(ExcSet.Lambda_s*1e9,'%10.3f'),...
                        'nm with PP=',num2str(ExcSet.PP_Sp_Wvl(j2),'%10.3f'),...
                        'nm, P_p=',num2str(10*log10(Exc.Pwr(j3,1))+30,'%10.2f'), ...
                        'dBm, P_s=',num2str(10*log10(Exc.Pwr(j3,3))+30,'%10.2f'), ...
                        'dBm, using signal phase ', num2str(Exc.Phi(j4,3)/pi,'%10.3f'),...
                        'rad by 7w power']);

                  options = odeset('RelTol', 1e-9, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                   options);  % Initial Conditions
                  %-------------------------------------------------------------
                end; % end of if-else SimMod.PSA_nPIA

              % Pwr_Var(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
              % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
              
                Raw.P_P1_Lin(j2,j4) =      (P(end, 1));
                Raw.P_P2_Lin(j2,j4) =      (P(end, 2));
                Raw.P_SI_Lin(j2,j4) =      (P(end, 3));
                Raw.P_P3_Lin(j2,j4) =      (P(end, 4));
                Raw.P_P4_Lin(j2,j4) =      (P(end, 5));
                Raw.P_S1_Lin(j2,j4) =      (P(end, 6));
                Raw.P_S2_Lin(j2,j4) =      (P(end, 7));
              
              %---------------------------------------------------------------------------
              else % SimMod.FLD_nPWR == 1,
                %% Calculating the propagation by solving ODEs of Field evolution
                disp(['Calculating PSA @',num2str(ExcSet.Lambda_s*1e9,'%10.3f'), ...
                      'nm with PP=',num2str(ExcSet.PP_Sp_Wvl(j2),'%10.3f'), ...
                      'nm, P_p=',num2str(10*log10(Exc.Pwr(j3,1))+30,'%10.2f'), ...
                      'dBm, P_s=',num2str(10*log10(Exc.Pwr(j3,3))+30,'%10.2f'), ...
                      'dBm, using signal phase ', num2str(Exc.Phi(j4,3)/pi,'%10.3f'),...
                      'rad by 7w field']); 

                option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                [z,A]  = ode45(@dA_dz_7w, [0 ExcSet.PropZ_Len], ...
                               [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                               option); % Initial Conditions 
              %%
              % Fld_Var(:,j3) = [A(end,1) A(end,2) A(end,3) A(end,4) A(end,5) A(end,6) A(end,7)];
              % Pwr_Var(:,j3) = (abs(Fld_Var(:,j3))).^2; 
              % Phi_Var(:,j3) = [angle(A(end,1)) angle(A(end,2)) angle(A(end,3)) angle(A(end,4))...
              %                  angle(A(end,5)) angle(A(end,6)) angle(A(end,7))];
              %%
                Raw.F_P1_Lin(j2,j4) =       A(end, 1);
                Raw.F_P2_Lin(j2,j4) =       A(end, 2);
                Raw.F_SI_Lin(j2,j4) =       A(end, 3);
                Raw.F_P3_Lin(j2,j4) =       A(end, 4);
                Raw.F_P4_Lin(j2,j4) =       A(end, 5);
                Raw.F_S1_Lin(j2,j4) =       A(end, 6);
                Raw.F_S2_Lin(j2,j4) =       A(end, 7);

                Raw.P_P1_Lin(j2,j4) =  (abs(A(end, 1))).^2;
                Raw.P_P2_Lin(j2,j4) =  (abs(A(end, 2))).^2;
                Raw.P_SI_Lin(j2,j4) =  (abs(A(end, 3))).^2;
                Raw.P_P3_Lin(j2,j4) =  (abs(A(end, 4))).^2;
                Raw.P_P4_Lin(j2,j4) =  (abs(A(end, 5))).^2;
                Raw.P_S1_Lin(j2,j4) =  (abs(A(end, 6))).^2;
                Raw.P_S2_Lin(j2,j4) =  (abs(A(end, 7))).^2;

                Raw.Phase_P1(j2,j4) = angle(A(end,1));                   % rad,
                Raw.Phase_P2(j2,j4) = angle(A(end,2));                   % rad,
                Raw.Phase_SI(j2,j4) = angle(A(end,3));                   % rad,
                Raw.Phase_P3(j2,j4) = angle(A(end,4));                   % rad,
                Raw.Phase_P4(j2,j4) = angle(A(end,5));                   % rad,
                Raw.Phase_S1(j2,j4) = angle(A(end,6));                   % rad,
                Raw.Phase_S2(j2,j4) = angle(A(end,7));                   % rad,

                Raw.P_P1_Log(j2,j4) = 10*log10( Raw.P_P1_Lin(j2,j4).*1000 ); % dBm, 
                Raw.P_P2_Log(j2,j4) = 10*log10( Raw.P_P2_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_SI_Log(j2,j4) = 10*log10( Raw.P_SI_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_P4_Log(j2,j4) = 10*log10( Raw.P_P4_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_P3_Log(j2,j4) = 10*log10( Raw.P_P3_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_S1_Log(j2,j4) = 10*log10( Raw.P_S1_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_S2_Log(j2,j4) = 10*log10( Raw.P_S2_Lin(j2,j4).*1000 ); % dBm,

                Raw.G_SI_Lin(j2,j4) = Raw.P_SI_Lin(j2,j4)./Exc.Pwr(j3,3);    % Rt ,
                Raw.G_SI_Log(j2,j4) = 10*log10( Raw.G_SI_Lin(j2,j4) );       % dB ,

              end; % end of if-else SimMod.FLD_nPWR
              %-----------------------------------------------------------------
          
            end; % end of for j5 = 1 : SimMod.PropL_num,
            %-------------------------------------------------------------------

          end; % end of for j4 = 1 : SimMod.Phase_num,
          %---------------------------------------------------------------------

        end; % end of for j3 = 1 : SimMod.Power_num,
        %-----------------------------------------------------------------------
      %%
      % F_P1_Lin(j2,:) = Fld_Var(1,:);
      % F_P2_Lin(j2,:) = Fld_Var(2,:);
      % F_SI_Lin(j2,:) = Fld_Var(3,:);
      % F_P3_Lin(j2,:) = Fld_Var(4,:);
      % F_P4_Lin(j2,:) = Fld_Var(5,:);
      % F_S1_Lin(j2,:) = Fld_Var(6,:);
      % F_S2_Lin(j2,:) = Fld_Var(7,:);

      % P_P1_Log(j2,:) = 10*log10( Pwr_Var(1,:).*1000 ); % dBm, 
      % P_P2_Log(j2,:) = 10*log10( Pwr_Var(2,:).*1000 ); % dBm, 
      % P_SI_Log(j2,:) = 10*log10( Pwr_Var(3,:).*1000 ); % dBm,
      % P_P4_Log(j2,:) = 10*log10( Pwr_Var(4,:).*1000 ); % dBm, 
      % P_P3_Log(j2,:) = 10*log10( Pwr_Var(5,:).*1000 ); % dBm, 
      % P_S1_Log(j2,:) = 10*log10( Pwr_Var(6,:).*1000 ); % dBm,
      % P_S2_Log(j2,:) = 10*log10( Pwr_Var(7,:).*1000 ); % dBm,

      % Phase_SI(j2,:) = Phi_Var(3,:);                   % rad, 
      
      % G_SI_Lin(j2,:) = Pwr_Var(3,:)./Pwr( 1,3);
      % G_SI_Log(j2,:) = 10*log10( G_SI_Lin(j1,:) );     % dB ,

      % [G_SI_max(j2,1),G_SI_max(j2,2)] = max(G_SI_Log(j2,:));
      %  G_SI_max(j2,3) = Phi(G_SI_max(j2,2),3);
      % [G_SI_min(j2,1),G_SI_min(j2,2)] = min(G_SI_Log(j2,:));
      %  G_SI_min(j2,3) = Phi(G_SI_max(j2,2),3);
        
      % G_SI_zro_m(j2,j1) = G_SI_Log(j2,1);
      % G_SI_max_m(j2,j1) = G_SI_max(j2,1);
      % G_SI_min_m(j2,j1) = G_SI_min(j2,1);
        
      % G_SI_max_Ph_m(j2,j1) = G_SI_max(j2,3);
      % G_SI_min_Ph_m(j2,j1) = G_SI_min(j2,3);
      %%
       [Max.G_SI_max(j2, 1),Max.G_SI_max(j2, 2)] = max(Raw.G_SI_Log(j2,:));
        Max.G_SI_max(j2, 3) = Exc.Phi(Max.G_SI_max(j2, 2),3);
         
        Max.P_P1_max(j2, 1) = Raw.P_P1_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_P2_max(j2, 1) = Raw.P_P2_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_SI_max(j2, 1) = Raw.P_SI_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_P3_max(j2, 1) = Raw.P_P3_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_P4_max(j2, 1) = Raw.P_P4_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_S1_max(j2, 1) = Raw.P_S1_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_S2_max(j2, 1) = Raw.P_S2_Log(j2,Max.G_SI_max(j2, 2));
         
       [Max.G_SI_min(j2, 1),Max.G_SI_min(j2, 2)] = min(Raw.G_SI_Log(j2,:));
        Max.G_SI_min(j2, 3) = Exc.Phi(Max.G_SI_min(j2, 2),3);

        Max.G_SI_zro_m(j2,j1) = Raw.G_SI_Log(j2,1);
        Max.G_SI_max_m(j2,j1) = Max.G_SI_max(j2,1);
        Max.G_SI_min_m(j2,j1) = Max.G_SI_min(j2,1);
        Max.G_SI_PSER_m(j2,j1) = Max.G_SI_max(j2,1)-Max.G_SI_min(j2,1);
        Max.G_SI_PSGA_m(j2,j1) = Max.G_SI_max(j2,1)+Max.G_SI_min(j2,1);
        Max.G_SI_max_Ph_m(j2,j1) = Max.G_SI_max(j2,3);
        Max.G_SI_min_Ph_m(j2,j1) = Max.G_SI_min(j2,3);

      end; % end of for j2 = 1 : SimMod.WlSep_num,
      %-------------------------------------------------------------------------
      
      save([num2str(ExcSet.Lambda_s),'.mat']);
    end; % end of for j1 = 1 : SimMod.Wlofs_num,
    %---------------------------------------------------------------------------
  % end of case 'Cen_Wavelength' under 'Max/Min Gain' of switch 
  %-----------------------------------------------------------------------------
  case 'PP_Separation',  
    fprintf(1,'with varying PP_Separation ... \n');
    
    for j1 = 1 : SimMod.Wlofs_num,
      
      % Loop j1 for Cen_Wavelength
      
      ExcSet.Lambda_s = NLM.ZDW + ExcSet.Wvlen_ofs(j1);
      
      %-------------------------------------------------------------------------
      for j2 = 1 : SimMod.WlSep_num,

        % Loop j2 for PP_Separation

        % Calculate beta_vector and delta_beta(22,1) by Freq detuning
        [beta_matrix(:,j2),ExcSet.PP_Sp_Wvl(j2)] = f_delta_beta(ExcSet.Lambda_s,ExcSet.Sc_PPS_fq(j2)); 
        ExcSet.PP_Sp_Wvl(j2) = ExcSet.PP_Sp_Wvl(j2)*2e9; % Transfer SI to nm
        d_kappa_Lin(:,j2)    = delta_beta;  

        %-----------------------------------------------------------------------
        for j3 = 1 : SimMod.Power_num,
          
          % Loop j3 for initial Pump/Signal_Power
          
          %---------------------------------------------------------------------
          for j4 = 1 : SimMod.Phase_num,

            % Loop j4 for initial Relative_Phase
            
            %-------------------------------------------------------------------
            for j5 = 1 : 1,
              
              % Create fields
              f_Create_Fld(j1,j2,j3,j4,j5);
              
              % Calculate ODEs of 7-wave
              if SimMod.FLD_nPWR == 0,
                %% Calculating the propagation by solving ODEs of Power/Phase evolution
                if SimMod.PSA_nPIA == 0,
                  % Calculating in PIA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PSA case']);         

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dE_dz_7w, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                   Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                   Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                   Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                   Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                   options);  % Initial Conditions    

                else % SimMod.PSA_nPIA == 1,
                  % Calculating in PSA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PIA case']); 

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                   options);  % Initial Conditions

                end; % end of if-else SimMod.PSA_nPIA

                Pwr_Var(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
              % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
              %---------------------------------------------------------------------------
              else % SimMod.FLD_nPWR == 1,
                %% Calculating the propagation by solving ODEs of Field evolution
                disp(['Calculating PSA @',num2str(ExcSet.Lambda_s*1e9,'%10.3f'),...
                      'nm with PP=',num2str(ExcSet.PP_Sp_Wvl(j2),'%10.3f'), ...
                      'nm, P_p=',num2str(10*log10(Exc.Pwr(j3,1))+30,'%10.2f'), ...
                      'dBm, P_s=',num2str(10*log10(Exc.Pwr(j3,3))+30,'%10.2f'), ...
                      'dBm, using signal phase ', num2str(Exc.Phi(j4,3)/pi,'%10.3f'),...
                      'rad by 7w field']); 

                option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                [z,A]  = ode45(@dA_dz_7w, [0 ExcSet.PropZ_Len], ...
                               [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                               option); % Initial Conditions
                               
                Raw.F_P1_Lin(j2,j4) =       A(end, 1);
                Raw.F_P2_Lin(j2,j4) =       A(end, 2);
                Raw.F_SI_Lin(j2,j4) =       A(end, 3);
                Raw.F_P3_Lin(j2,j4) =       A(end, 4);
                Raw.F_P4_Lin(j2,j4) =       A(end, 5);
                Raw.F_S1_Lin(j2,j4) =       A(end, 6);
                Raw.F_S2_Lin(j2,j4) =       A(end, 7);

                Raw.P_P1_Lin(j2,j4) =  (abs(A(end, 1))).^2;
                Raw.P_P2_Lin(j2,j4) =  (abs(A(end, 2))).^2;
                Raw.P_SI_Lin(j2,j4) =  (abs(A(end, 3))).^2;
                Raw.P_P3_Lin(j2,j4) =  (abs(A(end, 4))).^2;
                Raw.P_P4_Lin(j2,j4) =  (abs(A(end, 5))).^2;
                Raw.P_S1_Lin(j2,j4) =  (abs(A(end, 6))).^2;
                Raw.P_S2_Lin(j2,j4) =  (abs(A(end, 7))).^2;

                Raw.Phase_P1(j2,j4) = angle(A(end,1));                   % rad,
                Raw.Phase_P2(j2,j4) = angle(A(end,2));                   % rad,
                Raw.Phase_SI(j2,j4) = angle(A(end,3));                   % rad,
                Raw.Phase_P3(j2,j4) = angle(A(end,4));                   % rad,
                Raw.Phase_P4(j2,j4) = angle(A(end,5));                   % rad,
                Raw.Phase_S1(j2,j4) = angle(A(end,6));                   % rad,
                Raw.Phase_S2(j2,j4) = angle(A(end,7));                   % rad,

                Raw.P_P1_Log(j2,j4) = 10*log10( Raw.P_P1_Lin(j2,j4).*1000 ); % dBm, 
                Raw.P_P2_Log(j2,j4) = 10*log10( Raw.P_P2_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_SI_Log(j2,j4) = 10*log10( Raw.P_SI_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_P4_Log(j2,j4) = 10*log10( Raw.P_P4_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_P3_Log(j2,j4) = 10*log10( Raw.P_P3_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_S1_Log(j2,j4) = 10*log10( Raw.P_S1_Lin(j2,j4).*1000 ); % dBm,
                Raw.P_S2_Log(j2,j4) = 10*log10( Raw.P_S2_Lin(j2,j4).*1000 ); % dBm,
                
                Raw.G_SI_Lin(j2,j4) = Raw.P_SI_Lin(j2,j4)./Exc.Pwr(j3,3);    % Rt ,
                Raw.G_SI_Log(j2,j4) = 10*log10( Raw.G_SI_Lin(j2,j4) );       % dB ,

              end; % end of if-else SimMod.FLD_nPWR
              %-----------------------------------------------------------------
              
            end; % end of for j5 = 1 : SimMod.PropL_num,
            %-------------------------------------------------------------------

          end; % end of for j4 = 1 : SimMod.Phase_num,
          %---------------------------------------------------------------------
          
        end; % end of for j3 = 1 : SimMod.Power_num,
        %-----------------------------------------------------------------------

       [Max.G_SI_max(j2, 1),Max.G_SI_max(j2, 2)] = max(Raw.G_SI_Log(j2,:));
        Max.G_SI_max(j2, 3) = Exc.Phi(Max.G_SI_max(j2, 2),3);
         
        Max.P_P1_max(j2, 1) = Raw.P_P1_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_P2_max(j2, 1) = Raw.P_P2_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_SI_max(j2, 1) = Raw.P_SI_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_P3_max(j2, 1) = Raw.P_P3_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_P4_max(j2, 1) = Raw.P_P4_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_S1_max(j2, 1) = Raw.P_S1_Log(j2,Max.G_SI_max(j2, 2));
        Max.P_S2_max(j2, 1) = Raw.P_S2_Log(j2,Max.G_SI_max(j2, 2));
         
       [Max.G_SI_min(j2, 1),Max.G_SI_min(j2, 2)] = min(Raw.G_SI_Log(j2,:));
        Max.G_SI_min(j2, 3) = Exc.Phi(Max.G_SI_min(j2, 2),3);

        Max.G_SI_zro_m(j2,j1) = Raw.G_SI_Log(j2,1);
        Max.G_SI_max_m(j2,j1) = Max.G_SI_max(j2,1);
        Max.G_SI_min_m(j2,j1) = Max.G_SI_min(j2,1);
        Max.G_SI_PSER_m(j2,j1) = Max.G_SI_max(j2,1)-Max.G_SI_min(j2,1);
        Max.G_SI_PSGA_m(j2,j1) = Max.G_SI_max(j2,1)+Max.G_SI_min(j2,1);        
        Max.G_SI_max_Ph_m(j2,j1) = Max.G_SI_max(j2,3);
        Max.G_SI_min_Ph_m(j2,j1) = Max.G_SI_min(j2,3);

      end; % end of for j2 = 1 : SimMod.WlSep_num,
      %-------------------------------------------------------------------------
        
    end; % end of for j1 = 1 : SimMod.Wlofs_num,
    %---------------------------------------------------------------------------0,NLM.Len
  % end of case PP_Separation of switch
  %-----------------------------------------------------------------------------
  case 'Pump_Power',     
    fprintf(1,'with varying Pump_Power ... \n');
    
    for j1 = 1 : SimMod.Wlofs_num,
      
      % Loop j1 for Cen_Wavelength
      
      ExcSet.Lambda_s = NLM.ZDW + ExcSet.Wvlen_ofs(j1);
      
      %-------------------------------------------------------------------------
      for j2 = 1 : SimMod.WlSep_num,
        
        % Loop j2 for PP_Separation
        
        % Calculate beta_vector and delta_beta(22,1) by Freq detuning
       [beta_matrix(:,j2),ExcSet.Sc_PPS_wl(j2)] = f_delta_beta(ExcSet.Lambda_s,ExcSet.Dl_PPS_fq(j2)); 
        ExcSet.Sc_PPS_wl(j2) = ExcSet.Sc_PPS_wl(j2)*2e9; 
        d_kappa_Lin(:,j2)    = delta_beta; 
        
        %-----------------------------------------------------------------------
        for j3 = 1 : SimMod.Power_num,
          
          % Loop j3 for initial Pump/Signal_Power
          
          %---------------------------------------------------------------------
          for j4 = 1 : SimMod.Phase_num,

            % Loop j4 for initial Relative_Phase
              
            switch SimMod.Scan_cond, 
              case 'Normal',    
                %---------------------------------------------------------------
                for j5 = 1 : 1,
                  
                  % Create field 
                  f_Create_Fld(j1,j2,j3,j4,j5);
                  
                  if SimMod.FLD_nPWR == 0, 
                    
                  else % SimMod.FLD_nPWR == 1,
                    %% Calculating the propagation by solving ODEs of Field evolution
                    disp(['Calculating PSA @',...
                          'Cen=', num2str(ExcSet.Lambda_s*1e9,         '%10.3f'),'nm ',...
                          'PPS=', num2str(ExcSet.Dl_PPS_wl(j2),        '%10.3f'),'nm ',...
                          'P_p=', num2str(10*log10(Exc.Pwr(j3,1))+30+3,'%10.2f'),'dBm ',...
                          'P_s=', num2str(10*log10(Exc.Pwr(j3,3))+30  ,'%10.2f'),'dBm ',...
                          'Ph_s=',num2str(Exc.Phi(j4,3)/pi,            '%10.3f'),'rad ',...
                          'Len=', num2str(NLM.Len,                     '%10.1f'),'m ',...
                          'by 7w field']); 

                    option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                    [z,A]  = ode45(@dA_dz_7w, [0 ExcSet.PropZ_Len], ...
                                   [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                    Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                                   option); % Initial Conditions

                    Raw.F_P1_Lin(j3,j4) =       A(end,1);
                    Raw.F_P2_Lin(j3,j4) =       A(end,2);
                    Raw.F_SI_Lin(j3,j4) =       A(end,3);
                    Raw.F_P3_Lin(j3,j4) =       A(end,4);
                    Raw.F_P4_Lin(j3,j4) =       A(end,5);
                    Raw.F_S1_Lin(j3,j4) =       A(end,6);
                    Raw.F_S2_Lin(j3,j4) =       A(end,7);

                    Raw.P_P1_Lin(j3,j4) =  (abs(A(end,1))).^2;
                    Raw.P_P2_Lin(j3,j4) =  (abs(A(end,2))).^2;
                    Raw.P_SI_Lin(j3,j4) =  (abs(A(end,3))).^2;
                    Raw.P_P3_Lin(j3,j4) =  (abs(A(end,4))).^2;
                    Raw.P_P4_Lin(j3,j4) =  (abs(A(end,5))).^2;
                    Raw.P_S1_Lin(j3,j4) =  (abs(A(end,6))).^2;
                    Raw.P_S2_Lin(j3,j4) =  (abs(A(end,7))).^2;

                    Raw.Phase_P1(j3,j4) = angle(A(end,1));                   % rad,
                    Raw.Phase_P2(j3,j4) = angle(A(end,2));                   % rad,
                    Raw.Phase_SI(j3,j4) = angle(A(end,3));                   % rad,
                    Raw.Phase_P3(j3,j4) = angle(A(end,4));                   % rad,
                    Raw.Phase_P4(j3,j4) = angle(A(end,5));                   % rad,
                    Raw.Phase_S1(j3,j4) = angle(A(end,6));                   % rad,
                    Raw.Phase_S2(j3,j4) = angle(A(end,7));                   % rad,

                    Raw.P_P1_Log(j3,j4) = 10*log10( Raw.P_P1_Lin(j3,j4).*1000 ); % dBm, 
                    Raw.P_P2_Log(j3,j4) = 10*log10( Raw.P_P2_Lin(j3,j4).*1000 ); % dBm,
                    Raw.P_SI_Log(j3,j4) = 10*log10( Raw.P_SI_Lin(j3,j4).*1000 ); % dBm,
                    Raw.P_P4_Log(j3,j4) = 10*log10( Raw.P_P4_Lin(j3,j4).*1000 ); % dBm,
                    Raw.P_P3_Log(j3,j4) = 10*log10( Raw.P_P3_Lin(j3,j4).*1000 ); % dBm,
                    Raw.P_S1_Log(j3,j4) = 10*log10( Raw.P_S1_Lin(j3,j4).*1000 ); % dBm,
                    Raw.P_S2_Log(j3,j4) = 10*log10( Raw.P_S2_Lin(j3,j4).*1000 ); % dBm,

                    Raw.G_SI_Lin(j3,j4) = Raw.P_SI_Lin(j3,j4)./Exc.Pwr(j3,3);    % Rt ,
                    Raw.G_SI_Log(j3,j4) = 10*log10( Raw.G_SI_Lin(j3,j4) );       % dB ,

                  end; % end of if-else SimMod.FLD_nPWR
                  %-------------------------------------------------------------
                
                end; % end of for j5 = 1 : SimMod.PropL_num,
                %---------------------------------------------------------------
                
              % end of case 'Normal' in switch SimMod.Scan_cond
              %-----------------------------------------------------------------
              case 'Propagate', 
                %---------------------------------------------------------------
                for j5 = 1 : 1,
                  
                  % Create field 
                  f_Create_Fld(j1,j2,j3,j4,j5);

                  % Calculate ODEs of 7-wave
                  if SimMod.FLD_nPWR == 0,
                    %% Calculating the propagation by solving ODEs of Power/Phase evolution
                    if SimMod.PSA_nPIA == 0,
                      % Calculating in PIA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PSA case']);         

                      options = odeset('RelTol', 1e-9, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                       Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                       Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                       Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                       Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                       options);  % Initial Conditions    

                    else % SimMod.PSA_nPIA == 1,
                      % Calculating in PSA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PIA case']); 

                      options = odeset('RelTol', 1e-9, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                       options);  % Initial Conditions

                    end; % end of if-else SimMod.PSA_nPIA

                    Var.Pwr(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
                  % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
                  %-------------------------------------------------------------
                  else % SimMod.FLD_nPWR == 1,
                    %% Calculating the propagation by solving ODEs of Field evolution
                    disp(['Calculating PSA @',num2str(ExcSet.Lambda_s*1e9,'%4.3f'),'nm ',...
                          'PPS=',num2str(ExcSet.Dl_PPS_wl(j2),'%4.3f'),'nm ',...
                          'P_p=',num2str(10*log10(Exc.Pwr(j3,1))+30,'%4.2f'),'dBm ',...
                          'P_s=',num2str(10*log10(Exc.Pwr(j3,3))+30,'%4.2f'),'dBm ',...
                          'Phi_sig=',num2str(Exc.Phi(j4,3)/pi,'%1.4f'),'rad ',...
                          num2str(ExcSet.PropZ_Len(end),'%10.3f'),'m, by 7w field']);

                    option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                    [z,A]  = ode45 (@dA_dz_7w,[0, ExcSet.PropZ_Len], ...
                                    [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                     Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                                    option); % Initial Conditions
                                  
                    %% Resorting outputs
                    
                    Raw.F_P1_Lin(j4,:) =       A(end, 1); % V, A1 field linear
                    Raw.F_P2_Lin(j4,:) =       A(end, 2); % V, A2 field linear
                    Raw.F_SI_Lin(j4,:) =       A(end, 3); % V, A0 field linear
                    Raw.F_P3_Lin(j4,:) =       A(end, 4); 
                    Raw.F_P4_Lin(j4,:) =       A(end, 5);
                    Raw.F_S1_Lin(j4,:) =       A(end, 6);
                    Raw.F_S2_Lin(j4,:) =       A(end, 7);

                    Raw.P_P1_Lin(j4,:) =  (abs(A(end, 1))).^2; % W, A1 power lin
                    Raw.P_P2_Lin(j4,:) =  (abs(A(end, 2))).^2; % W, A2 power lin
                    Raw.P_SI_Lin(j4,:) =  (abs(A(end, 3))).^2; 
                    Raw.P_P3_Lin(j4,:) =  (abs(A(end, 4))).^2;
                    Raw.P_P4_Lin(j4,:) =  (abs(A(end, 5))).^2;
                    Raw.P_S1_Lin(j4,:) =  (abs(A(end, 6))).^2;
                    Raw.P_S2_Lin(j4,:) =  (abs(A(end, 7))).^2;

                    Raw.Phase_P1(j4,:) = angle(A(end,1));  % rad, A1 phase lin
                    Raw.Phase_P2(j4,:) = angle(A(end,2));  % rad, A2 phase lin
                    Raw.Phase_SI(j4,:) = angle(A(end,3));  % rad, A0 phase lin
                    Raw.Phase_P3(j4,:) = angle(A(end,4));  % rad,
                    Raw.Phase_P4(j4,:) = angle(A(end,5));  % rad,
                    Raw.Phase_S1(j4,:) = angle(A(end,6));  % rad,
                    Raw.Phase_S2(j4,:) = angle(A(end,7));  % rad,

                    Raw.P_P1_Log(j4,:) = 10*log10( Raw.P_P1_Lin(j4,:).*1000 ); % dBm, A1 power log
                    Raw.P_P2_Log(j4,:) = 10*log10( Raw.P_P2_Lin(j4,:).*1000 ); % dBm, 
                    Raw.P_SI_Log(j4,:) = 10*log10( Raw.P_SI_Lin(j4,:).*1000 ); % dBm, A0 power log
                    Raw.P_P3_Log(j4,:) = 10*log10( Raw.P_P3_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_P4_Log(j4,:) = 10*log10( Raw.P_P4_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_S1_Log(j4,:) = 10*log10( Raw.P_S1_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_S2_Log(j4,:) = 10*log10( Raw.P_S2_Lin(j4,:).*1000 ); % dBm,

                    Raw.G_SI_Lin(j4,:) = Raw.P_SI_Lin(j4,:)./Exc.Pwr(j3,3);    % Rt ,
                    Raw.G_SI_Log(j4,:) = 10*log10( Raw.G_SI_Lin(j4,:) );       % dB ,

                  end; % end of if-else SimMod.FLD_nPWR
                  %-------------------------------------------------------------
                end; % end of for j5 = 1 : 1 in 'Propagate',
                %---------------------------------------------------------------
                
              % end of case 'Propagate' in switch SimMod.Scan_cond
              %-----------------------------------------------------------------
              case 'Len_Z_var', 
                %---------------------------------------------------------------
                for j5 = 2 : SimMod.PropZ_num,

                  f_Create_Fld(j1,j2,j3,j4,j5);

                  % calculate DE
                  if SimMod.FLD_nPWR == 0,
                    %% Calculating the propagation by solving ODEs of Power/Phase evolution 
                    if SimMod.PSA_nPIA == 0,
                      % Calculating in PIA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PSA case']);         

                      options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                       Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                       Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                       Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                       Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                       options);  % Initial Conditions    

                    else % SimMod.PSA_nPIA == 1,
                      % Calculating in PSA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PIA case']); 

                      options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                       options);  % Initial Conditions

                    end; % end of if-else SimMod.PSA_nPIA

                    Var.Pwr(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
                  % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
                  %---------------------------------------------------------------------------
                  else % SimMod.FLD_nPWR == 1,
                    %% Calculating the propagation by solving ODEs of field evolution 
                    disp(['Calculating PSA @',...
                          'Cen=', num2str(ExcSet.Lambda_s*1e9,         '%10.3f'),'nm ',...
                          'PPS=', num2str(ExcSet.Dl_PPS_wl(j2),        '%10.3f'),'nm ',...
                          'P_p=', num2str(10*log10(Exc.Pwr(j3,1)+Exc.Pwr(j3,2))+30,'%10.2f'),'dBm ',...
                          'P_s=', num2str(10*log10(Exc.Pwr(j3,3))+30  ,'%10.2f'),'dBm ',...
                          'Ph_s=',num2str(Exc.Phi(j4,3)/pi,            '%10.3f'),'rad ',...
                          'Len=', num2str(NLM.Len,                     '%10.1f'),'m ',...
                          'by 7w field']); 

                    option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                    [z,A]  = ode45 (@dA_dz_7w, [0 ExcSet.PropZ_Len(j5)], ...
                                    [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                     Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                                    option); % Initial Conditions

                    Raw.F_P1_Lin(j4,j5) =       A(end, 1);
                    Raw.F_P2_Lin(j4,j5) =       A(end, 2);
                    Raw.F_SI_Lin(j4,j5) =       A(end, 3);
                    Raw.F_P3_Lin(j4,j5) =       A(end, 4);
                    Raw.F_P4_Lin(j4,j5) =       A(end, 5);
                    Raw.F_S1_Lin(j4,j5) =       A(end, 6);
                    Raw.F_S2_Lin(j4,j5) =       A(end, 7);

                    Raw.P_P1_Lin(j4,j5) =  (abs(A(end, 1))).^2;
                    Raw.P_P2_Lin(j4,j5) =  (abs(A(end, 2))).^2;
                    Raw.P_SI_Lin(j4,j5) =  (abs(A(end, 3))).^2;
                    Raw.P_P3_Lin(j4,j5) =  (abs(A(end, 4))).^2;
                    Raw.P_P4_Lin(j4,j5) =  (abs(A(end, 5))).^2;
                    Raw.P_S1_Lin(j4,j5) =  (abs(A(end, 6))).^2;
                    Raw.P_S2_Lin(j4,j5) =  (abs(A(end, 7))).^2;

                    Raw.Phase_P1(j4,j5) = angle(A(end,1));                   % rad,
                    Raw.Phase_P2(j4,j5) = angle(A(end,2));                   % rad,
                    Raw.Phase_SI(j4,j5) = angle(A(end,3));                   % rad,
                    Raw.Phase_P3(j4,j5) = angle(A(end,4));                   % rad,
                    Raw.Phase_P4(j4,j5) = angle(A(end,5));                   % rad,
                    Raw.Phase_S1(j4,j5) = angle(A(end,6));                   % rad,
                    Raw.Phase_S2(j4,j5) = angle(A(end,7));                   % rad,

                    Raw.P_P1_Log(j4,j5) = 10*log10( Raw.P_P1_Lin(j4,j5).*1000 ); % dBm, 
                    Raw.P_P2_Log(j4,j5) = 10*log10( Raw.P_P2_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_SI_Log(j4,j5) = 10*log10( Raw.P_SI_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_P3_Log(j4,j5) = 10*log10( Raw.P_P3_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_P4_Log(j4,j5) = 10*log10( Raw.P_P4_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_S1_Log(j4,j5) = 10*log10( Raw.P_S1_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_S2_Log(j4,j5) = 10*log10( Raw.P_S2_Lin(j4,j5).*1000 ); % dBm,

                    Raw.G_SI_Lin(j4,j5) = Raw.P_SI_Lin(j4,j5)./Exc.Pwr(j3,3);    % Rt ,
                    Raw.G_SI_Log(j4,j5) = 10*log10( Raw.G_SI_Lin(j4,j5) );       % dB ,

                  end; % end of if-else SimMod.FLD_nPWR
                  %-------------------------------------------------------------
                end; % end of for j5 = 1 : SimMod.PropZ_num,
                %---------------------------------------------------------------
                
                % Data processing to be completed
                
              % end of case 'Len_Z_var' in switch SimMod.Scan_cond
              %-----------------------------------------------------------------
            end; % end of switch SimMod.Scan_cond
            %-------------------------------------------------------------------
            
          end; % end of for j4 = 1 : SimMod.Phase_num,
          %---------------------------------------------------------------------
          
          switch SimMod.Scan_cond, 
            case 'Normal',   % 
             [Max.G_SI_max(j3, 1),Max.G_SI_max(j3, 2)] = max(Raw.G_SI_Log(j3,:));
              Max.G_SI_max(j3, 3) = Exc.Phi(Max.G_SI_max(j3, 2),3);

             [Max.G_SI_min(j3, 1),Max.G_SI_min(j3, 2)] = min(Raw.G_SI_Log(j3,:));
              Max.G_SI_min(j3, 3) = Exc.Phi(Max.G_SI_min(j3, 2),3);

              Max.G_SI_zero_m(j3,j1)   = Raw.G_SI_Log(j3,1);
              Max.G_SI_max_m (j3,j1)   = Max.G_SI_max(j3,1);
              Max.G_SI_min_m (j3,j1)   = Max.G_SI_min(j3,1);
              Max.G_SI_PSER_m(j3,j1)   = Max.G_SI_max(j3,1)-Max.G_SI_min(j3,1);
              Max.G_SI_PSGA_m(j3,j1)   = Max.G_SI_max(j3,1)+Max.G_SI_min(j3,1);
              Max.G_SI_max_Ph_m(j3,j1) = Max.G_SI_max(j3,3);
              Max.G_SI_min_Ph_m(j3,j1) = Max.G_SI_min(j3,3);
            %-------------------------------------------------------------------
            case 'Propagate',% 
              % Peak values
             [Max.G_SI_max(j3, 1),Max.G_SI_max(j3, 2)] = max(Raw.G_SI_Log(:,end)); % Gmax, ind
              Max.G_SI_max(j3, 3) = Exc.Phi(Max.G_SI_max(j3, 2),3);                % Phase

             [Max.G_SI_min(j3, 1),Max.G_SI_min(j3, 2)] = min(Raw.G_SI_Log(:,end)); % 
              Max.G_SI_min(j3, 3) = Exc.Phi(Max.G_SI_min(j3, 2),3);

              Max.G_SI_zro_m   (j3,j2) = Raw.G_SI_Log(1,end);
              Max.G_SI_max_m   (j3,j2) = Max.G_SI_max(j3,1);
              Max.G_SI_min_m   (j3,j2) = Max.G_SI_min(j3,1);
              Max.G_SI_PSER_m  (j3,j2) = Max.G_SI_max(j3,1)-Max.G_SI_min(j3,1);
              Max.G_SI_PSGA_m  (j3,j2) = Max.G_SI_max(j3,1)+Max.G_SI_min(j3,1);
              Max.G_SI_max_Ph_m(j3,j2) = Max.G_SI_max(j3,3);
              Max.G_SI_min_Ph_m(j3,j2) = Max.G_SI_min(j3,3);

              % Evolution for Gmax
              Max.F_P1_Lin_max_evo (j3, :) = Raw.F_P1_Lin(Max.G_SI_max(j3, 2),:);
              Max.F_P2_Lin_max_evo (j3, :) = Raw.F_P2_Lin(Max.G_SI_max(j3, 2),:);
              Max.F_SI_Lin_max_evo (j3, :) = Raw.F_SI_Lin(Max.G_SI_max(j3, 2),:);
              Max.F_P3_Lin_max_evo (j3, :) = Raw.F_P3_Lin(Max.G_SI_max(j3, 2),:);
              Max.F_P4_Lin_max_evo (j3, :) = Raw.F_P4_Lin(Max.G_SI_max(j3, 2),:);
              Max.F_S1_Lin_max_evo (j3, :) = Raw.F_S1_Lin(Max.G_SI_max(j3, 2),:);
              Max.F_S2_Lin_max_evo (j3, :) = Raw.F_S2_Lin(Max.G_SI_max(j3, 2),:);

              Max.P_P1_Lin_max_evo (j3, :) = Raw.P_P1_Lin(Max.G_SI_max(j3, 2),:);
              Max.P_P2_Lin_max_evo (j3, :) = Raw.P_P2_Lin(Max.G_SI_max(j3, 2),:);
              Max.P_SI_Lin_max_evo (j3, :) = Raw.P_SI_Lin(Max.G_SI_max(j3, 2),:);
              Max.P_P3_Lin_max_evo (j3, :) = Raw.P_P3_Lin(Max.G_SI_max(j3, 2),:);
              Max.P_P4_Lin_max_evo (j3, :) = Raw.P_P4_Lin(Max.G_SI_max(j3, 2),:);
              Max.P_S1_Lin_max_evo (j3, :) = Raw.P_S1_Lin(Max.G_SI_max(j3, 2),:);
              Max.P_S2_Lin_max_evo (j3, :) = Raw.P_S2_Lin(Max.G_SI_max(j3, 2),:);

              Max.Phase_P1_max_evo (j3, :) = Raw.Phase_P1(Max.G_SI_max(j3, 2),:);
              Max.Phase_P2_max_evo (j3, :) = Raw.Phase_P2(Max.G_SI_max(j3, 2),:);
              Max.Phase_SI_max_evo (j3, :) = Raw.Phase_SI(Max.G_SI_max(j3, 2),:);
              Max.Phase_P3_max_evo (j3, :) = Raw.Phase_P3(Max.G_SI_max(j3, 2),:);
              Max.Phase_P4_max_evo (j3, :) = Raw.Phase_P4(Max.G_SI_max(j3, 2),:);
              Max.Phase_S1_max_evo (j3, :) = Raw.Phase_S1(Max.G_SI_max(j3, 2),:);
              Max.Phase_S2_max_evo (j3, :) = Raw.Phase_S2(Max.G_SI_max(j3, 2),:);

              Max.P_P1_Log_max_evo (j3, :) = Raw.P_P1_Log(Max.G_SI_max(j3, 2),:);
              Max.P_P2_Log_max_evo (j3, :) = Raw.P_P2_Log(Max.G_SI_max(j3, 2),:);
              Max.P_SI_Log_max_evo (j3, :) = Raw.P_SI_Log(Max.G_SI_max(j3, 2),:);
              Max.P_P3_Log_max_evo (j3, :) = Raw.P_P3_Log(Max.G_SI_max(j3, 2),:);
              Max.P_P4_Log_max_evo (j3, :) = Raw.P_P4_Log(Max.G_SI_max(j3, 2),:);
              Max.P_S1_Log_max_evo (j3, :) = Raw.P_S1_Log(Max.G_SI_max(j3, 2),:);
              Max.P_S2_Log_max_evo (j3, :) = Raw.P_S2_Log(Max.G_SI_max(j3, 2),:);

              % Evolution for Gmin
              Max.F_P1_Lin_min_evo (j3, :) = Raw.F_P1_Lin(Max.G_SI_max(j3, 2),:);
              Max.F_P2_Lin_min_evo (j3, :) = Raw.F_P2_Lin(Max.G_SI_min(j3, 2),:);
              Max.F_SI_Lin_min_evo (j3, :) = Raw.F_SI_Lin(Max.G_SI_min(j3, 2),:);
              Max.F_P3_Lin_min_evo (j3, :) = Raw.F_P3_Lin(Max.G_SI_min(j3, 2),:);
              Max.F_P4_Lin_min_evo (j3, :) = Raw.F_P4_Lin(Max.G_SI_min(j3, 2),:);
              Max.F_S1_Lin_min_evo (j3, :) = Raw.F_S1_Lin(Max.G_SI_min(j3, 2),:);
              Max.F_S2_Lin_min_evo (j3, :) = Raw.F_S2_Lin(Max.G_SI_min(j3, 2),:);

              Max.P_P1_Lin_min_evo (j3, :) = Raw.P_P1_Lin(Max.G_SI_min(j3, 2),:);
              Max.P_P2_Lin_min_evo (j3, :) = Raw.P_P2_Lin(Max.G_SI_min(j3, 2),:);
              Max.P_SI_Lin_min_evo (j3, :) = Raw.P_SI_Lin(Max.G_SI_min(j3, 2),:);
              Max.P_P3_Lin_min_evo (j3, :) = Raw.P_P3_Lin(Max.G_SI_min(j3, 2),:);
              Max.P_P4_Lin_min_evo (j3, :) = Raw.P_P4_Lin(Max.G_SI_min(j3, 2),:);
              Max.P_S1_Lin_min_evo (j3, :) = Raw.P_S1_Lin(Max.G_SI_min(j3, 2),:);
              Max.P_S2_Lin_min_evo (j3, :) = Raw.P_S2_Lin(Max.G_SI_min(j3, 2),:);

              Max.Phase_P1_min_evo (j3, :) = Raw.Phase_P1(Max.G_SI_min(j3, 2),:);
              Max.Phase_P2_min_evo (j3, :) = Raw.Phase_P2(Max.G_SI_min(j3, 2),:);
              Max.Phase_SI_min_evo (j3, :) = Raw.Phase_SI(Max.G_SI_min(j3, 2),:);
              Max.Phase_P3_min_evo (j3, :) = Raw.Phase_P3(Max.G_SI_min(j3, 2),:);
              Max.Phase_P4_min_evo (j3, :) = Raw.Phase_P4(Max.G_SI_min(j3, 2),:);
              Max.Phase_S1_min_evo (j3, :) = Raw.Phase_S1(Max.G_SI_min(j3, 2),:);
              Max.Phase_S2_min_evo (j3, :) = Raw.Phase_S2(Max.G_SI_min(j3, 2),:);

              Max.P_P1_Log_min_evo (j3, :) = Raw.P_P1_Log(Max.G_SI_min(j3, 2),:);
              Max.P_P2_Log_min_evo (j3, :) = Raw.P_P2_Log(Max.G_SI_min(j3, 2),:);
              Max.P_SI_Log_min_evo (j3, :) = Raw.P_SI_Log(Max.G_SI_min(j3, 2),:);
              Max.P_P3_Log_min_evo (j3, :) = Raw.P_P3_Log(Max.G_SI_min(j3, 2),:);
              Max.P_P4_Log_min_evo (j3, :) = Raw.P_P4_Log(Max.G_SI_min(j3, 2),:);
              Max.P_S1_Log_min_evo (j3, :) = Raw.P_S1_Log(Max.G_SI_min(j3, 2),:);
              Max.P_S2_Log_min_evo (j3, :) = Raw.P_S2_Log(Max.G_SI_min(j3, 2),:);
            %-------------------------------------------------------------------
            case 'Len_Z_var',% 
            %-------------------------------------------------------------------
          end; % end of switch 
          
        end; % end of for j3 = 1 : SimMod.Power_num,
        %-----------------------------------------------------------------------

      end; % end of for j2 = 1 : SimMod.WlSep_num,
      %-------------------------------------------------------------------------
        
    end; % end of for j1 = 1 : SimMod.Wlofs_num,
    %---------------------------------------------------------------------------
  % end of case 'Pump_Power' in switch SimMod.Scan_mode
  %-----------------------------------------------------------------------------
  case 'Signal_Power',   
    fprintf(1,'with varying Signal_Power ... \n');
    
    for j1 = 1 : SimMod.Wlofs_num,
      
      % Loop j1 for Cen_Wavelength
      
      ExcSet.Lambda_s = NLM.ZDW + ExcSet.Wvlen_ofs(j1);
      
      %-------------------------------------------------------------------------
      for j2 = 1 : SimMod.WlSep_num,

        % Loop j2 for PP_Separation

        % Calculate beta_vector and delta_beta(22,1) by Freq detuning
       [beta_matrix(:,j2),ExcSet.Sc_PPS_wl(j2)] = f_delta_beta(ExcSet.Lambda_s,ExcSet.Dl_PPS_fq(j2)); 
        ExcSet.Sc_PPS_wl(j2) = ExcSet.Sc_PPS_wl(j2)*2e9; % Transfer SI to nm
        d_kappa_Lin(:,j2)    = delta_beta;  

        %-----------------------------------------------------------------------
        for j3 = 1 : SimMod.Power_num,
          
          % Loop j3 for initial Pump/Signal_Power
          
          for j4 = 1 : SimMod.Phase_num,

            % Loop j4 for initial Relative_Phase
            
            %-------------------------------------------------------------------
            for j5 = 1 : 1,
              
              % Create fields
              f_Create_Fld(j1,j2,j3,j4,j5);
              
              % Calculate ODEs of 7-wave
              if SimMod.FLD_nPWR == 0,
                %% Calculating the propagation by solving ODEs of Power/Phase evolution
                if SimMod.PSA_nPIA == 0,
                  % Calculating in PIA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PSA case']);         

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                   Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                   Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                   Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                   Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                   options);  % Initial Conditions    

                else % SimMod.PSA_nPIA == 1,
                  % Calculating in PSA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PIA case']); 

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                   options);  % Initial Conditions

                end; % end of if-else SimMod.PSA_nPIA

                Pwr_Var(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
              % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
              %---------------------------------------------------------------------------
              else % SimMod.FLD_nPWR == 1,
                %% Calculating the propagation by solving ODEs of Field evolution
                disp(['Calculating PSA @',...
                      'Cen=', num2str(ExcSet.Lambda_s*1e9,         '%10.3f'),'nm ',...
                      'PPS=', num2str(ExcSet.Dl_PPS_wl(j2),        '%10.3f'),'nm ',...
                      'P_p=', num2str(10*log10(Exc.Pwr(j3,1))+30+3,'%10.2f'),'dBm ',...
                      'P_s=', num2str(10*log10(Exc.Pwr(j3,3))+30  ,'%10.2f'),'dBm ',...
                      'Ph_s=',num2str(Exc.Phi(j4,1)/pi,            '%10.3f'),'rad ',...
                      'Len=', num2str(NLM.Len,                     '%10.1f'),'m ',...
                      'by 7w field']); 

                option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                [z,A]  = ode45(@dA_dz_7w, [0 ExcSet.PropZ_Len], ...
                               [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                               option); % Initial Conditions
                               
                Raw.F_P1_Lin(j3,j4) =       A(end, 1);
                Raw.F_P2_Lin(j3,j4) =       A(end, 2);
                Raw.F_SI_Lin(j3,j4) =       A(end, 3);
                Raw.F_P3_Lin(j3,j4) =       A(end, 4);
                Raw.F_P4_Lin(j3,j4) =       A(end, 5);
                Raw.F_S1_Lin(j3,j4) =       A(end, 6);
                Raw.F_S2_Lin(j3,j4) =       A(end, 7);

                Raw.P_P1_Lin(j3,j4) =  (abs(A(end, 1))).^2;
                Raw.P_P2_Lin(j3,j4) =  (abs(A(end, 2))).^2;
                Raw.P_SI_Lin(j3,j4) =  (abs(A(end, 3))).^2;
                Raw.P_P3_Lin(j3,j4) =  (abs(A(end, 4))).^2;
                Raw.P_P4_Lin(j3,j4) =  (abs(A(end, 5))).^2;
                Raw.P_S1_Lin(j3,j4) =  (abs(A(end, 6))).^2;
                Raw.P_S2_Lin(j3,j4) =  (abs(A(end, 7))).^2;

                Raw.Phase_P1(j3,j4) = angle(A(end,1));                   % rad,
                Raw.Phase_P2(j3,j4) = angle(A(end,2));                   % rad,
                Raw.Phase_SI(j3,j4) = angle(A(end,3));                   % rad,
                Raw.Phase_P3(j3,j4) = angle(A(end,4));                   % rad,
                Raw.Phase_P4(j3,j4) = angle(A(end,5));                   % rad,
                Raw.Phase_S1(j3,j4) = angle(A(end,6));                   % rad,
                Raw.Phase_S2(j3,j4) = angle(A(end,7));                   % rad,

                Raw.P_P1_Log(j3,j4) = 10*log10( Raw.P_P1_Lin(j3,j4).*1000 ); % dBm, 
                Raw.P_P2_Log(j3,j4) = 10*log10( Raw.P_P2_Lin(j3,j4).*1000 ); % dBm,
                Raw.P_SI_Log(j3,j4) = 10*log10( Raw.P_SI_Lin(j3,j4).*1000 ); % dBm,
                Raw.P_P4_Log(j3,j4) = 10*log10( Raw.P_P4_Lin(j3,j4).*1000 ); % dBm,
                Raw.P_P3_Log(j3,j4) = 10*log10( Raw.P_P3_Lin(j3,j4).*1000 ); % dBm,
                Raw.P_S1_Log(j3,j4) = 10*log10( Raw.P_S1_Lin(j3,j4).*1000 ); % dBm,
                Raw.P_S2_Log(j3,j4) = 10*log10( Raw.P_S2_Lin(j3,j4).*1000 ); % dBm,
                
                Raw.G_SI_Lin(j3,j4) = Raw.P_SI_Lin(j3,j4)./Exc.Pwr(j3,3);    % Rt ,
                Raw.G_SI_Log(j3,j4) = 10*log10( Raw.G_SI_Lin(j3,j4) );       % dB ,

              end; % end of if-else SimMod.FLD_nPWR
              %-----------------------------------------------------------------
            end; % end of for j5 = 1 : SimMod.PropZ_num,
            %-------------------------------------------------------------------

          end; % end of for j4 = 1 : SimMod.Phase_num,
          %---------------------------------------------------------------------
          %% Ext: Maximum gain and ... 
         [Max.G_SI_max(j3, 1),Max.G_SI_max(j3, 2)] = max(Raw.G_SI_Log(j3,:));
          Max.G_SI_max(j3, 3) = Exc.Phi(Max.G_SI_max(j3, 2),3);
          
          Max.P_P1_max(j3, 1) = Raw.P_P1_Log(j3,Max.G_SI_max(j3, 2));
          Max.P_P2_max(j3, 1) = Raw.P_P2_Log(j3,Max.G_SI_max(j3, 2));
          Max.P_SI_max(j3, 1) = Raw.P_SI_Log(j3,Max.G_SI_max(j3, 2));
          Max.P_P3_max(j3, 1) = Raw.P_P3_Log(j3,Max.G_SI_max(j3, 2));
          Max.P_P4_max(j3, 1) = Raw.P_P4_Log(j3,Max.G_SI_max(j3, 2));
          Max.P_S1_max(j3, 1) = Raw.P_S1_Log(j3,Max.G_SI_max(j3, 2));
          Max.P_S2_max(j3, 1) = Raw.P_S2_Log(j3,Max.G_SI_max(j3, 2));

         [Max.G_SI_min(j3, 1),Max.G_SI_min(j3, 2)] = min(Raw.G_SI_Log(j3,:));
          Max.G_SI_min(j3, 3) = Exc.Phi(Max.G_SI_min(j3, 2),3);
          
          Max.P_P1_min(j3, 1) = Raw.P_P1_Log(j3,Max.G_SI_min(j3, 2));
          Max.P_P2_min(j3, 1) = Raw.P_P2_Log(j3,Max.G_SI_min(j3, 2));
          Max.P_SI_min(j3, 1) = Raw.P_SI_Log(j3,Max.G_SI_min(j3, 2));
          Max.P_P3_min(j3, 1) = Raw.P_P3_Log(j3,Max.G_SI_min(j3, 2));
          Max.P_P4_min(j3, 1) = Raw.P_P4_Log(j3,Max.G_SI_min(j3, 2));
          Max.P_S1_min(j3, 1) = Raw.P_S1_Log(j3,Max.G_SI_min(j3, 2));
          Max.P_S2_min(j3, 1) = Raw.P_S2_Log(j3,Max.G_SI_min(j3, 2));

          Max.G_SI_zro_m(j3,j1) = Raw.G_SI_Log(j3,1);
          Max.G_SI_max_m(j3,j1) = Max.G_SI_max(j3,1);
          Max.G_SI_min_m(j3,j1) = Max.G_SI_min(j3,1);
          Max.G_SI_PSER_m(j3,j1) = Max.G_SI_max(j3,1)-Max.G_SI_min(j3,1);
          Max.G_SI_PSGA_m(j3,j1) = Max.G_SI_max(j3,1)+Max.G_SI_min(j3,1);
          Max.G_SI_max_Ph_m(j3,j1) = Max.G_SI_max(j3,3);
          Max.G_SI_min_Ph_m(j3,j1) = Max.G_SI_min(j3,3);
          
        end; % end of for j3 = 1 : SimMod.Power_num,
        %-----------------------------------------------------------------------

      end; % end of for j2 = 1 : SimMod.WlSep_num,
      %-------------------------------------------------------------------------
      
    end; % end of for j1 = 1 : SimMod.Wlofs_num,
    %---------------------------------------------------------------------------
  % end of case Signal_Power in switch SimMod.Scan_mode
  %-----------------------------------------------------------------------------
  case 'Relative_Phase', 
    fprintf(1,'with varying Relative_Phase ... \n');
    
    for j1 = 1 : SimMod.Wlofs_num,   
      
      % Loop j1 for Cen_Wavelength
      
      ExcSet.Lambda_s = NLM.ZDW + ExcSet.Wvlen_ofs(j1);
      
      %-------------------------------------------------------------------------
      for j2 = 1 : SimMod.WlSep_num,

        % Loop j2 for PP_Separation

        % Calculate beta_vector and delta_beta(22,1) by Freq detuning
        [beta_matrix(:,j2),ExcSet.PP_Sp_Wvl(j2)] = f_delta_beta(ExcSet.Lambda_s,ExcSet.Sc_PPS_fq(j2)); 
        ExcSet.PP_Sp_Wvl(j2) = ExcSet.PP_Sp_Wvl(j2)*2e9;
        d_kappa_Lin(:,j2)    = delta_beta; 

        %-----------------------------------------------------------------------
        for j3 = 1 : SimMod.Power_num,
          
          % Loop j3 for initial Pump/Signal_Power
          
          for j4 = 1 : SimMod.Phase_num,

            % Loop j4 for initial Relative_Phase
            
            %-------------------------------------------------------------------
            for j5 = 1 : 1,
              
              % Create fields
              f_Create_Fld(j1,j2,j3,j4,j5);

              % Calculate ODEs of 7-wave
              if SimMod.FLD_nPWR == 0,
                %% Calculating the propagation by solving ODEs of Power/Phase evolution
                if SimMod.PSA_nPIA == 0,
                  % Calculating in PIA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PSA case']);         

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                   Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                   Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                   Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                   Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                   options);  % Initial Conditions    

                else % SimMod.PSA_nPIA == 1,
                  % Calculating in PSA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PIA case']); 

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                   options);  % Initial Conditions

                end; % end of if-else SimMod.PSA_nPIA

                Pwr_Var(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
              % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
              %---------------------------------------------------------------------------
              else % SimMod.FLD_nPWR == 1,
                %% Calculating the propagation by solving ODEs of Field evolution
                disp(['Calculating PSA @',num2str(ExcSet.Lambda_s*1e9,'%10.3f'), ...
                      'nm with PP=',num2str(ExcSet.PP_Sp_Wvl(j2),'%10.3f'), ...
                      'nm, P_p=',num2str(10*log10(Exc.Pwr(j3,1))+30,'%10.2f'), ...
                      'dBm, P_s=',num2str(10*log10(Exc.Pwr(j3,3))+30,'%10.2f'), ...
                      'dBm, using signal phase ', num2str(Exc.Phi(j4,3),'%10.3f'),...
                      'rad by 7w field']); 

                option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                [z,A]  = ode45(@dA_dz_7w, [0,NLM.Len], ...
                               [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                               option); % Initial Conditions

              end; % end of if-else SimMod.FLD_nPWR
              %-----------------------------------------------------------------
            end; % end of for j5 = 1 : SimMod.PropZ_num,
            %-------------------------------------------------------------------

            Raw.F_P1_Lin(j2,j4) =       A(end, 1);
            Raw.F_P2_Lin(j2,j4) =       A(end, 2);
            Raw.F_SI_Lin(j2,j4) =       A(end, 3);
            Raw.F_P3_Lin(j2,j4) =       A(end, 4);
            Raw.F_P4_Lin(j2,j4) =       A(end, 5);
            Raw.F_S1_Lin(j2,j4) =       A(end, 6);
            Raw.F_S2_Lin(j2,j4) =       A(end, 7);

            Raw.P_P1_Lin(j2,j4) =  (abs(A(end, 1))).^2;
            Raw.P_P2_Lin(j2,j4) =  (abs(A(end, 2))).^2;
            Raw.P_SI_Lin(j2,j4) =  (abs(A(end, 3))).^2;
            Raw.P_P3_Lin(j2,j4) =  (abs(A(end, 4))).^2;
            Raw.P_P4_Lin(j2,j4) =  (abs(A(end, 5))).^2;
            Raw.P_S1_Lin(j2,j4) =  (abs(A(end, 6))).^2;
            Raw.P_S2_Lin(j2,j4) =  (abs(A(end, 7))).^2;

            Raw.Phase_P1(j2,j4) = angle(A(end,1));                   % rad,
            Raw.Phase_P2(j2,j4) = angle(A(end,2));                   % rad,
            Raw.Phase_SI(j2,j4) = angle(A(end,3));                   % rad,
            Raw.Phase_P3(j2,j4) = angle(A(end,4));                   % rad,
            Raw.Phase_P4(j2,j4) = angle(A(end,5));                   % rad,
            Raw.Phase_S1(j2,j4) = angle(A(end,6));                   % rad,
            Raw.Phase_S2(j2,j4) = angle(A(end,7));                   % rad,

            Raw.P_P1_Log(j2,j4) = 10*log10( Raw.P_P1_Lin(j2,j4).*1000 ); % dBm,
            Raw.P_P2_Log(j2,j4) = 10*log10( Raw.P_P2_Lin(j2,j4).*1000 ); % dBm,
            Raw.P_SI_Log(j2,j4) = 10*log10( Raw.P_SI_Lin(j2,j4).*1000 ); % dBm,
            Raw.P_P4_Log(j2,j4) = 10*log10( Raw.P_P4_Lin(j2,j4).*1000 ); % dBm,
            Raw.P_P3_Log(j2,j4) = 10*log10( Raw.P_P3_Lin(j2,j4).*1000 ); % dBm,
            Raw.P_S1_Log(j2,j4) = 10*log10( Raw.P_S1_Lin(j2,j4).*1000 ); % dBm,
            Raw.P_S2_Log(j2,j4) = 10*log10( Raw.P_S2_Lin(j2,j4).*1000 ); % dBm,

            Raw.G_SI_Lin(j2,j4) = Raw.P_SI_Lin(j2,j4)./Exc.Pwr(j3,3);    % Rt ,
            Raw.G_SI_Log(j2,j4) = 10*log10( Raw.G_SI_Lin(j2,j4) );       % dB ,
            
          end; % end of for j4 = 1 : SimMod.Phase_num,
          %---------------------------------------------------------------------
          
        end; % end of for j3 = 1 : SimMod.Power_num,
        %-----------------------------------------------------------------------
        
       [Max.G_SI_max(j2, 1),Max.G_SI_max(j2, 2)] = max(Raw.G_SI_Log(j2,:));
        Max.G_SI_max(j2, 3) = Exc.Phi(Max.G_SI_max(j2, 2),3);
        
       [Max.G_SI_min(j2, 1),Max.G_SI_min(j2, 2)] = min(Raw.G_SI_Log(j2,:));
        Max.G_SI_min(j2, 3) = Exc.Phi(Max.G_SI_min(j2, 2),3);
        
        Max.G_SI_zro_m(j2,j1) = Raw.G_SI_Log(j2,1);
        Max.G_SI_max_m(j2,j1) = Max.G_SI_max(j2,1);
        Max.G_SI_min_m(j2,j1) = Max.G_SI_min(j2,1);
        Max.G_SI_PSER_m(j2,j1) = Max.G_SI_max(j2,1)-Max.G_SI_min(j2,1);
        Max.G_SI_PSGA_m(j2,j1) = Max.G_SI_max(j2,1)+Max.G_SI_min(j2,1);
        Max.G_SI_max_Ph_m(j2,j1) = Max.G_SI_max(j2,3);
        Max.G_SI_min_Ph_m(j2,j1) = Max.G_SI_min(j2,3);
          
      end; % end of for j2 = 1 : SimMod.WlSep_num,
      %-------------------------------------------------------------------------
      
      save([num2str(ExcSet.Lambda_s),'.mat']);
    end; % end of for j1 = 1 : SimMod.Wlofs_num,
    %---------------------------------------------------------------------------
  % end of case 'Propagation_Z' of switch SimMod.Scan_mode
  %-----------------------------------------------------------------------------
  case 'Propagation_Z',  
    fprintf(1,'with varying Propagation_Z ... \n');
    
    for j1 = 1 : SimMod.Wlofs_num, 
      
      % Loop j1 for Cen_Wavelength
      
      ExcSet.Lambda_s = NLM.ZDW + ExcSet.Wvlen_ofs(j1);
      
      %-------------------------------------------------------------------------
      for j2 = 1 : SimMod.WlSep_num,

        % Loop j2 for PP_Separation

        % Calculate beta_vector and delta_beta(22,1) by Freq detuning
        [beta_matrix(:,j2),ExcSet.Sc_PPS_wl(j2)] = f_delta_beta(ExcSet.Lambda_s,ExcSet.Dl_PPS_fq(j2)); 
        ExcSet.Sc_PPS_wl(j2) = ExcSet.Sc_PPS_wl(j2)*2e9;
        d_kappa_Lin(:,j2)    = delta_beta;  

        %-----------------------------------------------------------------------
        for j3 = 1 : SimMod.Power_num,
          
          % Loop j3 for initial Pump/Signal_Power
          
          for j4 = 1 : SimMod.Phase_num,

            % Loop j4 for initial Relative_Phase 
            
            switch SimMod.Scan_cond, 
              
              case 'Propagate', 
                %---------------------------------------------------------------
                for j5 = 1 : 1,

                  f_Create_Fld(j1,j2,j3,j4,j5);

                  % Calculate ODEs of 7-wave
                  if SimMod.FLD_nPWR == 0,
                    %% Calculating the propagation by solving ODEs of Power/Phase evolution
                    if SimMod.PSA_nPIA == 0,
                      % Calculating in PIA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PSA case']);         

                      options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                       Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                       Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                       Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                       Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                       options);  % Initial Conditions    

                    else % SimMod.PSA_nPIA == 1,
                      % Calculating in PSA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PIA case']); 

                      options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                       options);  % Initial Conditions

                    end; % end of if-else SimMod.PSA_nPIA

                    Var.Pwr(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
                  % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
                  %-------------------------------------------------------------
                  else % SimMod.FLD_nPWR == 1,
                    %% Calculating the propagation by solving ODEs of Field evolution
                    disp(['Calculating PSA @',...
                      'Cen=', num2str(ExcSet.Lambda_s*1e9,         '%10.3f'),'nm ',...
                      'PPS=', num2str(ExcSet.Dl_PPS_wl(j2),        '%10.3f'),'nm ',...
                      'P_p=', num2str(10*log10(Exc.Pwr(j3,1))+30+3,'%10.2f'),'dBm ',...
                      'P_s=', num2str(10*log10(Exc.Pwr(j3,3))+30  ,'%10.2f'),'dBm ',...
                      'Ph_s=',num2str(Exc.Phi(j4,3)/pi,            '%10.3f'),'rad ',...
                      'Len=', num2str(ExcSet.PropZ_Len(end),       '%10.1f'),'m ',...
                      'by 7w field']);

                    option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                    [z,A]  = ode45 (@dA_dz_7w, [0,ExcSet.PropZ_Len], ...
                                    [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                     Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                                    option); % Initial Conditions

                    Raw.F_P1_Lin(j4,:) =       A(end, 1);
                    Raw.F_P2_Lin(j4,:) =       A(end, 2);
                    Raw.F_SI_Lin(j4,:) =       A(end, 3);
                    Raw.F_P3_Lin(j4,:) =       A(end, 4);
                    Raw.F_P4_Lin(j4,:) =       A(end, 5);
                    Raw.F_S1_Lin(j4,:) =       A(end, 6);
                    Raw.F_S2_Lin(j4,:) =       A(end, 7);

                    Raw.P_P1_Lin(j4,:) =  (abs(A(end, 1))).^2;
                    Raw.P_P2_Lin(j4,:) =  (abs(A(end, 2))).^2;
                    Raw.P_SI_Lin(j4,:) =  (abs(A(end, 3))).^2;
                    Raw.P_P3_Lin(j4,:) =  (abs(A(end, 4))).^2;
                    Raw.P_P4_Lin(j4,:) =  (abs(A(end, 5))).^2;
                    Raw.P_S1_Lin(j4,:) =  (abs(A(end, 6))).^2;
                    Raw.P_S2_Lin(j4,:) =  (abs(A(end, 7))).^2;

                    Raw.Phase_P1(j4,:) = angle(A(end,1));                   % rad,
                    Raw.Phase_P2(j4,:) = angle(A(end,2));                   % rad,
                    Raw.Phase_SI(j4,:) = angle(A(end,3));                   % rad,
                    Raw.Phase_P3(j4,:) = angle(A(end,4));                   % rad,
                    Raw.Phase_P4(j4,:) = angle(A(end,5));                   % rad,
                    Raw.Phase_S1(j4,:) = angle(A(end,6));                   % rad,
                    Raw.Phase_S2(j4,:) = angle(A(end,7));                   % rad,

                    Raw.P_P1_Log(j4,:) = 10*log10( Raw.P_P1_Lin(j4,:).*1000 ); % dBm, 
                    Raw.P_P2_Log(j4,:) = 10*log10( Raw.P_P2_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_SI_Log(j4,:) = 10*log10( Raw.P_SI_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_P3_Log(j4,:) = 10*log10( Raw.P_P3_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_P4_Log(j4,:) = 10*log10( Raw.P_P4_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_S1_Log(j4,:) = 10*log10( Raw.P_S1_Lin(j4,:).*1000 ); % dBm,
                    Raw.P_S2_Log(j4,:) = 10*log10( Raw.P_S2_Lin(j4,:).*1000 ); % dBm,

                    Raw.G_SI_Lin(j4,:) = Raw.P_SI_Lin(j4,:)./Exc.Pwr(j3,3);    % Rt ,
                    Raw.G_SI_Log(j4,:) = 10*log10( Raw.G_SI_Lin(j4,:) );       % dB ,

                  end; % end of if-else SimMod.FLD_nPWR
                  %-------------------------------------------------------------

                end; % end of for j5 = 1 : 1 in 'Propagate',
                %---------------------------------------------------------------
              % end of case 'Propagate' in switch SimMod.Scan_cond
              %-----------------------------------------------------------------
              
              case 'Len_Z_var', 
                
                %-----------------------------------------------------------------
                for j5 = 2 : SimMod.PropZ_num,

                  f_Create_Fld(j1,j2,j3,j4,j5);

                  % calculate DE
                  if SimMod.FLD_nPWR == 0,
                    %% Calculating the propagation by solving ODEs of Power/Phase evolution
                    if SimMod.PSA_nPIA == 0,
                      % Calculating in PIA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PSA case']);         

                      options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                          1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                       Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                       Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                       Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                       Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                       options);  % Initial Conditions    

                    else % SimMod.PSA_nPIA == 1,
                      % Calculating in PSA mode
                      disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                            'nm using Power in 7-wave PIA case']); 

                      options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                      [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                       [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                       options);  % Initial Conditions

                    end; % end of if-else SimMod.PSA_nPIA

                    Var.Pwr(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
                  % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
                  %---------------------------------------------------------------------------
                  else % SimMod.FLD_nPWR == 1,
                    %% Calculating the propagation by solving ODEs of Field evolution
                    disp(['Calculating PSA @',num2str(ExcSet.Lambda_s*1e9,'%10.3f'), ...
                          'nm with PP=',num2str(ExcSet.PP_Sp_Wvl(j2),'%10.3f'), ...
                          'nm, P_p=',num2str(10*log10(Exc.Pwr(j3,1))+30,'%10.2f'), ...
                          'dBm, P_s=',num2str(10*log10(Exc.Pwr(j3,3))+30,'%10.2f'), ...
                          'dBm, using signal phase ', num2str(Exc.Phi(j4,3)/pi,'%10.3f'),...
                          'rad, in ', num2str(ExcSet.PropZ_Len(j5),'%10.3f'),...
                          'm, by 7w field']);

                    option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                    [z,A]  = ode45 (@dA_dz_7w, [0 ExcSet.PropZ_Len(j5)], ...
                                    [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                     Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                                    option); % Initial Conditions

                    Raw.F_P1_Lin(j4,j5) =       A(end, 1);
                    Raw.F_P2_Lin(j4,j5) =       A(end, 2);
                    Raw.F_SI_Lin(j4,j5) =       A(end, 3);
                    Raw.F_P3_Lin(j4,j5) =       A(end, 4);
                    Raw.F_P4_Lin(j4,j5) =       A(end, 5);
                    Raw.F_S1_Lin(j4,j5) =       A(end, 6);
                    Raw.F_S2_Lin(j4,j5) =       A(end, 7);

                    Raw.P_P1_Lin(j4,j5) =  (abs(A(end, 1))).^2;
                    Raw.P_P2_Lin(j4,j5) =  (abs(A(end, 2))).^2;
                    Raw.P_SI_Lin(j4,j5) =  (abs(A(end, 3))).^2;
                    Raw.P_P3_Lin(j4,j5) =  (abs(A(end, 4))).^2;
                    Raw.P_P4_Lin(j4,j5) =  (abs(A(end, 5))).^2;
                    Raw.P_S1_Lin(j4,j5) =  (abs(A(end, 6))).^2;
                    Raw.P_S2_Lin(j4,j5) =  (abs(A(end, 7))).^2;

                    Raw.Phase_P1(j4,j5) = angle(A(end,1));                   % rad,
                    Raw.Phase_P2(j4,j5) = angle(A(end,2));                   % rad,
                    Raw.Phase_SI(j4,j5) = angle(A(end,3));                   % rad,
                    Raw.Phase_P3(j4,j5) = angle(A(end,4));                   % rad,
                    Raw.Phase_P4(j4,j5) = angle(A(end,5));                   % rad,
                    Raw.Phase_S1(j4,j5) = angle(A(end,6));                   % rad,
                    Raw.Phase_S2(j4,j5) = angle(A(end,7));                   % rad,

                    Raw.P_P1_Log(j4,j5) = 10*log10( Raw.P_P1_Lin(j4,j5).*1000 ); % dBm, 
                    Raw.P_P2_Log(j4,j5) = 10*log10( Raw.P_P2_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_SI_Log(j4,j5) = 10*log10( Raw.P_SI_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_P3_Log(j4,j5) = 10*log10( Raw.P_P3_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_P4_Log(j4,j5) = 10*log10( Raw.P_P4_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_S1_Log(j4,j5) = 10*log10( Raw.P_S1_Lin(j4,j5).*1000 ); % dBm,
                    Raw.P_S2_Log(j4,j5) = 10*log10( Raw.P_S2_Lin(j4,j5).*1000 ); % dBm,

                    Raw.G_SI_Lin(j4,j5) = Raw.P_SI_Lin(j4,j5)./Exc.Pwr(j3,3);    % Rt ,
                    Raw.G_SI_Log(j4,j5) = 10*log10( Raw.G_SI_Lin(j4,j5) );       % dB ,

                  end; % end of if-else SimMod.FLD_nPWR
                  %-------------------------------------------------------------

                end; % end of for j5 = 1 : SimMod.PropL_num,
                %---------------------------------------------------------------
              % end of case 'Len_Z_var' in switch SimMod.Scan_cond
              %-----------------------------------------------------------------
              
            end; % end of switch SimMod.Scan_cond
            %-------------------------------------------------------------------

          end; % end of for j4 = 1 : SimMod.Phase_num,
          %---------------------------------------------------------------------

        end; % end of for j3 = 1 : SimMod.Power_num,
        %-----------------------------------------------------------------------

       [Max.G_SI_max(j2, 1),Max.G_SI_max(j2, 2)] = max(Raw.G_SI_Log(:,end)); % Max gain and corresponding index
        Max.G_SI_max(j2, 3) = Exc.Phi(Max.G_SI_max(j2, 2),3);                % relative phase for max gain
         
       [Max.G_SI_min(j2, 1),Max.G_SI_min(j2, 2)] = min(Raw.G_SI_Log(:,end)); % Min gain and corresponding index
        Max.G_SI_min(j2, 3) = Exc.Phi(Max.G_SI_min(j2, 2),3);                % relative phase for min gain
        
        Max.G_SI_max_z(j2, :) = Raw.G_SI_Log(Max.G_SI_max(j2, 2),:);
        
        Max.P_P1_max_z(j2, :) = Raw.P_P1_Log(Max.G_SI_max(j2, 2),:);
        Max.P_P2_max_z(j2, :) = Raw.P_P2_Log(Max.G_SI_max(j2, 2),:);
        Max.P_SI_max_z(j2, :) = Raw.P_SI_Log(Max.G_SI_max(j2, 2),:);
        Max.P_P3_max_z(j2, :) = Raw.P_P3_Log(Max.G_SI_max(j2, 2),:);
        Max.P_P4_max_z(j2, :) = Raw.P_P4_Log(Max.G_SI_max(j2, 2),:);
        Max.P_S1_max_z(j2, :) = Raw.P_S1_Log(Max.G_SI_max(j2, 2),:);
        Max.P_S2_max_z(j2, :) = Raw.P_S2_Log(Max.G_SI_max(j2, 2),:);
        
        Max.G_SI_min_z(j2, :) = Raw.G_SI_Log(Max.G_SI_min(j2, 2),:);
        
        Max.P_P1_min_z(j2, :) = Raw.P_P1_Log(Max.G_SI_min(j2, 2),:);
        Max.P_P2_min_z(j2, :) = Raw.P_P2_Log(Max.G_SI_min(j2, 2),:);
        Max.P_SI_min_z(j2, :) = Raw.P_SI_Log(Max.G_SI_min(j2, 2),:);
        Max.P_P3_min_z(j2, :) = Raw.P_P3_Log(Max.G_SI_min(j2, 2),:);
        Max.P_P4_min_z(j2, :) = Raw.P_P4_Log(Max.G_SI_min(j2, 2),:);
        Max.P_S1_min_z(j2, :) = Raw.P_S1_Log(Max.G_SI_min(j2, 2),:);
        Max.P_S2_min_z(j2, :) = Raw.P_S2_Log(Max.G_SI_min(j2, 2),:);

%         Max.G_SI_zro_m(j2,j1) = Raw.G_SI_Log(j2,1);
%         Max.G_SI_max_m(j2,j1) = Max.G_SI_max(j2,1);
%         Max.G_SI_min_m(j2,j1) = Max.G_SI_min(j2,1);
%         Max.G_SI_PSER_m(j2,j1) = Max.G_SI_max(j2,1)-Max.G_SI_min(j2,1);
%         Max.G_SI_PSGA_m(j2,j1) = Max.G_SI_max(j2,1)+Max.G_SI_min(j2,1);
%         Max.G_SI_max_Ph_m(j2,j1) = Max.G_SI_max(j2,3);
%         Max.G_SI_min_Ph_m(j2,j1) = Max.G_SI_min(j2,3);
        
      end; % end of for j2 = 1 : SimMod.WlSep_num,
      %-------------------------------------------------------------------------
        
    % save([num2str(ExcSet.Lambda_s),'.mat']);
      
    end; % end of for j1 = 1 : SimMod.Wlofs_num,
    %---------------------------------------------------------------------------
  % end of case 'Propagation_Z' in switch SimMod.Scan_mode
  %-----------------------------------------------------------------------------
  case 'Constellation',  
    fprintf(1,'with varying Signal_Power ... \n');
    
    for j1 = 1 : SimMod.Wlofs_num,
      
      % Loop j1 for Cen_Wavelength
      
      ExcSet.Lambda_s = NLM.ZDW + ExcSet.Wvlen_ofs(j1);
      
      %-------------------------------------------------------------------------
      for j2 = 1 : SimMod.WlSep_num,

        % Loop j2 for PP_Separation

        % Calculate beta_vector and delta_beta(22,1) by Freq detuning
        [beta_matrix(:,j2),ExcSet.Sc_PPS_wl(j2)] = f_delta_beta(ExcSet.Lambda_s,ExcSet.Dl_PPS_fq(j2)); 
        ExcSet.dL_PPS_wl(j2) = ExcSet.Sc_PPS_wl(j2)*2e9; % Transfer SI to nm
        d_kappa_Lin(:,j2)    = delta_beta;  

        %-----------------------------------------------------------------------
        for j3 = 1 : SimMod.Power_num,
          
          % Loop j3 for initial Pump/Signal_Power
          
          for j4 = 1 : SimMod.Noise_num,

            % Loop j4 for initial Relative_Phase
            
            %-------------------------------------------------------------------
            for j5 = 1 : 1,
              
              % Create fields
              f_Nis_Gen(SimMod.Noise_num,j3); % 
              f_Create_Fld(j1,j2,j3,j4,j5);   % 
              
              % Calculate ODEs of 7-wave
              if SimMod.FLD_nPWR == 0,
                %% Calculating the propagation by solving ODEs of Power/Phase evolution
                if SimMod.PSA_nPIA == 0,
                  % Calculating in PIA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PSA case']);         

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 ...
                                      1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dE_dz_7wave, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2 ...
                                   Theta( 1) Theta( 2) Theta( 3) Theta( 4) Theta( 5) Theta( 6) ...
                                   Theta( 7) Theta( 8) Theta( 9) Theta(10) Theta(11) Theta(12) ...
                                   Theta(13) Theta(14) Theta(15) Theta(16) Theta(17) Theta(18) ...
                                   Theta(19) Theta(20) Theta(21) Theta(22) ], ...
                                   options);  % Initial Conditions    

                else % SimMod.PSA_nPIA == 1,
                  % Calculating in PSA mode
                  disp(['Calculate at Pump-Signal separation: ',num2str(L_Detuning(i1),'%10.3f'), ...
                        'nm using Power in 7-wave PIA case']); 

                  options = odeset('RelTol', 1e-6, 'AbsTol', [1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                  [z,P]   = ode45 (@dA_dz_7w, [0,NLM.Len], ...
                                   [A_P1 A_P2 A_SI A_P4 A_P3 A_S1 A_S2], ...
                                   options);  % Initial Conditions

                end; % end of if-else SimMod.PSA_nPIA

                Pwr_Var(:,j2) = [P(end,1) P(end,2) P(end,3) P(end,4) P(end,5) P(end,6) P(end,7)];
              % Pwr_Var_Log   = 10*log10( P_Var.*1000 );
              %---------------------------------------------------------------------------
              else % SimMod.FLD_nPWR == 1,
                %% Calculating the propagation by solving ODEs of Field evolution
                disp(['Calculating PSA @',num2str(ExcSet.Lambda_s*1e9,'%10.3f'), ...
                      'nm with PPS=',num2str(ExcSet.dL_PPS_wl(j2),'%10.3f'), ...
                      'nm, P_p=',   num2str(10*log10(abs(Exc.Fld(1,1))^2)+30+3,'%10.2f'), ...
                      'dBm, P_s=',  num2str(10*log10(abs(Exc.Fld(1,3))^2)+30+0,'%10.2f'), ...
                      'dBm, using signal phase ', num2str(angle(Exc.Fld(1,3))/pi,'%10.3f'),...
                      'rad by 7w field']); 

                option = odeset('RelTol',1e-9, 'AbsTol',[1e-12 1e-12 1e-12 1e-12 1e-12 1e-12 1e-12]); % 
                [z,A]  = ode45(@dA_dz_7w, [0 ExcSet.PropZ_Len], ...
                               [Exc.Fld(1,1) Exc.Fld(1,2) Exc.Fld(1,3) Exc.Fld(1,4)...
                                Exc.Fld(1,5) Exc.Fld(1,6) Exc.Fld(1,7)],...
                               option); % Initial Conditions
                               
                Raw.F_P1_Lin(j3,j4) =       A(end, 1);
                Raw.F_P2_Lin(j3,j4) =       A(end, 2);
                Raw.F_SI_Lin(j3,j4) =       A(end, 3);
                Raw.F_P3_Lin(j3,j4) =       A(end, 4);
                Raw.F_P4_Lin(j3,j4) =       A(end, 5);
                Raw.F_S1_Lin(j3,j4) =       A(end, 6);
                Raw.F_S2_Lin(j3,j4) =       A(end, 7);

                Raw.P_P1_Lin(j3,j4) =  (abs(A(end, 1))).^2;
                Raw.P_P2_Lin(j3,j4) =  (abs(A(end, 2))).^2;
                Raw.P_SI_Lin(j3,j4) =  (abs(A(end, 3))).^2;
                Raw.P_P3_Lin(j3,j4) =  (abs(A(end, 4))).^2;
                Raw.P_P4_Lin(j3,j4) =  (abs(A(end, 5))).^2;
                Raw.P_S1_Lin(j3,j4) =  (abs(A(end, 6))).^2;
                Raw.P_S2_Lin(j3,j4) =  (abs(A(end, 7))).^2;

                Raw.Phase_P1(j3,j4) = angle(A(end,1));                   % rad,
                Raw.Phase_P2(j3,j4) = angle(A(end,2));                   % rad,
                Raw.Phase_SI(j3,j4) = angle(A(end,3));                   % rad,
                Raw.Phase_P3(j3,j4) = angle(A(end,4));                   % rad,
                Raw.Phase_P4(j3,j4) = angle(A(end,5));                   % rad,
                Raw.Phase_S1(j3,j4) = angle(A(end,6));                   % rad,
                Raw.Phase_S2(j3,j4) = angle(A(end,7));                   % rad,

              end; % end of if-else SimMod.FLD_nPWR
              %-----------------------------------------------------------------
            end; % end of for j5 = 1 : SimMod.PropZ_num,
            %-------------------------------------------------------------------

          end; % end of for j4 = 1 : SimMod.Phase_num,
          %---------------------------------------------------------------------
          
        end; % end of for j3 = 1 : SimMod.Power_num,
        %-----------------------------------------------------------------------

      end; % end of for j2 = 1 : SimMod.WlSep_num,
      %-------------------------------------------------------------------------
      
    end; % end of for j1 = 1 : SimMod.Wlofs_num,
    %---------------------------------------------------------------------------
  % end of case Signal_Power in switch SimMod.Scan_mode
  %----------------------------------------------------------------------------- 
  otherwise, 
    %
  %-----------------------------------------------------------------------------
end; % 

fprintf(1,'[COMPLETED]\n');
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------

Time_Elapsed = toc/60;
disp(['Elapsed time: ',num2str(Time_Elapsed,'%10.2f'),'minutes']);

%-------------------------------------------------------------------------------

%% Storage of the results 

if (SimMod.Wlofs_num == 1) && (SimMod.WlSep_num == 1),
  save([   num2str(SimMod.Wave_num, '%1.0f'),'w_',...
           SimMod.Scan_mode,...
    '_ofs',num2str(ExcSet.Wvlen_ofs*1e9,'%3.1f'),...
    '_pps',num2str(ExcSet.Dl_PPS_wl,'%3.1f'),...
    '_pwr',num2str(SimMod.Power_num,'%3.0f'),...
    '_phs',num2str(SimMod.Phase_num,'%3.0f'),...
    '_z',  num2str(SimMod.Len_Mul,  '%2.0f'),...
    'x',   num2str(SimMod.PropZ_num,'%4.0f'),...
    '.mat']); 
else 
  save([   num2str(SimMod.Wave_num, '%1.0f'),'w_',...
           SimMod.Scan_mode,...
    '_ofs',num2str(SimMod.Wlofs_num,'%3.0f'),...
    '_pps',num2str(SimMod.WlSep_num,'%3.0f'),...
    '_pwr',num2str(SimMod.Power_num,'%3.0f'),...
    '_phs',num2str(SimMod.Phase_num,'%3.0f'),...
    '_z',  num2str(SimMod.Len_Mul,  '%2.0f'),...
    'x',   num2str(SimMod.PropZ_num,'%4.0f'),...
    '.mat']); 
end; 
%-------------------------------------------------------------------------------

%% Visualization of the results 
f_Visualization( ...
  1, ...% sig_g,    Fig.01 Signal gain vs. PP_Separation
  1, ...% pwr_evo,  Fig.02 Power evolution vs. PP_Separation
  1, ...% kappa,    Fig.03 Total phase misExcSet.Dl_PPS_wlmatch vs. PP_Separation
  1, ...% sig_g_NPS,Fig.01 Signal gain vs. nonlineariy phase shift (rPL)
  1, ...% PSER_PSGA,Fig.06 Signal gain vs. Nonlineariy phase shift (rPL)
  1, ...% PtA_TC,   Fig.05 Phase-to-amplitude transfer function
  1, ...% PtP_TC,   Fig.07 Phase-to-phase transfer function
  1  ...% PtT_TC_c, Fig.08 Complex output signal trajectory vs. Input signal phase
  );
%-------------------------------------------------------------------------------

%% End of file 
%-------------------------------------------------------------------------------
% end of file Prop_CWODE_7w.m
%-------------------------------------------------------------------------------
%% using Git

!git init
!git add --all
!git commit -m "xdd0707"
!git push origin master
