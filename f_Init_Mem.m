%-------------------------------------------------------------------------------
% Beginning of file f_Init_Mem.m
%-------------------------------------------------------------------------------
%
% Project : Calculation of nonlinear propagation in fiber for FOPA based PSA 
%           using 7-wave model
% FileName: f_Init_Mem.m
% Function: Initialize storage memories
% Version : 
%   v0.0 @2015.11.30 created  by W.L. Xie 
%        Change of the architecture of the program, create this code script
%
%   v0.1 @2015.03.31 released by W.L. Xie
%        2015.09.15 Change HNLF parameters into classes HNLF.*
% Describe: Initialize the variables for simulation
%   
%-------------------------------------------------------------------------------

function [ ] = f_Init_Mem(  )
% j1: Center_Wavelength
% j2: PP_Separation
% j3: Pump/Signal_Power
% j4: Relative_Phase
% j5: Propagation_Len

global SimMod;

global Exc;  
% Exc.Pwr
% Exc.Phi
% Exc.Fld

global ExcSet;
% ExcSet.Lambda_s
% ExcSet.Lambda_c
% ExcSet.Sc_PP_Frq
% ExcSet.PP_Sp_Frq
% ExcSet.Sc_PP_Wvl
% ExcSet.PP_Sp_Wvl

%% global variables 

global Var;

% Var.Fld;   %
% Var.Phi;   %
% Var.Pwr;   % 

global Raw;

% Raw.F_P1_Lin;  % Field of P1 in linear scale
% Raw.F_P2_Lin;  %
% Raw.F_SI_Lin;  %
% Raw.F_P3_Lin;  %
% Raw.F_P4_Lin;  %
% Raw.F_S1_Lin;  %
% Raw.F_S2_Lin;  %

% Raw.P_P1_Lin;  % Power of P1 in linear scale (W)
% Raw.P_P2_Lin;  %
% Raw.P_SI_Lin;  %
% Raw.P_P3_Lin;  %
% Raw.P_P4_Lin;  %
% Raw.P_S1_Lin;  %
% Raw.P_S2_Lin;  %

% Raw.P_P1_Log;  % Power of P1 in log scale (dBm)
% Raw.P_P2_Log;  %
% Raw.P_SI_Log;  %
% Raw.P_P3_Log;  %
% Raw.P_P4_Log;  %
% Raw.P_S1_Log;  %
% Raw.P_S2_Log;  %

% Raw.Phase_P1;  % Phase of P1 (rad)
% Raw.Phase_P2;  %
% Raw.Phase_SI;  %
% Raw.Phase_P3;  %
% Raw.Phase_P4;  %
% Raw.Phase_S1;  %
% Raw.Phase_S2;  %

% Raw.G_SI_Lin;  % Gain in linear scale
% Raw.G_SI_Log;  % Gain in log scale (dB)

global Max;

% Max.G_SI_max;  % (a+b)/2+((a-b)/2)*cos(2*(x-c))
% Max.G_SI_min;  % (Gmax+Gmin)/2+((Gmax-Gmin)/2)*cos(2*(x-PHmax))
% 1: Gain maxima
% 2: relative phase index corresponds to the gain maxima 
% 3: relative phase value corresponds to the gain maxima 

% Max.G_SI_zro_m;
% Max.G_SI_max_m;
% Max.G_SI_min_m;
% Max.G_SI_PSER_m;   % Phase-sensitive extinction ratio
% Max.G_SI_PSGA_m;   % Phase-sensitive gain asymmetry
% Max.G_SI_max_Ph_m; % 
% Max.G_SI_min_Ph_m; % 

%% Initial memory   

switch SimMod.Scan_mode,
  case 'Cen_Wavelength',
    
    Var.Fld      = zeros(7,               SimMod.Phase_num);
    Var.Pwr      = zeros(7,               SimMod.Phase_num);
    Var.Phi      = zeros(7,               SimMod.Phase_num);
    
    Raw.F_P1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P3_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P4_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_S1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_S2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.P_P1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P3_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P4_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.P_P1_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P2_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_SI_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P3_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P4_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S1_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S2_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.Phase_P1 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P2 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_SI = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P3 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P4 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_S1 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_S2 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.G_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.G_SI_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Max.G_SI_max = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_min = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    
    Max.G_SI_zro_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_max_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_min_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_PSER_m   = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_PSGA_m   = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_max_Ph_m = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_min_Ph_m = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
  % end of case Signal_Power in switch SimMod.Scan_mode
  case 'PP_Separation', 
    
    Var.Fld      = zeros(7,               SimMod.Phase_num);
    Var.Pwr      = zeros(7,               SimMod.Phase_num);
    Var.Phi      = zeros(7,               SimMod.Phase_num);
    
    Raw.F_P1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P3_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P4_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_S1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_S2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.P_P1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P3_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P4_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.P_P1_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P2_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_SI_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P3_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P4_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S1_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S2_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.Phase_P1 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P2 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_SI = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P3 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P4 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_S1 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_S2 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.G_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.G_SI_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Max.G_SI_max = zeros(SimMod.WlSep_num,               7);
    Max.G_SI_min = zeros(SimMod.WlSep_num,               7);
    
    Max.G_SI_zro_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_max_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_min_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_PSER_m   = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_PSGA_m   = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_max_Ph_m = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
    Max.G_SI_min_Ph_m = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
  % end of case Signal_Power in switch SimMod.Scan_mode
  case 'Pump_Power',    
    switch SimMod.Scan_cond, 
      case 'Normal',    
        Var.Fld        = zeros(7,               SimMod.Phase_num);
        Var.Phi        = zeros(7,               SimMod.Phase_num);
        Var.Pwr        = zeros(7,               SimMod.Phase_num);

        Raw.F_P1_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.F_P2_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.F_SI_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.F_P3_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.F_P4_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.F_S1_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.F_S2_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);

        Raw.P_P1_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.P_P2_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.P_SI_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.P_P3_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.P_P4_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.P_S1_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.P_S2_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);

        Raw.Phase_P1   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.Phase_P2   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.Phase_SI   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.Phase_P3   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.Phase_P4   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.Phase_S1   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.Phase_S2   = zeros(SimMod.Power_num,SimMod.Phase_num);

        Raw.G_SI_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
        Raw.G_SI_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);

        Max.G_SI_max   = zeros(SimMod.Power_num,               7);
        Max.G_SI_min   = zeros(SimMod.Power_num,               7);

        Max.G_SI_zro_m    = zeros(SimMod.Power_num,SimMod.Wlofs_num);
        Max.G_SI_max_m    = zeros(SimMod.Power_num,SimMod.Wlofs_num);
        Max.G_SI_min_m    = zeros(SimMod.Power_num,SimMod.Wlofs_num);
        Max.G_SI_PSER_m   = zeros(SimMod.Power_num,SimMod.Wlofs_num);
        Max.G_SI_PSGA_m   = zeros(SimMod.Power_num,SimMod.Wlofs_num);
        Max.G_SI_max_Ph_m = zeros(SimMod.Power_num,SimMod.Wlofs_num);
        Max.G_SI_min_Ph_m = zeros(SimMod.Power_num,SimMod.Wlofs_num);

      % end of case 'Normal' in switch SimMod.Scan_cond
      %-------------------------------------------------------------------------
      case 'Propagate', 
        

      % end of case 'Propagate' in switch SimMod.Scan_cond
      %-------------------------------------------------------------------------
      case 'Len_Z_var', 
        

      % end of case 'Len_Z_var' in switch SimMod.Scan_cond
      %-------------------------------------------------------------------------
    end; % end of switch SimMod.Scan_cond
    %---------------------------------------------------------------------------
  % end of case Pump_Power in switch SimMod.Scan_mode
  %-----------------------------------------------------------------------------
  case 'Signal_Power',  
    
    Var.Fld        = zeros(7,               SimMod.Phase_num);
    Var.Phi        = zeros(7,               SimMod.Phase_num);
    Var.Pwr        = zeros(7,               SimMod.Phase_num);
    
    Raw.F_P1_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.F_P2_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.F_SI_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.F_P3_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.F_P4_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.F_S1_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.F_S2_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    
    Raw.P_P1_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.P_P2_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.P_SI_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.P_P3_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.P_P4_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.P_S1_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.P_S2_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    
    Raw.Phase_P1   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.Phase_P2   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.Phase_SI   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.Phase_P3   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.Phase_P4   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.Phase_S1   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.Phase_S2   = zeros(SimMod.Power_num,SimMod.Phase_num);
    
    Raw.G_SI_Lin   = zeros(SimMod.Power_num,SimMod.Phase_num);
    Raw.G_SI_Log   = zeros(SimMod.Power_num,SimMod.Phase_num);
    
    Max.G_SI_max   = zeros(SimMod.Power_num,               7);
    Max.G_SI_min   = zeros(SimMod.Power_num,               7);
    
    Max.G_SI_zro_m    = zeros(SimMod.Power_num,SimMod.Wlofs_num);
    Max.G_SI_max_m    = zeros(SimMod.Power_num,SimMod.Wlofs_num);
    Max.G_SI_min_m    = zeros(SimMod.Power_num,SimMod.Wlofs_num);
    Max.G_SI_PSER_m   = zeros(SimMod.Power_num,SimMod.Wlofs_num);
    Max.G_SI_PSGA_m   = zeros(SimMod.Power_num,SimMod.Wlofs_num);
    Max.G_SI_max_Ph_m = zeros(SimMod.Power_num,SimMod.Wlofs_num);
    Max.G_SI_min_Ph_m = zeros(SimMod.Power_num,SimMod.Wlofs_num);
  % end of case Signal_Power in switch SimMod.Scan_mode
  case 'Relative_Phase',
    
    Var.Fld      = zeros(7,               SimMod.Phase_num);
    Var.Phi      = zeros(7,               SimMod.Phase_num);
    Var.Pwr      = zeros(7,               SimMod.Phase_num);
    
    Raw.F_P1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P3_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_P4_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_S1_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.F_S2_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.P_P1_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P2_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_SI_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P3_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_P4_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S1_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.P_S2_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.Phase_P1 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P2 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_SI = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P3 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_P4 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_S1 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.Phase_S2 = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Raw.G_SI_Lin = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    Raw.G_SI_Log = zeros(SimMod.WlSep_num,SimMod.Phase_num);
    
    Max.G_SI_max   = zeros(SimMod.WlSep_num,               7);
    Max.G_SI_min   = zeros(SimMod.WlSep_num,               7);
    
    Max.G_SI_zro_m    = zeros(SimMod.WlSep_num,SimMod.WlSep_num);
    Max.G_SI_max_m    = zeros(SimMod.WlSep_num,SimMod.WlSep_num);
    Max.G_SI_min_m    = zeros(SimMod.WlSep_num,SimMod.WlSep_num);
    Max.G_SI_PSER_m   = zeros(SimMod.WlSep_num,SimMod.WlSep_num);
    Max.G_SI_PSGA_m   = zeros(SimMod.WlSep_num,SimMod.WlSep_num);
    Max.G_SI_max_Ph_m = zeros(SimMod.WlSep_num,SimMod.WlSep_num);
    Max.G_SI_min_Ph_m = zeros(SimMod.WlSep_num,SimMod.WlSep_num);
  % end of case Signal_Power in switch SimMod.Scan_mode
  case 'Propagation_Z', 
    
    Var.Fld        = zeros(7,               SimMod.Phase_num);
    Var.Phi        = zeros(7,               SimMod.Phase_num);
    Var.Pwr        = zeros(7,               SimMod.Phase_num);
    
    Raw.F_P1_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.F_P2_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.F_SI_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.F_P3_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.F_P4_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.F_S1_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.F_S2_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    
    Raw.P_P1_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_P2_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_SI_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_P3_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_P4_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_S1_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_S2_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    
    Raw.P_P1_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_P2_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_SI_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_P3_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_P4_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_S1_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.P_S2_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    
    Raw.Phase_P1   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.Phase_P2   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.Phase_SI   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.Phase_P3   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.Phase_P4   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.Phase_S1   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.Phase_S2   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    
    Raw.G_SI_Lin   = zeros(SimMod.Phase_num,SimMod.PropZ_num);
    Raw.G_SI_Log   = zeros(SimMod.Phase_num,SimMod.PropZ_num);

    Max.G_SI_max   = zeros(SimMod.WlSep_num,               3);
    Max.G_SI_min   = zeros(SimMod.WlSep_num,               3);
    
    Max.G_SI_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    
    Max.P_P1_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_P2_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_SI_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_P3_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_P4_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_S1_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_S2_max_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    
    Max.G_SI_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
     
    Max.P_P1_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_P2_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_SI_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_P3_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_P4_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_S1_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    Max.P_S2_min_z = zeros(SimMod.WlSep_num,SimMod.PropZ_num);
    
%     Max.G_SI_zro_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
%     Max.G_SI_max_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
%     Max.G_SI_min_m    = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
%     Max.G_SI_PSER_m   = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
%     Max.G_SI_PSGA_m   = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
%     Max.G_SI_max_Ph_m = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
%     Max.G_SI_min_Ph_m = zeros(SimMod.WlSep_num,SimMod.Wlofs_num);
  % end of case Propagation_Len in switch SimMod.Scan_mode
  otherwise,            
    %
end; % end of switch
%-------------------------------------------------------------------------------

end % end of function f_Init_Mem
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% end of file f_Init_Mem.m
%-------------------------------------------------------------------------------