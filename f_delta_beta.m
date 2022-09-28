%-------------------------------------------------------------------------------
% Beginning of file f_delta_beta.m
%-------------------------------------------------------------------------------
%
% Optical Parametric Amplification: 7-wave Analysis
% 
%                P1        P2
%                |         |
%      P4        |         |         P3
%      |         |   S&I   |         |
%      |    S1   |    |    |    S2   |
%      |    |    |    |    |    |    |
%      |    |    |    |    |    |    |
% ------------------------------------------> Wavelength(nm)
%      A4   A6   A1   A3   A2   A7   A5    
%      d    g    a    c    b    h    f
% 
% Data   : 2015.03.16
% Version: Calculation of dispersion parameters for general case
%          Beta_0, Beta_1, Beta_2, Beta_3 and Beta_4
% 

%%
function [beta_vector, Lambda_d] = f_delta_beta(Lambda_0, PP_Detuning)
% beta_vector: beta vector under Taylor Power Series expansion
% Lambda_d   : (Lambda_P2 - Lambda_P1)/2
% Lambda_0   : center wavelength of all waves
% PP_Detuning: pump-pump detuning in Frequency

%% Global constant 
global c;          %

%% Simulation paramters 
global SimMod;

%% Global varibles 
global NLM;      % HNLF parameters
%      NLM.SN
%      NLM.ZDW;  % m, zero dispersion wavelength
%      NLM.ZDF;  % Hz, zero dispersion frequency
%      NLM.D0;   % ps/(km*nm^1), HNLF dispersion
%      NLM.D1;   % ps/(km*nm^2), HNLF dispersion slope
%      NLM.D2;   % ps/(km*nm^3), slope of HNLF dispersion slope

global ExcSet;      % 

global beta_matrix; % 
global delta_beta;  % delta beta vector
global d_kappa_Lin; %

%% Local varibles 

%-------------------------------------------------------------------------------

%% Allocate wave parameters 

F_Detune  = PP_Detuning/2;       % Hz, frequency of wave detuning
Freq_0    = c/Lambda_0;          % Hz, frequency of center wave

Freq_P4   = Freq_0 + 3*F_Detune; % pi*Hz, angular frequency of f0
Freq_S1   = Freq_0 + 2*F_Detune; % pi*Hz, angular frequency of f0 
Freq_P1   = Freq_0 + 1*F_Detune; % pi*Hz, angular frequency of f0 
Freq_SI   = Freq_0;
Freq_P2   = Freq_0 - 1*F_Detune; % pi*Hz, angular frequency of f0 
Freq_S2   = Freq_0 - 2*F_Detune; % pi*Hz, angular frequency of f0 
Freq_P3   = Freq_0 - 3*F_Detune; % pi*Hz, angular frequency of f0 

Freq_c    = (Freq_P1 + Freq_P2)/2;% Hz, center frequency of P1 & P2
Freq_d    = (Freq_P1 - Freq_P2)/2;% Hz, frequency detuning among waves

Lambda_P4 = c/Freq_P4;           % m, wavelength of Pump_4
Lambda_S1 = c/Freq_S1;           % m, wavelength of Signal_1
Lambda_P1 = c/Freq_P1;           % m, wavelength of Pump_1
Lambda_SI = Lambda_0;            % m, wavelength of Signal&Idler
Lambda_P2 = c/Freq_P2;           % m, wavelength of Pump_2
Lambda_S2 = c/Freq_S2;           % m, wavelength of Signal_2  
Lambda_P3 = c/Freq_P3;           % m, wavelength of Pump_3  

Lambda_c  = (Lambda_P2 + Lambda_P1)/2;
Lambda_d  = (Lambda_P2 - Lambda_P1)/2;
%-------------------------------------------------------------------------------

%% HNLF parameter 

NLM_D0  = (NLM.D0);
NLM_D0  = (NLM.D1)*(Lambda_0 - NLM.ZDW); % s/m^2, HNLF dispersion
NLM_D1  = (NLM.D1);                       % s/m^3, HNLF dispersion slope
NLM_D2  = (NLM.D2);                       % s/m^4, slope of HNLF dispersion slope
%-------------------------------------------------------------------------------
                         
%% Calculate propagation constant/wavevector   
if NLM.SN == 19, % Photonic Crystal Waveguide  
  
  [Beta_0, Beta_1, Beta_2, Beta_3, Beta_4] = f_Cal_beta(Freq_SI);
  beta_vector = [Beta_0, Beta_1, Beta_2, Beta_3, Beta_4];
  
  Beta_P4= (1/1) *(((2*pi)^1)*((Freq_P4-Freq_0)^1))*Beta_1 ...
          +(1/2) *(((2*pi)^2)*((Freq_P4-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P4-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P4-Freq_0)^4))*Beta_4;  
  beta_d = Beta_P4;

  Beta_S1= (1/1) *(((2*pi)^1)*((Freq_S1-Freq_0)^1))*Beta_1 ...
          +(1/2) *(((2*pi)^2)*((Freq_S1-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_S1-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_S1-Freq_0)^4))*Beta_4;     
  beta_g = Beta_S1;

  Beta_P1= (1/1) *(((2*pi)^1)*((Freq_P1-Freq_0)^1))*Beta_1 ...
          +(1/2) *(((2*pi)^2)*((Freq_P1-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P1-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P1-Freq_0)^4))*Beta_4;
  beta_a = Beta_P1;

  Beta_SI= (1/1) *(((2*pi)^1)*((Freq_SI-Freq_0)^1))*Beta_1 ...
          +(1/2) *(((2*pi)^2)*((Freq_SI-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_SI-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_SI-Freq_0)^4))*Beta_4;
  beta_c = Beta_SI;

  Beta_P2= (1/1) *(((2*pi)^1)*((Freq_P2-Freq_0)^1))*Beta_1 ...
          +(1/2) *(((2*pi)^2)*((Freq_P2-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P2-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P2-Freq_0)^4))*Beta_4;
  beta_b = Beta_P2;

  Beta_S2= (1/1) *(((2*pi)^1)*((Freq_S2-Freq_0)^1))*Beta_1 ...
          +(1/2) *(((2*pi)^2)*((Freq_S2-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_S2-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_S2-Freq_0)^4))*Beta_4;
  beta_h = Beta_S2;

  Beta_P3= (1/1) *(((2*pi)^1)*((Freq_P3-Freq_0)^1))*Beta_1 ...
          +(1/2) *(((2*pi)^2)*((Freq_P3-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P3-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P3-Freq_0)^4))*Beta_4;
  beta_f = Beta_P3; 
  
else % Calculate beta vector 
  
  % Beta_0 : s^0/m, proprogayion constant
    Beta_0 =  0; % 
  % B_0    =  (1)   *(((2*pi)^0)*((f_s-f_0)^0))*Beta_0;

  % Beta_1 : s^1/m, 1/vg, inverse of group velocity
    Beta_1 =  0; % 
  % B_1    =  (1)   *(((2*pi)^1)*((f_s-f_0)^1))*Beta_1;

  % Beta_2 : s^2/m, dispersion
    Beta_2 = -(Lambda_0^2/((2*pi*c)^1))*NLM_D0;
  % Beta_2 = -(HNLF.ZDW^2/((2*pi*c)^1))*HNLF_D0;
  % B_2    =  (1/2) *(((2*pi)^2)*((f_s-f_0)^2))*Beta_2;

  % Beta_3 : s^3/m, dispersion slope
    Beta_3 =  (Lambda_0^4/((2*pi*c)^2))*((2/Lambda_0)*NLM_D0+NLM_D1);
  % Beta_3 =  (HNLF.ZDW^4/((2*pi*c)^2))*((2/HNLF.ZDW)*HNLF_D0+HNLF_D1);
  % B_3    =  (1/6) *(((2*pi)^3)*((f_s-f_0)^3))*Beta_3;

  % Beta_4 : s^4/m, slope of dispersion slope
    Beta_4 = -(Lambda_0^4/((2*pi*c)^3))*(6*NLM_D0+6*Lambda_0*NLM_D1+(Lambda_0^2)*NLM_D2);
  % Beta_4 = -(HNLF.ZDW^4/((2*pi*c)^3))*(6*HNLF_D0+6*HNLF.ZDW*HNLF_D1+(HNLF.ZDW^2)*HNLF_D2);
  % B_4    =  (1/24)*(((2*pi)^4)*((f_s-f_0)^4))*Beta_4;

  beta_vector = [Beta_0, Beta_1, Beta_2, Beta_3, Beta_4];

  %% Calculate beta for each wave 

  Beta_P4= (1/2) *(((2*pi)^2)*((Freq_P4-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P4-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P4-Freq_0)^4))*Beta_4;  
  beta_d = Beta_P4;

  Beta_S1= (1/2) *(((2*pi)^2)*((Freq_S1-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_S1-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_S1-Freq_0)^4))*Beta_4;     
  beta_g = Beta_S1;

  Beta_P1= (1/2) *(((2*pi)^2)*((Freq_P1-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P1-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P1-Freq_0)^4))*Beta_4;
  beta_a = Beta_P1;

  Beta_SI= (1/2) *(((2*pi)^2)*((Freq_SI-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_SI-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_SI-Freq_0)^4))*Beta_4;
  beta_c = Beta_SI;

  Beta_P2= (1/2) *(((2*pi)^2)*((Freq_P2-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P2-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P2-Freq_0)^4))*Beta_4;
  beta_b = Beta_P2;

  Beta_S2= (1/2) *(((2*pi)^2)*((Freq_S2-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_S2-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_S2-Freq_0)^4))*Beta_4;
  beta_h = Beta_S2;

  Beta_P3= (1/2) *(((2*pi)^2)*((Freq_P3-Freq_0)^2))*Beta_2 ...
          +(1/6) *(((2*pi)^3)*((Freq_P3-Freq_0)^3))*Beta_3 ...
          +(1/24)*(((2*pi)^4)*((Freq_P3-Freq_0)^4))*Beta_4;
  beta_f = Beta_P3; 

end;
%-------------------------------------------------------------------------------

%% Calculate delta_beta for 13 NDFWM and 9 DFWM processes

delta_beta(01) =        2*beta_b - beta_a - beta_f; % delta_beta(01) =    2b - a - f
delta_beta(02) = beta_b + beta_c - beta_a - beta_h; % delta_beta(02) = b + c - a - h
delta_beta(03) =       -2*beta_a + beta_b + beta_d; % delta_beta(03) =   -2a + b + d 
delta_beta(04) = beta_b + beta_g - beta_a - beta_c; % delta_beta(04) = b + g - a - c
delta_beta(05) =        2*beta_c - beta_a - beta_b; % delta_beta(05) =    2c - a - b
delta_beta(06) = beta_c + beta_d - beta_a - beta_g; % delta_beta(06) = c + d - a - g
delta_beta(07) =       -2*beta_a + beta_c + beta_g; % delta_beta(07) =   -2a + c + g 
delta_beta(08) = beta_c + beta_h - beta_a - beta_f; % delta_beta(08) = c + h - a - f
delta_beta(09) = beta_d + beta_f - beta_a - beta_b; % delta_beta(09) = d + f - a - b
delta_beta(10) = beta_d + beta_h - beta_a - beta_c; % delta_beta(10) = d + h - a - c
delta_beta(11) = beta_f + beta_g - beta_a - beta_h; % delta_beta(11) = f + g - a - h
delta_beta(12) =        2*beta_g - beta_a - beta_d; % delta_beta(12) =    2g - a - d
delta_beta(13) = beta_g + beta_h - beta_a - beta_b; % delta_beta(13) = g + h - a - b
delta_beta(14) =        2*beta_h - beta_b - beta_f; % delta_beta(14) =    2h - b - f
delta_beta(15) = beta_c + beta_g - beta_b - beta_d; % delta_beta(15) = c + g - b - d
delta_beta(16) =       -2*beta_b + beta_c + beta_h; % delta_beta(16) =   -2b + c + h 
delta_beta(17) = beta_d + beta_h - beta_b - beta_g; % delta_beta(17) = d + h - b - g
delta_beta(18) = beta_f + beta_g - beta_b - beta_c; % delta_beta(18) = f + g - b - c
delta_beta(19) = beta_b + beta_h - beta_c - beta_f; % delta_beta(19) = b + h - c - f
delta_beta(20) =       -2*beta_c + beta_d + beta_f; % delta_beta(20) =   -2c + d + f
delta_beta(21) =       -2*beta_c + beta_g + beta_h; % delta_beta(21) =   -2c + g + h
delta_beta(22) = beta_g + beta_h - beta_d - beta_f; % delta_beta(22) = g + h - d - f
%-------------------------------------------------------------------------------

%% Main switch
switch( SimMod.Wave_num );
  case( 03 ) % 3-wave case, delta_beta among three entering wave
    delta_beta(01) = 0;
    delta_beta(02) = 0;
    delta_beta(03) = 0;
    delta_beta(04) = 0;
  % delta_beta(05) = delta_beta(05); % Prime 3-wave DG-FWM
    delta_beta(05) = 2*(1/ 2)*((2*pi)^2)*((Freq_SI-Freq_c)^2-(Freq_d)^2)*Beta_2 ...
                   + 2*(1/24)*((2*pi)^4)*((Freq_SI-Freq_c)^4-(Freq_d)^4)*Beta_4;
    delta_beta(06) = 0;
    delta_beta(07) = 0;
    delta_beta(08) = 0;
    delta_beta(09) = 0;
    delta_beta(10) = 0;
    delta_beta(11) = 0;
    delta_beta(12) = 0;
    delta_beta(13) = 0;
    delta_beta(14) = 0;
    delta_beta(15) = 0;
    delta_beta(16) = 0;
    delta_beta(17) = 0;
    delta_beta(18) = 0;
    delta_beta(19) = 0;
    delta_beta(20) = 0;
    delta_beta(21) = 0;
    delta_beta(22) = 0;
  %-----------------------------------------------------------------------------
  case( 05 ) % 5-wave case, 
  % delta_beta(01) = 0;
  % delta_beta(02) = 0;
  % delta_beta(03) = 0;
  % delta_beta(04) = 0; %
  % delta_beta(05) = 0; % Prime 3-wave DG-FWM
  % delta_beta(06) = 0; % 
  % delta_beta(07) = 0; %
  % delta_beta(08) = 0; % 
  % delta_beta(09) = 0; %
  % delta_beta(10) = 0;
  % delta_beta(11) = 0;
  % delta_beta(12) = 0;
  % delta_beta(13) = 0; %
  % delta_beta(14) = 0;
  % delta_beta(15) = 0;
  % delta_beta(16) = 0; %
  % delta_beta(17) = 0;
  % delta_beta(18) = 0;
  % delta_beta(19) = 0;
  % delta_beta(20) = 0;
  % delta_beta(21) = 0; %
  % delta_beta(22) = 0;
  %-----------------------------------------------------------------------------
  case( 07 ) % 7-wave case, all the phase-mismatch should be taken into account
  % delta_beta(01) = 0;
  % delta_beta(02) = 0;
  % delta_beta(03) = 0;
  % delta_beta(04) = 0;
  % delta_beta(05) = 0;
  % delta_beta(06) = 0;
  % delta_beta(07) = 0;
  % delta_beta(08) = 0;
  % delta_beta(09) = 0;
  % delta_beta(10) = 0;
  % delta_beta(11) = 0;
  % delta_beta(12) = 0;
  % delta_beta(13) = 0;
  % delta_beta(14) = 0;
  % delta_beta(15) = 0;
  % delta_beta(16) = 0;
  % delta_beta(17) = 0;
  % delta_beta(18) = 0;
  % delta_beta(19) = 0;
  % delta_beta(20) = 0;
  % delta_beta(21) = 0;
  % delta_beta(22) = 0;
  %-----------------------------------------------------------------------------
  otherwise 
  %-----------------------------------------------------------------------------  
end % end of switch
%-------------------------------------------------------------------------------

end % end of function f_delta_beta
%-------------------------------------------------------------------------------

%%
%-------------------------------------------------------------------------------
% end of file f_delta_beta.m
%-------------------------------------------------------------------------------