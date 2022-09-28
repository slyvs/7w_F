%-------------------------------------------------------------------------------
% Beginning of file f_NLM.m
%-------------------------------------------------------------------------------
% 
% Project : Calculation of nonlinear propagation in fiber for FOPA based PSA 
%           using 7-wave model
% FileName: f_NLM.m
% Function: Parameters of Nonlinear Medium (NLM)
% Version : 
%   v0.0 @2015.03.30 created  by W.L. Xie 
%        Add HNLF of LAC-GLOP
%   v0.1 @2015.09.15 released by W.L. Xie
%        * Change HNLF parameters into classes HNLF.*
%   v0.2 @2016.07.25 released by W.L. Xie, on the flight to Paris
%        * Change the function and class from HNLF to NLM
% Describe: Parameter of High Non-linear Fiber (HNLF)
%           04: 01+02+03  OFS HNLF standard @ ZDW=1547.50nm
%           07: 06+05     OFS HNLF Spine @ ZDW=1566.00nm
%           10: only for simulation test
%           11: 2014.ECOC Maxime Baillot, 6-wave
%           12: 2009.CLEO Mingyi Gao, 9-wave
%           13: 2012.OFC  Mingyi Gao, 7-wave
%           14: 2012.OL   Mingyi Gao, 7-wave
%           15: 2008.OE   Bismuth-Oxide-based HNLF
%           16: 2011.OL   C.Lun, 3-w analytical
%           17: 2012.OE   C.Lun, 3-w numerical
%           18: 2014,JLT  D.Liu, 7-w numerical
%           19: 2015.     Photonic Crystal Waveguide from A.F's PRA paper
%           20: 2016.     GaP SL PCW, 
%-------------------------------------------------------------------------------

function f_NLM( n )

%% Global constant 
global c;

%% Global varibles 
global NLM;          % HNLF parameters
%      NLM.Len;      % m, length of HNLF
%      NLM.ZDW;      % m, zero dispersion wavelength
%      NLM.ZDF;      % m, zero dispersion frequency
%      NLM.D0;       % ps/nm^1/km, HNLF dispersion
%      NLM.D1;       % ps/nm^2/km, HNLF dispersion slope
%      NLM.D2;       % ps/nm^3/km, slope of HNLF dispersion slope
%      NLM.Loss;     % dB, total loss @1550
%      NLM.PMD;      % ps, total PMD C-band
%      NLM.Gamma;    % 1/(W*km), Gamma, Nonlinear coefficient
%      NLM.EffAr;    % Effective area
%      NLM.Alpha_dB; % dB/km, Attenuation of HNLF   
%      NLM.Alpha;    % m^-1, linear attenuation coefficience of HNLF

%% NLM parameters 
switch( n );
  case( 00 ) % Simulation ONLY
    %% HNLF ID: Simulation
    NLM.Len   = 100;          % m, length of HNLF
    NLM.ZDW   = 1565.00e-9;   % m, zero dispersion wavelength,default:1547.00e-9
  % HNLF.D1    = 0.017;        % ps/km/nm^2, 
    NLM.D0    = 0;
    NLM.D1    = 0.027e3;      %  s/    m^3, Default:0.017e3
    NLM.D2    = 0.000;        % ps/nm^3/km, slope of HNLF dispersion slope,       
  % HNLF.Gamma = 11.3 ;        % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 15.0e-3;      % 1/(W* m), Gamma, Nonlinear coefficient,default:11.3e-3
    NLM.alpha_dB = 0.0;                % dB/km, HNLF_Loss/(HNLF_Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10)); 
  case( 01 ) % 
    %% HNLF ID: 131233340001
    NLM.Len   = 204;           % m, length of HNLF
    NLM.ZDW   = 1548.00e-9;    % m, zero dispersion wavelength
    NLM.D0_15 = 0.04;          % ps/nm/km, HNLF dispersion @ 1550nm
    NLM.D1_15 = 0.017;         % ps/nm^2/km, HNLF dispersion slope @ 1550nm
    NLM.D0    = 0.017e3*(Lambda_0 - NLM.ZDW); % s/m^1, 
  % HNLF.D1    = 0.017;         % ps/m/nm^2, 
    NLM.D1    = 0.017e3;       %  s/   m^2, 
    NLM.D2    = 0.000;         % ps/km/nm^3, slope of HNLF dispersion slope
  % HNLF.D2    =-0.003e-3;      % ps/nm^3/m, slope of HNLF dispersion slope, 
                                % calculated using D0 = ax^2 + bx + c
                                % a = -0.0015, b = 0.023, c = 0
    NLM.Loss  = 0.77;          % dB, total loss @1550
    NLM.PMD   = 0.04;          % ps, total PMD C-band
    NLM.Gamma = 11.3;          % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.alpha_dB = 0.900;              % dB/km, HNLF_Loss/(HNLF_Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
  case( 02 ) % 
    %% HNLF ID: 131233340002
    NLM.Len   = 300;           % m, length of HNLF
    NLM.ZDW   = 1548.00e-9;    % m, zero dispersion wavelength
    NLM.D0    = 0.03;          % ps/km/nm  , HNLF dispersion
    NLM.D1    = 0.017;         % ps/km/nm^2, HNLF dispersion slope
    NLM.D2    = 0.002;         % ps/km/nm^3, slope of HNLF dispersion slope, 
                                % calculated using D0 = ax^2 + bx + c
                                % a = 0.001, b = 0.013, c = 0
    NLM.Loss  = 0.75;          % dB, total loss @1550
    NLM.PMD   = 0.04;          % ps, total PMD C-band
    NLM.Gamma = 11.3;          % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.alpha_dB = 0.900;              % dB/km, HNLF_Loss/(HNLF_Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
  case( 03 ) % 
    %% HNLF ID: 131233340003
    NLM.Len   = 507;           % m, length of HNLF
    NLM.ZDW  = 1547.00e-9;    % m, zero dispersion wavelength
  % HNLF.D0    = 0.06;          % ps/km/nm  , HNLF dispersion
    NLM.D0    = 0.04e-6;       % ps/    m^2, HNLF dispersion
  % HNLF.D1    = 0.017;         % ps/km/nm^2, HNLF dispersion slope
    NLM.D1    = 0.017e3;       %  s/    m^3, HNLF dispersion slope
  % HNLF.D2    = 0.000;         % ps/km/nm^3, slope of HNLF dispersion slope
    NLM.D2    = 0.000e12;      % ps/    m^4, slope of HNLF dispersion slope
    NLM.Loss  = 1.05;          % dB, total loss @1550
    NLM.PMD   = 0.04;          % ps, total PMD C-band
  % HNLF.Gamma = 11.3;          % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 11.3e-3;       % 1/(W* m), Gamma, Nonlinear coefficient     
    NLM.alpha_dB = 0.900;              % dB/km, HNLF_Loss/(HNLF_Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
  case( 04 ) % 
    %% HNLF: (01)131233340001 + (02)131233340002 + (03)131233340003
    NLM.Len   = 204+300+507;   % m, length of HNLF
    NLM.ZDW   = 1547.50e-9;    % m, zero dispersion wavelength
    NLM.ZDF   = c/NLM.ZDW;    % Hz, zero dispersion frequency
  % NLM.D0    = 0.04;          % ps/km/nm^1, HNLF dispersion @1550nm
    NLM.D0    = 0.04e-6;       %  s/    m^2, HNLF dispersion @1550nm
  % NLM.D1    = 0.017;         % ps/km/nm^2, HNLF dispersion slope
    NLM.D1    = 0.017e3;       %  s/    m^3, HNLF dispersion slope
  % NLM.D2    = 0.005;         % ps/km/nm^3, slope of HNLF dispersion slope
    NLM.D2    = 0.000e12;      %  s/    m^4, slope of HNLF dispersion slope
    NLM.Loss  = 0.77+0.75+1.05;% dB, total loss @1550
    NLM.PMD   = 0.04;          % ps, total PMD C-band
  % NLM.Gamma = 11.3 ;         % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 11.3e-3;       % 1/(W* m), Gamma, Nonlinear coefficient     
    NLM.alpha_dB = 0.900;              % dB/km, HNLF.Loss/(HNLF.Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
  case( 05 ) % 
    %% HNLF ID: 14251113090001
    NLM.Len   = 255;           % m, length of HNLF
    NLM.ZDW   = 1566.00e-9;    % m, zero dispersion wavelength
    NLM.D0    = -1.30;         % ps/nm/km, HNLF dispersion
    NLM.D1    = 0.083;         % ps/nm^2/km, HNLF dispersion slope
    NLM.Loss  = 0.64;          % dB, total loss @1550
    NLM.PMD   = 0.05;          % ps, total PMD C-band @1550nm
    NLM.Gamma = 8.7;           % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.EffAr = 16.1;          % um^2, effective area @1550nm
    NLM.alpha_dB = 0.960;      % dB/km, HNLF_Loss/(HNLF_Len/1e3);
  case( 06 ) % 
    %% HNLF ID: 14251113090002
    NLM.Len   = 256;           % m, length of HNLF
    NLM.ZDW   = 1566.00e-9;    % m, zero dispersion wavelength
    NLM.D0    = -1.30;         % ps/nm/km, HNLF dispersion
    NLM.D1    = 0.083;         % ps/nm^2/km, HNLF dispersion slope
    NLM.Loss  = 0.61;          % dB, total loss @1550
    NLM.PMD   = 0.05;          % ps, total PMD C-band @1550nm
    NLM.Gamma = 8.7;           % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.EffAr = 16.1;          % um^2, effective area @1550nm
    NLM.alpha_dB   = 0.960;    % dB/km, HNLF_Loss/(HNLF_Len/1e3);
  case( 07 ) % 
    %% HNLF: (05)14251113090001 + (06)14251113090002
    NLM.Len   = 255+256;       % m, length of HNLF
    NLM.ZDW   = 1566.00e-9;    % m, zero dispersion wavelength
  % HNLF.D0    =-1.30;          % ps/km/nm^1, HNLF dispersion @1550nm
    NLM.D0    =-1.30e-6;       %  s/    m^2, HNLF dispersion @1550nm
  % HNLF.D1    = 0.083;         % ps/km/nm^2, 
    NLM.D1    = 0.083e3;       %  s/    m^3, HNLF dispersion slope
  % HNLF.D2    = 0.000;         % ps/km/nm^3, slope of HNLF dispersion slope 
    NLM.D2    = 0.000e12;      %  s/    m^4, slope of HNLF dispersion slope
    NLM.Loss  = 0.64+0.61;     % dB, total loss @1550
    NLM.PMD   = 0.05;          % ps, total PMD C-band 
  % HNLF.Gamma = 8.7 ;          % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 8.7e-3;        % 1/(W* m), Gamma, Nonlinear coefficient     
    NLM.alpha_dB = 0.960;              % dB/km, HNLF.Loss/(HNLF.Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));       
  case( 09 ) % 
    %% HNLF ID: XXXXXXXXXXXXXX
    NLM.ZDW   = 1565.00e-9;    % m, zero dispersion wavelength
  % HNLF.D0    = 2.0e-6;        % ps/(nm  *m), dispersion
    NLM.D0    = DD*(Lambda_0 - NLM.ZDW); % ps/(nm*km), dispersion calculated by slope
  % HNLF.D1    = 0.027;         % ps/(nm^2*m), dispersion slope 1st derivative
    NLM.D1    = 0.027e3;       %  s/      m^3, dispersion slope 1st derivative
    NLM.D2    = 0.000e0;       % ps/(nm^3*m), slope of dispersion slope, second derivative 
  case( 10 ) % 
    %% HNLF ID: XXXXXXXXXXXXXX
    NLM.Len   = 100;
    NLM.ZDW   = 1550.00e-9; % m, zero dispersion wavelength
  % HNLF.D0    = 0.017e3*(Lambda_0 - HNLF.ZDW); % s/m^1, 
  % HNLF.D0    = 0.04;          % ps/nm/km^1, HNLF dispersion
    NLM.D0    = 0.04e-6;       % ps/    m^2, HNLF dispersion
  % HNLF.D1    = 0.017;         % ps/km/nm^2, 
    NLM.D1    = 0.017e3;       %  s/    m^3, HNLF dispersion slope
  % HNLF.D2    = 0.005;         % ps/km/nm^3, slope of HNLF dispersion slope
    NLM.D2    = 0.000e12;      % ps/    m^4, slope of HNLF dispersion slope
    NLM.Loss  = 0.10;          % dB, total loss @1550
    NLM.PMD   = 0.05;          % ps, total PMD C-band 
  % HNLF.Gamma = 10.0 ;         % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 10.0e-3;       % 1/(W* m), Gamma, Nonlinear coefficient     
    NLM.alpha_dB = 0.000;              % dB/km, HNLF.Loss/(HNLF.Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));    
  case( 11 ) % 
    %% 2014.ECOC Maxime Baillot, 6-wave
    NLM.Len   = 10.40e3;       % m, length of HNLF
    NLM.ZDW   = 1548.00e-9;    % m, zero dispersion wavelength
  % HNLF.D1    = 0.077;         % ps/(km*nm^2), HNLF dispersion slope
    NLM.D1    = 0.077e3;       %  s/      m^3, HNLF dispersion slope in International Unit
    NLM.D2    = 0.000;         % ps/(km*nm^3), slope of HNLF dispersion slope
    NLM.Gamma = 3.30e-3;       % 1/(W*m), Gamma, Nonlinear coefficient
    NLM.alpha_dB = 0.200;              % 0.200 dB/km, HNLF_Loss/(HNLF_Len/1e3);      
    NLM.alpha    = NLM.alpha_dB/4343; %  losses in m-1, 10^(-Alpha_dB/(10*1e3));
  case( 12 ) % 
    %% 2009.CLEO Mingyi Gao, 9-wave 
    NLM.Len   = 1000;          % m, length of HNLF
    NLM.ZDW   = 1565.00e-9;    % m, zero dispersion wavelength
    NLM.ZDF   = c/NLM.ZDW;    % Hz, zero dispersion frequency
  % HNLF.D0    = 0.000;         % ps/km/nm^1, HNLF dispersion @1550nm
    NLM.D0    = 0.000e-6;      %  s/    m^2, HNLF dispersion @1550nm
  % HNLF.D1    = 0.027;         % ps/km/nm^2, HNLF dispersion slope
    NLM.D1    = 0.027e3;       %  s/    m^3, HNLF dispersion slope
  % HNLF.D2    = 0.000;         % ps/km/nm^3, slope of HNLF dispersion slope
    NLM.D2    = 0.000e12;      %  s/    m^4, slope of HNLF dispersion slope
  % HNLF.Gamma = 15.0 ;         % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 15.0e-3;       % 1/(W* m), Gamma, Nonlinear coefficient
    NLM.alpha_dB = 0.001;              % dB/km, HNLF_Loss/(HNLF_Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
  case( 13 ) % 
    %% 2012.OFC Mingyi Gao, 7-wave 
    NLM.Len   = 600;           % m, length of HNLF
    NLM.ZDW   = 1541.67e-9;    % m, zero dispersion wavelength
    NLM.ZDF   = c/NLM.ZDW;    % Hz, zero dispersion frequency
  % HNLF.D0    = 0.494;         % ps/km/nm^1, HNLF dispersion @1550nm
    NLM.D0    = 0.494e-6;      % ps/    m^2, HNLF dispersion @1550nm
  % HNLF.D1    = 0.025;         % ps/km/nm^2, HNLF dispersion slope
    NLM.D1    = 0.025e3;       %  s/    m^3, HNLF dispersion slope
    NLM.D2    = 0.000;         % ps/km/nm^3, slope of HNLF dispersion slope, 
  % HNLF.Gamma = 15.0 ;         % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 10.0e-3;       % 1/(W* m), Gamma, Nonlinear coefficient
    NLM.alpha_dB = 0.670;              % dB/km, HNLF_Loss/(HNLF_Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
  case( 14 ) % 
    %% 2012.OL Mingyi Gao, 7-wave 
  % Lambda_S   = 1564.00e-9;    % m, signal wavelength 
    NLM.Len   = 600;           % m, length of HNLF
    NLM.ZDW   = 1541.67e-9;    % m, zero dispersion wavelength
  % HNLF.D0    = 0.482;         % ps/km/nm^1,  
    NLM.D0    = 0.482e-6;      %  s/    m^2,  
  % HNLF.D0    = 0;             %  s/    m^2,  
  % HNLF.D1    = 0.025;         % ps/km/nm^2, 
    NLM.D1    = 0.025e3;       %  s/    m^3, 
    NLM.D2    = 0.000;         % ps/km/nm^3, slope of HNLF dispersion slope, 
  % HNLF.Gamma = 10.0 ;         % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 10.0e-3;       % 1/(W* m), Gamma, Nonlinear coefficient
    NLM.alpha_dB = 0.670;              % dB/km, HNLF_Loss/(HNLF_Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
  case( 19 ) % 
    %% Photonic Crystal Waveguide Afraido
  % HNLF.Len   = 1.5e-3;        % m, length of HNLF
    NLM.Len   = 1.0e-3;        % m, length of HNLF
    NLM.ZDW   = 1552.30e-9;    % m, zero dispersion wavelength
    NLM.ZDF   = c/NLM.ZDW;    % Hz, zero dispersion frequency
  % HNLF.D0    = 0.04;          % ps/nm/km^1, HNLF dispersion
    NLM.D0    = 0.000e-6;      % ps/    m^2, HNLF dispersion
  % HNLF.D1    = 0.017;         % ps/km/nm^2, 
    NLM.D1    = 0.000e3;       %  s/    m^3, HNLF dispersion slope
  % HNLF.D2    = 0.005;         % ps/km/nm^3, slope of HNLF dispersion slope
    NLM.D2    = 0.000e12;      % ps/    m^4, slope of HNLF dispersion slope
  % HNLF.Loss  = 0.77+0.75+1.05;% dB, total loss @1550
  % HNLF.PMD   = 0.04;          % ps, total PMD C-band
  % HNLF.Gamma = 2000e3;        % 1/(W*km), Gamma, Nonlinear coefficient
    NLM.Gamma = 2000.0;        % 1/(W* m), Gamma, Nonlinear coefficient     
    NLM.alpha_dB = 30.0e5;          % dB/km, HNLF.Loss/(HNLF.Len/1e3);
    NLM.alpha_dB = 00.0e5;          % dB/km, HNLF.Loss/(HNLF.Len/1e3);
    NLM.alpha    = NLM.alpha_dB/4343; % (-1/50)*log(10^-(Alpha_dB/10));
end; % end of switch f_HNLF(n)
  
end % end of function f_NLM
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% end of file f_NLM.m
%-------------------------------------------------------------------------------