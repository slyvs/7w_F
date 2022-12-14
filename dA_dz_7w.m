%-------------------------------------------------------------------------------
% Beginning of file dA_dz_7w.m 
%-------------------------------------------------------------------------------
%
% Optical Parametric Amplification: 7-wave model, analysis in numerical fashion 
% 
%                P1        P2
%                |         |
%      P3        |         |         P4
%      |         |   S&I   |         |
%      |    S1   |    |    |    S2   |
%      |    |    |    |    |    |    |
%      |    |    |    |    |    |    |
% ------------------------------------------> wavelength(nm)
%      A4   A6   A1   A3   A2   A7   A5    
%      d    g    a    c    b    h    f
%     (A5) (A3) (A1) (A0) (A2) (A4) (A6) 
% 
% Project : Calculation of nonlinear propagation in fiber for OPA based PSA 
%           using 7-wave model
% FileName: dA_dz_7w.m
% Function: Numerically calculate the nonlinear propagation of OPA based PIA/PSA 
%           in nonlinear medium using coupled wave ordinary differential equation  
%           of 7-wave model, including:
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
% Differential equation for field evolution of  7 waves
% Differential equation for relative phases of 22 FWM processes
%-------------------------------------------------------------------------------
%% Fun: 
function dA = dA_dz_7wave(z,A)
% z            : Length of medium
%  A(01)~ A(07): Power evolution of 7 waves 
% dA(01)~dA(07): Differential for power of 7 waves 
%  A(08)~ A(29): Relative phase and Phase-matching of 13 NDFWN and 9 DFWM 
% dA(08)~dA(29): Differential of relative phases and phase-matching

%% Simulation paramters 
global SimMod;   % 
global ExcSet;   % 
global Exc;

%% Global varibles 
global NLM;
global delta_beta;
global dPhi;
global delta_theta; % rad, Relative phase of 13 degenerate and 9 non-degenerate FWM

%% 
% delta_theta(01) = 2*E( 9) - E( 8) -   E(12)        ; %     degenerate FWM
% delta_theta(02) =   E( 9) + E(10) -   E( 8) - E(14); % non-degenerate FWM
% delta_theta(03) =   E( 9) + E(11) - 2*E( 8);         %     degenerate FWM
% delta_theta(04) =   E( 9) + E(13) -   E( 8) - E(10); % non-degenerate FWM
% delta_theta(05) = 2*E(10) - E( 8) -   E( 9)        ; %     degenerate FWM
% delta_theta(06) =   E(10) + E(11) -   E( 8) - E(13); % non-degenerate FWM
% delta_theta(07) =   E(10) + E(13) - 2*E( 8)        ; %     degenerate FWM
% delta_theta(08) =   E(10) + E(14) -   E( 8) - E(12); % non-degenerate FWM
% delta_theta(09) =   E(11) + E(12) -   E( 8) - E( 9); % non-degenerate FWM
% delta_theta(10) =   E(11) + E(14) -   E( 8) - E(10); % non-degenerate FWM
% delta_theta(11) =   E(12) + E(13) -   E( 8) - E(14); % non-degenerate FWM
% delta_theta(12) = 2*E(13) - E( 8) -   E(11)        ; %     degenerate FWM
% delta_theta(13) =   E(13) + E(14) -   E( 8) - E( 9); % non-degenerate FWM
% delta_theta(14) = 2*E(14) - E( 9) -   E(12)        ; %     degenerate FWM
% delta_theta(15) =   E(10) + E(13) -   E( 9) - E(11); % non-degenerate FWM
% delta_theta(16) =   E(10) + E(14) - 2*E( 9)        ; %     degenerate FWM
% delta_theta(17) =   E(11) + E(14) -   E( 9) - E(13); % non-degenerate FWM
% delta_theta(18) =   E(12) + E(13) -   E( 9) - E(10); % non-degenerate FWM
% delta_theta(19) =   E( 9) + E(14) -   E(10) - E(12); % non-degenerate FWM
% delta_theta(20) =   E(11) + E(12) - 2*E(10)        ; %     degenerate FWM
% delta_theta(21) =   E(13) + E(14) - 2*E(10)        ; %     degenerate FWM
% delta_theta(22) =   E(13) + E(14) -   E(11) - E(12); % non-degenerate FWM

%% -----------------------------------------------------------------------------
if SimMod.FLD_nPWR == 1, % Calculate propagation using Field

  %% CODE for fields evolution of all the interacting waves
  delta_theta       = zeros(22,1); % rad, Initial the relative phase of 22 FWMs

  % Initial coupled ordinary differential equation for complex field of 7 waves
  dA    = zeros(7,1); % 

  % Coupled ordinary differential equation of Pump1 (A1) field
  dA(01)= - (NLM.alpha/2)*A(01) ...
          + 1i*NLM.Gamma*( ...
            A(01)*(abs(A(01))^2+2*(abs(A(02))^2+abs(A(03))^2 ...
                                  +abs(A(04)*ExcSet.HP)^2+abs(A(05)*ExcSet.HP)^2 ...
                                  +abs(A(06)*ExcSet.HS)^2+abs(A(07)*ExcSet.HS)^2))...
        +   A(02)*A(02)*conj(A(05))*exp( (delta_theta(01)+delta_beta(01)*z)*1i )*ExcSet.HP...
        +   A(03)*A(03)*conj(A(02))*exp( (delta_theta(05)+delta_beta(05)*z)*1i )...
        +   A(06)*A(06)*conj(A(04))*exp( (delta_theta(12)+delta_beta(12)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(03)*conj(A(07))*exp( (delta_theta(02)+delta_beta(02)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(04)*conj(A(01))*exp( (delta_theta(03)+delta_beta(03)*z)*1i )*ExcSet.HP...
        + 2*A(04)*A(05)*conj(A(02))*exp( (delta_theta(09)+delta_beta(09)*z)*1i )*ExcSet.HP...
        + 2*A(03)*A(04)*conj(A(06))*exp( (delta_theta(06)+delta_beta(06)*z)*1i )*ExcSet.HS...
        + 2*A(04)*A(07)*conj(A(03))*exp( (delta_theta(10)+delta_beta(10)*z)*1i )*ExcSet.HS...
        + 2*A(03)*A(07)*conj(A(05))*exp( (delta_theta(08)+delta_beta(08)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(06)*conj(A(03))*exp( (delta_theta(04)+delta_beta(04)*z)*1i )*ExcSet.HS...
        + 2*A(03)*A(06)*conj(A(01))*exp( (delta_theta(07)+delta_beta(07)*z)*1i )*ExcSet.HS...
        + 2*A(05)*A(06)*conj(A(07))*exp( (delta_theta(11)+delta_beta(11)*z)*1i )*ExcSet.HS...
        + 2*A(06)*A(07)*conj(A(02))*exp( (delta_theta(13)+delta_beta(13)*z)*1i )*ExcSet.HS...
                       );

  % Coupled ordinary differential equation of Pump2 (A2) field
  dA(02)= - (NLM.alpha/2)*A(02) ...
          + 1i*NLM.Gamma*( ...
            A(02)*(abs(A(02))^2+2*(abs(A(01)          )^2+abs(A(03)          )^2 ...
                                  +abs(A(04)*ExcSet.HP)^2+abs(A(05)*ExcSet.HP)^2 ...
                                  +abs(A(06)*ExcSet.HS)^2+abs(A(07)*ExcSet.HS)^2))...
        +   A(01)*A(01)*conj(A(04))*exp(-(delta_theta(03)+delta_beta(03)*z)*1i )*ExcSet.HP...
        +   A(03)*A(03)*conj(A(01))*exp( (delta_theta(05)+delta_beta(05)*z)*1i )...
        +   A(07)*A(07)*conj(A(05))*exp( (delta_theta(14)+delta_beta(14)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(03)*conj(A(06))*exp(-(delta_theta(04)+delta_beta(04)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(05)*conj(A(02))*exp(-(delta_theta(01)+delta_beta(01)*z)*1i )*ExcSet.HP...
        + 2*A(01)*A(07)*conj(A(03))*exp(-(delta_theta(02)+delta_beta(02)*z)*1i )*ExcSet.HS...
        + 2*A(03)*A(05)*conj(A(07))*exp(-(delta_theta(19)+delta_beta(19)*z)*1i )*ExcSet.HS...
        + 2*A(03)*A(06)*conj(A(04))*exp( (delta_theta(15)+delta_beta(15)*z)*1i )*ExcSet.HS...
        + 2*A(03)*A(07)*conj(A(02))*exp( (delta_theta(16)+delta_beta(16)*z)*1i )*ExcSet.HS...
        + 2*A(05)*A(06)*conj(A(03))*exp( (delta_theta(18)+delta_beta(18)*z)*1i )*ExcSet.HS...
        + 2*A(04)*A(05)*conj(A(01))*exp( (delta_theta(09)+delta_beta(09)*z)*1i )*ExcSet.HP...
        + 2*A(04)*A(07)*conj(A(06))*exp( (delta_theta(17)+delta_beta(17)*z)*1i )*ExcSet.HS...
        + 2*A(06)*A(07)*conj(A(01))*exp( (delta_theta(13)+delta_beta(13)*z)*1i )*ExcSet.HS...
                       );

  % Coupled ordinary differential equation of Signal&ddIlder (A3) field
  dA(03)= - (NLM.alpha/2)*A(03) ...
          + 1i*NLM.Gamma*( ...
            A(03)*(abs(A(03))^2+2*(abs(A(01)          )^2+abs(A(02)          )^2 ...
                                  +abs(A(04)*ExcSet.HP)^2+abs(A(05)*ExcSet.HP)^2 ...
                                  +abs(A(06)*ExcSet.HS)^2+abs(A(07)*ExcSet.HS)^2))...
        +   A(01)*A(01)*conj(A(06))*exp(-(delta_theta(07)+delta_beta(07)*z)*1i )*ExcSet.HS... 
        +   A(02)*A(02)*conj(A(07))*exp(-(delta_theta(16)+delta_beta(16)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(02)*conj(A(03))*exp(-(delta_theta(05)+delta_beta(05)*z)*1i )... 
        + 2*A(01)*A(05)*conj(A(07))*exp(-(delta_theta(08)+delta_beta(08)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(06)*conj(A(04))*exp(-(delta_theta(06)+delta_beta(06)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(07)*conj(A(02))*exp(-(delta_theta(02)+delta_beta(02)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(04)*conj(A(06))*exp(-(delta_theta(15)+delta_beta(15)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(06)*conj(A(01))*exp( (delta_theta(04)+delta_beta(04)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(07)*conj(A(05))*exp( (delta_theta(19)+delta_beta(19)*z)*1i )*ExcSet.HS...
        + 2*A(04)*A(07)*conj(A(01))*exp( (delta_theta(10)+delta_beta(10)*z)*1i )*ExcSet.HS...
        + 2*A(04)*A(05)*conj(A(03))*exp( (delta_theta(20)+delta_beta(20)*z)*1i )*ExcSet.HP...
        + 2*A(06)*A(07)*conj(A(03))*exp( (delta_theta(21)+delta_beta(21)*z)*1i )*ExcSet.HS...
        + 2*A(05)*A(06)*conj(A(02))*exp( (delta_theta(18)+delta_beta(18)*z)*1i )*ExcSet.HS...
                       );   

  % Coupled ordinary differential equation of Pump3 (A4) field
  dA(04)= - (NLM.alpha/2)*A(04)*ExcSet.HP ...
          + 1i*NLM.Gamma*( ...
            A(04)*(abs(A(04))^2+2*(abs(A(01))^2+abs(A(02))^2 ...
                                  +abs(A(03)          )^2+abs(A(05)*ExcSet.HP)^2 ...
                                  +abs(A(06)*ExcSet.HS)^2+abs(A(07)*ExcSet.HS)^2))...
        +   A(01)*A(01)*conj(A(02))*exp(-(delta_theta(03)+delta_beta(03)*z)*1i )*ExcSet.HP...
        +   A(03)*A(03)*conj(A(05))*exp(-(delta_theta(20)+delta_beta(20)*z)*1i )*ExcSet.HP...
        +   A(06)*A(06)*conj(A(01))*exp( (delta_theta(12)+delta_beta(12)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(02)*conj(A(05))*exp(-(delta_theta(09)+delta_beta(09)*z)*1i )*ExcSet.HP...
        + 2*A(01)*A(03)*conj(A(07))*exp(-(delta_theta(10)+delta_beta(10)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(06)*conj(A(03))*exp(-(delta_theta(06)+delta_beta(06)*z)*1i )*ExcSet.HS...
        + 2*A(03)*A(06)*conj(A(02))*exp( (delta_theta(15)+delta_beta(15)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(06)*conj(A(07))*exp(-(delta_theta(17)+delta_beta(17)*z)*1i )*ExcSet.HS...
        + 2*A(06)*A(07)*conj(A(05))*exp( (delta_theta(22)+delta_beta(22)*z)*1i )*ExcSet.HS...
                       )*ExcSet.HP;

  % Coupled ordinary differential equation of Pump4 (A5) field
  dA(05)= - (NLM.alpha/2)*A(05)*ExcSet.HP ...
          + 1i*NLM.Gamma*( ...
            A(05)*(abs(A(05))^2+2*(abs(A(01)          )^2+abs(A(02)          )^2 ...
                                  +abs(A(03)          )^2+abs(A(04)*ExcSet.HP)^2 ...
                                  +abs(A(06)*ExcSet.HS)^2+abs(A(07)*ExcSet.HS)^2))...
        +   A(02)*A(02)*conj(A(01))*exp( (delta_theta(01)+delta_beta(01)*z)*1i )*ExcSet.HP...
        +   A(03)*A(03)*conj(A(04))*exp(-(delta_theta(20)+delta_beta(20)*z)*1i )*ExcSet.HP...
        +   A(07)*A(07)*conj(A(02))*exp( (delta_theta(14)+delta_beta(14)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(02)*conj(A(04))*exp(-(delta_theta(09)+delta_beta(09)*z)*1i )*ExcSet.HP...
        + 2*A(01)*A(07)*conj(A(06))*exp(-(delta_theta(11)+delta_beta(11)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(03)*conj(A(06))*exp(-(delta_theta(18)+delta_beta(18)*z)*1i )*ExcSet.HS...
        + 2*A(06)*A(07)*conj(A(04))*exp( (delta_theta(22)+delta_beta(22)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(07)*conj(A(03))*exp( (delta_theta(19)+delta_beta(19)*z)*1i )*ExcSet.HS...
        + 2*A(03)*A(07)*conj(A(01))*exp( (delta_theta(08)+delta_beta(08)*z)*1i )*ExcSet.HS...
                       )*ExcSet.HP; 

  % Coupled ordinary differential equation of Signal1 (A6) field 
  dA(06)= - (NLM.alpha/2)*A(06)*ExcSet.HS ...
          + 1i*NLM.Gamma*( ...
            A(06)*(abs(A(06))^2+2*(abs(A(01)          )^2+abs(A(02)          )^2 ...
                                  +abs(A(03)*ExcSet.HP)^2+abs(A(04)*ExcSet.HP)^2 ...
                                  +abs(A(05)*ExcSet.HS)^2+abs(A(07)*ExcSet.HS)^2))...
        +   A(01)*A(01)*conj(A(03))*exp(-(delta_theta(07)+delta_beta(07)*z)*1i )*ExcSet.HS... 
        +   A(03)*A(03)*conj(A(07))*exp(-(delta_theta(21)+delta_beta(21)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(02)*conj(A(07))*exp(-(delta_theta(13)+delta_beta(13)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(03)*conj(A(02))*exp(-(delta_theta(04)+delta_beta(04)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(04)*conj(A(06))*exp(-(delta_theta(12)+delta_beta(12)*z)*1i )*ExcSet.HP...
        + 2*A(01)*A(07)*conj(A(05))*exp(-(delta_theta(11)+delta_beta(11)*z)*1i )*ExcSet.HP...
        + 2*A(04)*A(05)*conj(A(07))*exp(-(delta_theta(22)+delta_beta(22)*z)*1i )*ExcSet.HP...
        + 2*A(04)*A(07)*conj(A(02))*exp( (delta_theta(17)+delta_beta(17)*z)*1i )*ExcSet.HP...
        + 2*A(03)*A(04)*conj(A(01))*exp( (delta_theta(06)+delta_beta(06)*z)*1i )*ExcSet.HP...
        + 2*A(02)*A(04)*conj(A(03))*exp(-(delta_theta(15)+delta_beta(15)*z)*1i )*ExcSet.HP...
        + 2*A(02)*A(03)*conj(A(05))*exp(-(delta_theta(18)+delta_beta(18)*z)*1i )*ExcSet.HP...
                       )*ExcSet.HS;

  % Coupled ordinary differential equation of Signal2 (A7) field
  dA(07)= - (NLM.alpha/2)*A(07)*ExcSet.HS ...
          + 1i*NLM.Gamma*( ...
            A(07)*(abs(A(07))^2+2*(abs(A(01)          )^2+abs(A(02)          )^2 ...
                                  +abs(A(03)*ExcSet.HP)^2+abs(A(04)*ExcSet.HP)^2 ...
                                  +abs(A(05)*ExcSet.HS)^2+abs(A(06)*ExcSet.HS)^2))...
        +   A(03)*A(03)*conj(A(06))*exp(-(delta_theta(21)+delta_beta(21)*z)*1i )*ExcSet.HS...
        +   A(02)*A(02)*conj(A(03))*exp(-(delta_theta(16)+delta_beta(16)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(02)*conj(A(06))*exp(-(delta_theta(13)+delta_beta(13)*z)*1i )*ExcSet.HS...
        + 2*A(01)*A(03)*conj(A(04))*exp(-(delta_theta(10)+delta_beta(10)*z)*1i )*ExcSet.HP...
        + 2*A(01)*A(05)*conj(A(03))*exp(-(delta_theta(08)+delta_beta(08)*z)*1i )*ExcSet.HP...
        + 2*A(02)*A(03)*conj(A(01))*exp( (delta_theta(02)+delta_beta(02)*z)*1i )*ExcSet.HS...
        + 2*A(02)*A(05)*conj(A(07))*exp(-(delta_theta(14)+delta_beta(14)*z)*1i )*ExcSet.HP...
        + 2*A(02)*A(06)*conj(A(04))*exp(-(delta_theta(17)+delta_beta(17)*z)*1i )*ExcSet.HP...
        + 2*A(03)*A(05)*conj(A(02))*exp(-(delta_theta(19)+delta_beta(19)*z)*1i )*ExcSet.HP...
        + 2*A(04)*A(05)*conj(A(06))*exp(-(delta_theta(22)+delta_beta(22)*z)*1i )*ExcSet.HP...
        + 2*A(05)*A(06)*conj(A(01))*exp( (delta_theta(11)+delta_beta(11)*z)*1i )*ExcSet.HP...
                       )*ExcSet.HS; 
  % end of if-else SimMod.FLD_nPWR == 1
  %-----------------------------------------------------------------------------
  
else % SimMod.FLD_nPWR == 0, calculate propagation using Power and Phase

  if SimMod.PSA_nPIA == 1, % Calculate propagation in PSA mode

    % Initial differential equation for power of 7 waves
    % Initial differential equation for relative phase of 13+9 FWM processes
    dA    = zeros(29,1);

    % Differential equation of power P1 of A1/Pump1
    dA( 1)= - HNLF.alpha*E(1) ...
          - 2*HNLF.Gamma*( ...
              E(2)*sqrt(E(1)*E(5)          )*sin( (E(08)))...
          +   E(3)*sqrt(E(1)*E(2)          )*sin( (E(12)))...
          +   E(6)*sqrt(E(1)*E(4)          )*sin( (E(19)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(7))*sin( (E(09)))...
          + 2*E(1)*sqrt(E(2)*E(4)          )*sin( (E(10)))...
          + 2     *sqrt(E(1)*E(2)*E(4)*E(5))*sin( (E(16)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(6))*sin( (E(13)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(7))*sin( (E(17)))...
          + 2     *sqrt(E(1)*E(3)*E(5)*E(7))*sin( (E(15)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(6))*sin( (E(11)))...
          + 2*E(1)*sqrt(E(3)*E(6)          )*sin( (E(14)))...
          + 2     *sqrt(E(1)*E(5)*E(6)*E(7))*sin( (E(18)))...
          + 2     *sqrt(E(1)*E(2)*E(6)*E(7))*sin( (E(20)))...
                    );

    % Differential equation of power P2 of A2/Pump2
    dA( 2)= - HNLF.alpha*E(2) ...
          - 2*HNLF.Gamma*( ...
              E(1)*sqrt(E(2)*E(4)          )*sin(-(E(10)))...
          +   E(3)*sqrt(E(1)*E(2)          )*sin( (E(12)))...
          +   E(7)*sqrt(E(2)*E(5)          )*sin( (E(21)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(6))*sin(-(E(11)))...
          + 2*E(2)*sqrt(E(1)*E(5)          )*sin(-(E(08)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(7))*sin(-(E(09)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(7))*sin(-(E(26)))...
          + 2     *sqrt(E(2)*E(3)*E(4)*E(6))*sin( (E(22)))...
          + 2*E(2)*sqrt(E(3)*E(7)          )*sin( (E(23)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(6))*sin( (E(25)))...
          + 2     *sqrt(E(1)*E(2)*E(4)*E(5))*sin( (E(16)))...
          + 2     *sqrt(E(2)*E(4)*E(6)*E(7))*sin( (E(24)))...
          + 2     *sqrt(E(1)*E(2)*E(6)*E(7))*sin( (E(20)))...
                    ); 

    % Differential equation of power S&I of A3/Signal&Ilder                
    dA( 3)= - HNLF.alpha*E(3) ...
          - 2*HNLF.Gamma*( ...
              E(1)*sqrt(E(3)*E(6)          )*sin(-(E(14)))... 
          +   E(2)*sqrt(E(3)*E(7)          )*sin(-(E(23)))...
          + 2*E(3)*sqrt(E(1)*E(2)          )*sin(-(E(12)))... 
          + 2     *sqrt(E(1)*E(3)*E(5)*E(7))*sin(-(E(15)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(6))*sin(-(E(13)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(7))*sin(-(E(09)))...
          + 2     *sqrt(E(2)*E(3)*E(4)*E(6))*sin(-(E(22)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(6))*sin( (E(11)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(7))*sin( (E(26)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(7))*sin( (E(17)))...
          + 2*E(3)*sqrt(E(4)*E(5)          )*sin( (E(27)))...
          + 2*E(3)*sqrt(E(6)*E(7)          )*sin( (E(28)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(6))*sin( (E(25)))...
                    );   

    % Differential equation of power P3 of A4/Pump3
    dA( 4)= - HNLF.alpha*E(4) ...
          - 2*HNLF.Gamma*( ...
              E(1)*sqrt(E(2)*E(4)          )*sin(-(E(10)))...
          +   E(3)*sqrt(E(4)*E(5)          )*sin(-(E(27)))...
          +   E(6)*sqrt(E(1)*E(4)          )*sin( (E(19)))...
          + 2     *sqrt(E(1)*E(2)*E(4)*E(5))*sin(-(E(16)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(7))*sin(-(E(17)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(6))*sin(-(E(13)))...
          + 2     *sqrt(E(2)*E(3)*E(4)*E(6))*sin( (E(22)))...
          + 2     *sqrt(E(2)*E(4)*E(6)*E(7))*sin(-(E(24)))...
          + 2     *sqrt(E(4)*E(5)*E(6)*E(7))*sin( (E(29)))...
                    );

    % Differential equation of power P4 of A5/Pump4                
    dA( 5)= - HNLF.alpha*E(5) ...
          - 2*HNLF.Gamma*( ...
              E(2)*sqrt(E(1)*E(5)          )*sin( (E(08)))...
          +   E(3)*sqrt(E(4)*E(5)          )*sin(-(E(27)))...
          +   E(7)*sqrt(E(2)*E(5)          )*sin( (E(21)))...
          + 2     *sqrt(E(1)*E(2)*E(4)*E(5))*sin(-(E(16)))...
          + 2     *sqrt(E(1)*E(5)*E(6)*E(7))*sin(-(E(18)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(6))*sin(-(E(25)))...
          + 2     *sqrt(E(4)*E(5)*E(6)*E(7))*sin( (E(29)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(7))*sin( (E(26)))...
          + 2     *sqrt(E(1)*E(3)*E(5)*E(7))*sin( (E(15)))...
                    ); 

    % Differential equation of power S1 of A6/Signal1                
    dA( 6)= - HNLF.alpha*E(6) ...
          - 2*HNLF.Gamma*( ...
          +   E(1)*sqrt(E(3)*E(6)          )*sin(-(E(14)))... 
          +   E(3)*sqrt(E(6)*E(7)          )*sin(-(E(28)))...
          + 2     *sqrt(E(1)*E(2)*E(6)*E(7))*sin(-(E(20)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(6))*sin(-(E(11)))...
          + 2*E(6)*sqrt(E(1)*E(4)          )*sin(-(E(19)))...
          + 2     *sqrt(E(1)*E(5)*E(6)*E(7))*sin(-(E(18)))...
          + 2     *sqrt(E(4)*E(5)*E(6)*E(7))*sin(-(E(29)))...
          + 2     *sqrt(E(2)*E(4)*E(6)*E(7))*sin( (E(24)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(6))*sin( (E(13)))...
          + 2     *sqrt(E(2)*E(3)*E(4)*E(6))*sin(-(E(22)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(6))*sin(-(E(25)))...
                    );

    % Differential equation of power S2 of A7/Signal2 
    dA( 7)= - HNLF.alpha*E(7) ...
          - 2*HNLF.Gamma*( ...
              E(3)*sqrt(E(6)*E(7)          )*sin(-(E(28)))...
          +   E(2)*sqrt(E(3)*E(7)          )*sin(-(E(23)))...
          + 2     *sqrt(E(1)*E(2)*E(6)*E(7))*sin(-(E(20)))...
          + 2     *sqrt(E(1)*E(3)*E(4)*E(7))*sin(-(E(17)))...
          + 2     *sqrt(E(1)*E(3)*E(5)*E(7))*sin(-(E(15)))...
          + 2     *sqrt(E(1)*E(2)*E(3)*E(7))*sin( (E(09)))...
          + 2*E(7)*sqrt(E(2)*E(5)          )*sin(-(E(21)))...
          + 2     *sqrt(E(2)*E(4)*E(6)*E(7))*sin(-(E(24)))...
          + 2     *sqrt(E(2)*E(3)*E(5)*E(7))*sin(-(E(26)))...
          + 2     *sqrt(E(4)*E(5)*E(6)*E(7))*sin(-(E(29)))...
          + 2     *sqrt(E(1)*E(5)*E(6)*E(7))*sin( (E(18)))...
                    );

    %% 22 FWM: 13 NDFWM and 9 DFWM
    % Differential equation of phase-matching(01) + relative-phase(01) 
    % FWM Phase = ( 2b - f - a )z + 2*dPhi(2) - dPhi(1) - dPhi(5)
    dA(08)= delta_beta(01) ...
          +         2*dPhi(2) - dPhi(1) - dPhi(5);

    % Differential equation of phase-matching(02) + relative-phase(02) 
    % FWM Phase = ( 2b - f - a )z + dPhi(2) + dPhi(3) - dPhi(1) - dPhi(7)
    dA(09)= delta_beta(02) ...
          + dPhi(2) + dPhi(3) - dPhi(1) - dPhi(7);

    % Differential equation of phase-matching(03) + relative-phase(03) 
    % FWM Phase = ( 2b - f - a )z + dPhi(2) + dPhi(3) - dPhi(1) - dPhi(7)
    dA(10)= delta_beta(03) ...
          + dPhi(2) + dPhi(4)        -  2*dPhi(1); 

    % Differential equation of relative phase (04)
    dA(11)= delta_beta(04) ...
          + dPhi(2) + dPhi(6) - dPhi(1) - dPhi(3);

    % Differential equation of relative phase (05)
    dA(12)= delta_beta(05) ...
          +         2*dPhi(3) - dPhi(1) - dPhi(2);

    % Differential equation of relative phase (06)
    dA(13)= delta_beta(06) ...
          + dPhi(3) + dPhi(4) - dPhi(1) - dPhi(6);

    % Differential equation of relative phase (07)
    dA(14)= delta_beta(07) ...
          + dPhi(3) + dPhi(6)        -  2*dPhi(1); 

    % Differential equation of relative phase (08)
    dA(15)= delta_beta(08) ...
          + dPhi(3) + dPhi(7) - dPhi(1) - dPhi(5);

    % Differential equation of relative phase (09)
    dA(16)= delta_beta(09) ...
          + dPhi(4) + dPhi(5) - dPhi(1) - dPhi(2);

    % Differential equation of relative phase (10)
    dA(17)= delta_beta(10) ...
          + dPhi(4) + dPhi(7) - dPhi(1) - dPhi(3);

    % Differential equation of relative phase (11)
    dA(18)= delta_beta(11) ...
          + dPhi(5) + dPhi(6) - dPhi(1) - dPhi(7);

    % Differential equation of relative phase (12)
    dA(19)= delta_beta(12) ...
          +         2*dPhi(6) - dPhi(1) - dPhi(4);

    % Differential equation of relative phase (13)
    dA(20)= delta_beta(13) ...
          + dPhi(6) + dPhi(7) - dPhi(1) - dPhi(2);

    % Differential equation of relative phase (14)
    dA(21)= delta_beta(14) ...
          +         2*dPhi(7) - dPhi(2) - dPhi(5);

    % Differential equation of relative phase (15)
    dA(22)= delta_beta(15) ...
          + dPhi(3) + dPhi(6) - dPhi(2) - dPhi(4);

    % Differential equation of relative phase (16)
    dA(23)= delta_beta(16) ...
          + dPhi(3) + dPhi(7)        -  2*dPhi(2);

    % Differential equation of relative phase (17)
    dA(24)= delta_beta(17) ...
          + dPhi(4) + dPhi(7) - dPhi(2) - dPhi(6);

    % Differential equation of relative phase (18)
    dA(25)= delta_beta(18) ...
          + dPhi(5) + dPhi(6) - dPhi(2) - dPhi(3);

    % Differential equation of relative phase (19)
    dA(26)= delta_beta(19) ...
          + dPhi(2) + dPhi(7) - dPhi(3) - dPhi(5);

    % Differential equation of relative phase (20)
    dA(27)= delta_beta(20) ...
          + dPhi(4) + dPhi(5)        -  2*dPhi(3);

    % Differential equation of relative phase (21)
    dA(28)= delta_beta(21) ...
          + dPhi(6) + dPhi(7)        -  2*dPhi(3);

    % Differential equation of relative phase (22)
    dA(29)= delta_beta(22) ...
          + dPhi(6) + dPhi(7) - dPhi(4) - dPhi(5); 

  else % SimMod.PSA_nPIA == 0, calculate in PIA mode

  end; % end of if-else SimMod.PSA_nPIA
  %-----------------------------------------------------------------------------
end % End of if-else SimMod.FLD_nPWR
%-------------------------------------------------------------------------------
    
end % End of function dA_dz_7w.m
%-------------------------------------------------------------------------------

%%
%-------------------------------------------------------------------------------
% end of file dA_dz_7w.m
%-------------------------------------------------------------------------------