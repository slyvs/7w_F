%-------------------------------------------------------------------------------
% Beginning of file f_Cal_beta.m
%-------------------------------------------------------------------------------
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
function [beta0, beta1, beta2, beta3, beta4] = f_Cal_beta(frq)
% freq : vector of beta 
% beta : Beta Taylor Power Series expand

%% Global constant 
global c;          %

%% Simulation paramters 
global SimMod;
global HNLF;          % HNLF parameters
global ExcSet;   % 

%% Local varibles 

load('PCW.mat');
%-------------------------------------------------------------------------------

delta_frq       = frq - v_fq_model_interp;
[min_df, min_n] = min( abs(delta_frq) );
adj_frq         = v_fq_model_interp(min_n);
beta0           = 0;
beta1           = beta1_model_interp(min_n);
beta2           = beta2_model_interp(min_n);
beta3           = beta3_model_interp(min_n);
beta4           = beta4_model_interp(min_n);

end % end of function f_Cal_beta
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% end of file f_Cal_beta.m
%-------------------------------------------------------------------------------