%-------------------------------------------------------------------------------
% Beginning of file f_Nis_Gen.m
%-------------------------------------------------------------------------------
%
% Signal generation 
% 

function [ ] = f_Nis_Gen(nis_num,j3)

global SimMod;
global Exc;
global Nis;

%% 
% Nis.wgn1 = randn(Nis.Num,1); 
% Nis.wgn2 = randn(Nis.Num,1); 

Nis.sig_SNR     = 20; % [dB]

% Nis.sig_pow_mea_dBm = + 9; % dBm
% Nis.sig_pow_mea     = 10.^((Nis.sig_pow_mea_dBm - 30)/10);
%-------------------------------------------------------------------------------
%% 
% Nis.pow_nis_mea = 0; 
% Nis.pow_nis_var = Exc.Pwr(j3,3)/(10^(Nis.sig_SNR/10)); % Exc.Pwr(j3,3) 
% Nis.pow_nis_wgn = zeros(Nis.Num,1); 
% Nis.pow_nis_wgn = Nis.pow_nis_mea + sqrt(Nis.pow_nis_var).*Nis.wgn1; % N(pow_mea,pow_var) 
% Nis.pow_nis_wgn(end-SimMod.Phase_num+1:end,1) = 0; 
% 
% Nis.pha_nis_mea = 0; 
% Nis.pha_nis_var = (pi/2)/250; % (10^(Nis.sig_SNR/10))
% Nis.pha_nis_wgn = zeros(Nis.Num,1); 
% Nis.pha_nis_wgn = Nis.pha_nis_mea + sqrt(Nis.pha_nis_var).*Nis.wgn2; % N(pha_mea,pha_var) 
% Nis.pha_nis_wgn(end-SimMod.Phase_num+1:end,1) = 0; 
%-------------------------------------------------------------------------------
%% AWGN in channel 
Nis.pow_nis_var = Exc.Pwr(j3,3)/(10^(Nis.sig_SNR/10)); % Exc.Pwr(j3,3) 
Nis.nis_awgn    = zeros(Nis.Num,1); 
Nis.nis_awgn    = sqrt(Nis.pow_nis_var).*Nis.wgn1 + 1j*sqrt(Nis.pow_nis_var).*Nis.wgn2; 
Nis.nis_awgn(end-SimMod.Phase_num+1:end,1) = 0; 
%-------------------------------------------------------------------------------

end % end of function f_Nis_Gen
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% end of file f_Nis_Gen.m
%-------------------------------------------------------------------------------