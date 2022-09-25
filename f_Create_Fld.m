%-------------------------------------------------------------------------------
% Beginning of file f_Create_Fld.m
%-------------------------------------------------------------------------------
%

function [ ] = f_Create_Fld(j1,j2,j3,j4,j5) 
% j1 = 1 : SimMod.WlOfs_num, 
% j2 = 1 : SimMod.WlSep_num, 
% j3 = 1 : SimMod.Power_num, 
% j4 = 1 : SimMod.Phase_num, 
% j5 = 1 : SimMod.PropL_num, 

%% Global variables 
global SimMod; 
global ExcSet; 
global Exc; 
global Nis;

% Exc.Pwr;
% Exc.Phi;
% Exc.Fld;

% Exc.Pwr_P1(j3,j4);
% Exc.Pwr_P2(j3,j4);
% Exc.Pwr_SI(j3,j4);
% Exc.Pwr_P3(j3,j4);
% Exc.Pwr_P4(j3,j4);
% Exc.Pwr_S1(j3,j4);
% Exc.Pwr_S2(j3,j4);
% 
% Exc.Phs_P1(j3,j4);
% Exc.Phs_P2(j3,j4);
% Exc.Phs_SI(j3,j4);
% Exc.Phs_P3(j3,j4);
% Exc.Phs_P4(j3,j4);
% Exc.Phs_S1(j3,j4);
% Exc.Phs_S2(j3,j4);
% 
% Exc.Fld_P1(j3,j4);
% Exc.Fld_P2(j3,j4);
% Exc.Fld_SI(j3,j4);
% Exc.Fld_P3(j3,j4);
% Exc.Fld_P4(j3,j4);
% Exc.Fld_S1(j3,j4);
% Exc.Fld_S2(j3,j4);
%-------------------------------------------------------------------------------

%% Cal: Fields according to initial powers and phases, w/o Noise 
Exc.Fld(1,1) = sqrt(Exc.Pwr(j3,1))*cos(Exc.Phi(j4,1)) + ...
            1i*sqrt(Exc.Pwr(j3,1))*sin(Exc.Phi(j4,1)); % V, Initial field of Pump1  (P1), non-degenerate 
Exc.Fld(1,2) = sqrt(Exc.Pwr(j3,2))*cos(Exc.Phi(j4,2)) + ...
            1i*sqrt(Exc.Pwr(j3,2))*sin(Exc.Phi(j4,2)); % V, Initial field of Pump2  (P2), non-degenerate
Exc.Fld(1,3) = sqrt(Exc.Pwr(j3,3))*cos(Exc.Phi(j4,3)) + ...
            1i*sqrt(Exc.Pwr(j3,3))*sin(Exc.Phi(j4,3)); % V, Initial field of Signal0(S0), degenerate signal & Idler 
Exc.Fld(1,4) = sqrt(Exc.Pwr(j3,4))*cos(Exc.Phi(j4,4)) + ...
            1i*sqrt(Exc.Pwr(j3,4))*sin(Exc.Phi(j4,4)); % V, Initial field of Pump4, (P4), A(4)
Exc.Fld(1,5) = sqrt(Exc.Pwr(j3,5))*cos(Exc.Phi(j4,5)) + ...
            1i*sqrt(Exc.Pwr(j3,5))*sin(Exc.Phi(j4,5)); % V, Initial field of Pump3, (P3), A(5)
Exc.Fld(1,6) = sqrt(Exc.Pwr(j3,6))*cos(Exc.Phi(j4,6)) + ...
            1i*sqrt(Exc.Pwr(j3,6))*sin(Exc.Phi(j4,6)); % V, Initial field of Signal1(S1), A(6)
Exc.Fld(1,7) = sqrt(Exc.Pwr(j3,7))*cos(Exc.Phi(j4,7)) + ...
            1i*sqrt(Exc.Pwr(j3,7))*sin(Exc.Phi(j4,7)); % V, Initial field of Signal2(S2), A(7)
%-------------------------------------------------------------------------------
%% Cal: Fields according to initial powers and phases, w/  non-additive Noise 
% j4: Noise_num
% Exc.Fld(1,1) = sqrt(Exc.Pwr(j3,1)                      )*cos(Exc.Phi(j4,1)                      ) + ...
%             1i*sqrt(Exc.Pwr(j3,1)                      )*sin(Exc.Phi(j4,1)                      ); % V, Initial field of Pump1  (P1), non-degenerate 
% Exc.Fld(1,2) = sqrt(Exc.Pwr(j3,2)                      )*cos(Exc.Phi(j4,2)                      ) + ...
%             1i*sqrt(Exc.Pwr(j3,2)                      )*sin(Exc.Phi(j4,2)                      ); % V, Initial field of Pump2  (P2), non-degenerate
% Exc.Fld(1,3) = sqrt(Exc.Pwr(j3,3)+Nis.pow_nis_wgn(j4,1))*cos(Exc.Phi(j4,3)+Nis.pha_nis_wgn(j4,1)) + ...
%             1i*sqrt(Exc.Pwr(j3,3)+Nis.pow_nis_wgn(j4,1))*sin(Exc.Phi(j4,3)+Nis.pha_nis_wgn(j4,1)); % + ...
% Exc.Fld(1,4) = sqrt(Exc.Pwr(j3,4)                      )*cos(Exc.Phi(j4,4)                      ) + ...
%             1i*sqrt(Exc.Pwr(j3,4)                      )*sin(Exc.Phi(j4,4)                      ); % V, Initial field of Pump4, (P4), A(4)
% Exc.Fld(1,5) = sqrt(Exc.Pwr(j3,5)                      )*cos(Exc.Phi(j4,5)                      ) + ...
%             1i*sqrt(Exc.Pwr(j3,5)                      )*sin(Exc.Phi(j4,5)                      ); % V, Initial field of Pump3, (P3), A(5)
% Exc.Fld(1,6) = sqrt(Exc.Pwr(j3,6)                      )*cos(Exc.Phi(j4,6)                      ) + ...
%             1i*sqrt(Exc.Pwr(j3,6)                      )*sin(Exc.Phi(j4,6)                      ); % V, Initial field of Signal1(S1), A(6)
% Exc.Fld(1,7) = sqrt(Exc.Pwr(j3,7)                      )*cos(Exc.Phi(j4,7)                      ) + ...
%             1i*sqrt(Exc.Pwr(j3,7)                      )*sin(Exc.Phi(j4,7)                      ); % V, Initial field of Signal2(S2), A(7)
%-------------------------------------------------------------------------------
%% Cal: Fields according to initial powers and phases, w/  additive white gaussian Noise 
% j4: Noise_num 
% Exc.Fld(1,1) = sqrt(Exc.Pwr(j3,1))*cos(Exc.Phi(j4,1)) + ...
%             1i*sqrt(Exc.Pwr(j3,1))*sin(Exc.Phi(j4,1)); % V, Initial field of Pump1  (P1), non-degenerate 
% Exc.Fld(1,2) = sqrt(Exc.Pwr(j3,2))*cos(Exc.Phi(j4,2)) + ...
%             1i*sqrt(Exc.Pwr(j3,2))*sin(Exc.Phi(j4,2)); % V, Initial field of Pump2  (P2), non-degenerate
% Exc.Fld(1,3) = sqrt(Exc.Pwr(j3,3))*cos(Exc.Phi(j4,3)) + ...
%             1i*sqrt(Exc.Pwr(j3,3))*sin(Exc.Phi(j4,3)) + ...
%             Nis.nis_awgn(j4,1); % + ...
% Exc.Fld(1,4) = sqrt(Exc.Pwr(j3,4))*cos(Exc.Phi(j4,4)) + ...
%             1i*sqrt(Exc.Pwr(j3,4))*sin(Exc.Phi(j4,4)); % V, Initial field of Pump4, (P4), A(4)
% Exc.Fld(1,5) = sqrt(Exc.Pwr(j3,5))*cos(Exc.Phi(j4,5)) + ...
%             1i*sqrt(Exc.Pwr(j3,5))*sin(Exc.Phi(j4,5)); % V, Initial field of Pump3, (P3), A(5)
% Exc.Fld(1,6) = sqrt(Exc.Pwr(j3,6))*cos(Exc.Phi(j4,6)) + ...
%             1i*sqrt(Exc.Pwr(j3,6))*sin(Exc.Phi(j4,6)); % V, Initial field of Signal1(S1), A(6)
% Exc.Fld(1,7) = sqrt(Exc.Pwr(j3,7))*cos(Exc.Phi(j4,7)) + ...
%             1i*sqrt(Exc.Pwr(j3,7))*sin(Exc.Phi(j4,7)); % V, Initial field of Signal2(S2), A(7)
%-------------------------------------------------------------------------------
%===============================================================================
%% Switch SimMod.Wave_num (N/A) 
% switch SimMod.Wave_num,  
%   case  3,
%     Exc.Fld(1,1) = sqrt(Exc.Pwr(j3,1))*cos(Exc.Phi(j4,1)) + ...
%             1i*sqrt(Exc.Pwr(j3,1))*sin(Exc.Phi(j4,1)); % V, Initial field of Pump1  (P1), non-degenerate
%     Exc.Fld(1,2) = sqrt(Exc.Pwr(j3,2))*cos(Exc.Phi(j4,2)) + ...
%                 1i*sqrt(Exc.Pwr(j3,2))*sin(Exc.Phi(j4,2)); % V, Initial field of Pump2  (P2), non-degenerate
%     Exc.Fld(1,3) = sqrt(Exc.Pwr(j3,3))*cos(Exc.Phi(j4,3)) + ...
%                 1i*sqrt(Exc.Pwr(j3,3))*sin(Exc.Phi(j4,3)); % V, Initial field of Signal0(S0), degenerate signal & Idler    
%   %-----------------------------------------------------------------------------
%   case  5,
%     Exc.Fld(1,1) = sqrt(Exc.Pwr(j3,1))*cos(Exc.Phi(j4,1)) + ...
%                 1i*sqrt(Exc.Pwr(j3,1))*sin(Exc.Phi(j4,1)); % V, Initial field of Pump1  (P1), non-degenerate 
%     Exc.Fld(1,2) = sqrt(Exc.Pwr(j3,2))*cos(Exc.Phi(j4,2)) + ...
%                 1i*sqrt(Exc.Pwr(j3,2))*sin(Exc.Phi(j4,2)); % V, Initial field of Pump2  (P2), non-degenerate
%     Exc.Fld(1,3) = sqrt(Exc.Pwr(j3,3))*cos(Exc.Phi(j4,3)) + ...
%                 1i*sqrt(Exc.Pwr(j3,3))*sin(Exc.Phi(j4,3)); % V, Initial field of Signal0(S0), degenerate signal & Idler 
%     Exc.Fld(1,4) = sqrt(Exc.Pwr(j3,4))*cos(Exc.Phi(j4,4)) + ...
%                 1i*sqrt(Exc.Pwr(j3,4))*sin(Exc.Phi(j4,4)); % V, Initial field of Pump4, (P4), A(4)
%     Exc.Fld(1,5) = sqrt(Exc.Pwr(j3,5))*cos(Exc.Phi(j4,5)) + ...
%                 1i*sqrt(Exc.Pwr(j3,5))*sin(Exc.Phi(j4,5)); % V, Initial field of Pump3, (P3), A(5)
%     Exc.Fld(1,6) = sqrt(Exc.Pwr(j3,6))*cos(Exc.Phi(j4,6)) + ...
%                 1i*sqrt(Exc.Pwr(j3,6))*sin(Exc.Phi(j4,6)); % V, Initial field of Signal1(S1), A(6)
%     Exc.Fld(1,7) = sqrt(Exc.Pwr(j3,7))*cos(Exc.Phi(j4,7)) + ...
%                 1i*sqrt(Exc.Pwr(j3,7))*sin(Exc.Phi(j4,7)); % V, Initial field of Signal2(S2), A(7)
%   %-----------------------------------------------------------------------------
%   case  7,
%     Exc.Fld(1,1) = sqrt(Exc.Pwr(j3,1))*cos(Exc.Phi(j4,1)) + ...
%                 1i*sqrt(Exc.Pwr(j3,1))*sin(Exc.Phi(j4,1)); % V, Initial field of Pump1  (P1), non-degenerate 
%     Exc.Fld(1,2) = sqrt(Exc.Pwr(j3,2))*cos(Exc.Phi(j4,2)) + ...
%                 1i*sqrt(Exc.Pwr(j3,2))*sin(Exc.Phi(j4,2)); % V, Initial field of Pump2  (P2), non-degenerate
%     Exc.Fld(1,3) = sqrt(Exc.Pwr(j3,3))*cos(Exc.Phi(j4,3)) + ...
%                 1i*sqrt(Exc.Pwr(j3,3))*sin(Exc.Phi(j4,3)); % V, Initial field of Signal0(S0), degenerate signal & Idler 
%     Exc.Fld(1,4) = sqrt(Exc.Pwr(j3,4))*cos(Exc.Phi(j4,4)) + ...
%                 1i*sqrt(Exc.Pwr(j3,4))*sin(Exc.Phi(j4,4)); % V, Initial field of Pump4, (P4), A(4)
%     Exc.Fld(1,5) = sqrt(Exc.Pwr(j3,5))*cos(Exc.Phi(j4,5)) + ...
%                 1i*sqrt(Exc.Pwr(j3,5))*sin(Exc.Phi(j4,5)); % V, Initial field of Pump3, (P3), A(5)
%     Exc.Fld(1,6) = sqrt(Exc.Pwr(j3,6))*cos(Exc.Phi(j4,6)) + ...
%                 1i*sqrt(Exc.Pwr(j3,6))*sin(Exc.Phi(j4,6)); % V, Initial field of Signal1(S1), A(6)
%     Exc.Fld(1,7) = sqrt(Exc.Pwr(j3,7))*cos(Exc.Phi(j4,7)) + ...
%                 1i*sqrt(Exc.Pwr(j3,7))*sin(Exc.Phi(j4,7)); % V, Initial field of Signal2(S2), A(7)
%   %-----------------------------------------------------------------------------
% end; % end of switch SimMod.Wave_num
% %-------------------------------------------------------------------------------
% %===============================================================================

%% Storage the initial fleld, power, and phase 
Exc.Fld_P1(j3,j4) = Exc.Fld(1,1);
Exc.Fld_P2(j3,j4) = Exc.Fld(1,2);
Exc.Fld_SI(j3,j4) = Exc.Fld(1,3);
Exc.Fld_P3(j3,j4) = Exc.Fld(1,4);
Exc.Fld_P4(j3,j4) = Exc.Fld(1,5);
Exc.Fld_S1(j3,j4) = Exc.Fld(1,6);
Exc.Fld_S2(j3,j4) = Exc.Fld(1,7);

Exc.Pwr_P1(j3,j4) = abs(Exc.Fld(1,1))^2;
Exc.Pwr_P2(j3,j4) = abs(Exc.Fld(1,2))^2;
Exc.Pwr_SI(j3,j4) = abs(Exc.Fld(1,3))^2;
Exc.Pwr_P3(j3,j4) = abs(Exc.Fld(1,4))^2;
Exc.Pwr_P4(j3,j4) = abs(Exc.Fld(1,5))^2;
Exc.Pwr_S1(j3,j4) = abs(Exc.Fld(1,6))^2;
Exc.Pwr_S2(j3,j4) = abs(Exc.Fld(1,7))^2;

Exc.Phs_P1(j3,j4) = angle(Exc.Fld(1,1)); % mod( ,pi);
Exc.Phs_P2(j3,j4) = angle(Exc.Fld(1,2)); % mod( ,pi);
Exc.Phs_SI(j3,j4) = angle(Exc.Fld(1,3)); % mod( ,pi);
Exc.Phs_P3(j3,j4) = angle(Exc.Fld(1,4)); % mod( ,pi);
Exc.Phs_P4(j3,j4) = angle(Exc.Fld(1,5)); % mod( ,pi);
Exc.Phs_S1(j3,j4) = angle(Exc.Fld(1,6)); % mod( ,pi);
Exc.Phs_S2(j3,j4) = angle(Exc.Fld(1,7)); % mod( ,pi);
%-------------------------------------------------------------------------------

end % end of function f_Create_Fld
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% end of file f_Create_Fld.m
%-------------------------------------------------------------------------------