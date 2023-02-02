function n=c_efw_scp2ne(Vps,flag,coef)
%C_EFW_SCP2NE calculates the plasma density (Ne) for given EFW sc potential
%
% Ne = c_efw_scp2ne( Vps, [flag],[coef] );
%
%   Function that calculates the plasma density (Ne) for given EFW
%   probes to spacecraft potential values Vps.
%   Unfortnately we decided to take the easy simple linear fit.
%
%   Ne =  c_efw_scp2ne( Vps);
%   Vps = c_efw_scp2ne( Ne, -1 );
%
%   Vps - if it has two columns, then converts the second column
%   Vps - if negative values assumes those to be the negative of the spacecraft potential
%   flag
%        -1   convert from Ne to Vps
%   coef - multiplication factor for the final density, default = 1
%

%   Created by the good team of Andris Vaivads and Jan-Erik Wahlund at
%   the Swedish Institute of Space Physics, Uppsala division, 2002.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Vps_ref  = [  2  3  5 8  10  15   25   35];
%  Ne_ref   = [100 36 10 4 2.5 1.3 0.45 0.17];
%

% CLUSTER SC potential and Ne based on Andris
% Harri Laakso WHISPER data comparison.
% 2002.06.13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vps_ref  = [  2  3  5 8  10  15   25   35];
Ne_ref   = [100 36 10 4 2.5 1.3 0.45 0.17];

n           = Vps;
colind      = min(2,size(Vps,2)):size(Vps,2);

switch nargin
  case 1
    % Vps -> Ne
    flag=1;coef=1;
  case 2
    coef=1;
end

if flag==-1
  % Ne -> Vps
  warning off
  n(:,colind) = interp1( log10(Ne_ref), Vps_ref, log10(Vps(:,colind)/coef), 'linear', 'extrap' );
  warning on
else    % Vps -> Ne
  n(:,colind) = 10.^interp1( Vps_ref, log10(Ne_ref), abs(real(Vps(:,colind))), 'linear', 'extrap' )*coef;
  
end

return;
