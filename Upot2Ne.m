function Ne = Upot2Ne( Upot );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function Ne = Upot2Ne( Upot );
%
%   Function that calculates the plasma density (Ne) for given EFW 
%   probes to spacecraft potential values.
%   Unfortnately we decided to take the easy simple linear fit.
%
%   Created by the good team of Andris Vaivads and Jan-Erik Wahlund at 
%   the Swedish Institute of Space Physics, Uppsala division, 2002.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CLUSTER SC potential and Ne based on Andris 
% Harri Laakso WHISPER data comparison.
% 2002.06.13
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Upot_ref = [  2  3  5 8  10  15   25   35];
  Ne_ref   = [100 36 10 4 2.5 1.3 0.45 0.17];

  Upot = abs( real( Upot) );
  Ne   = 10.^interp1( Upot_ref, log10(Ne_ref), Upot, 'linear', 'extrap' );

  return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
