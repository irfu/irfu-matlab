function dE = mms_sdp_despin(e12,e34,phaseDeg)
%MMS_SDP_DESPIN  despin SDP data
%
% function dE = mms_sdp_despin(e12,e34,phaseDeg)
%
% Despin SDP data using phaseDeg [degrees]
%
% Input:
%   phaseDeg - spin phase in degrees
%   e12, e34 - electric field components in mV/m
%
%   Note: phaseDeg, e12 and e34 must be of the same size
%
% Output: 
%   Despun E (X and Y DSL)

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end

if ~(all(size(e12)==size(e34)) && all(size(e12)==size(phaseDeg)))
  errS = 'e12, e34 and phaseDeg must be of the same size';
  irf.log('critical',errS), error(errS)
end
  
phi_12 = MMS_CONST.Phaseshift.e12;
phi_34 = MMS_CONST.Phaseshift.e34; % Angles when phase=0 (X BSC direction)

compE = (e12*exp(-1i*phi_12) + e34*exp(-1i*phi_34)).*exp(1i*phaseDeg*pi/180);
dE = [real(compE) imag(compE)];

end