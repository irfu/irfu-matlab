function delta_off = mms_sdp_comp_delta_off
%MMS_SDP_COMP_DELTA_OFF  compute DELTA_OFF for datamanager
%
%  delta_off = mms_sdp_comp_delta_off()
%
%  Compute delta offsets

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

global MMS_CONST; if isempty(MMS_CONST), MMS_CONST = mms_constants(); end


spinfits = mms_sdp_datamanager('spinfits');
if mms_is_error(spinfits)
  errStr='Bad DCE input, cannot proceed.';
  irf.log('critical',errStr); error(errStr);
end
    
Delta_p12_p34 = double(spinfits.sfit.e12(:,2:3)) - ...
  double(spinfits.sfit.e34(:,2:3));

delta_off = median(Delta_p12_p34(:,1)) + median(Delta_p12_p34(:,2))*1j;