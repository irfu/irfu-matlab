function utcNew = utc_validate_and_pad(utc)
% UTC_VALIDATE  validate UTC string and pad with zeros
%
% utcNew = utc_validate_and_pad(utc)
%
%  Validates that string is compliand with the UTC format:
%  yyyy-mm-ddThh:mm:ss.[mmmuuunnnZ] and pads with missing zeros and Z.

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

MAX_NUM_IDX = 29; idxDotIDX_DOT = 20;

if ~GenericTimeArray.validate_utc_time_str(utc)
  errStr = 'UTC must be a string: yyyy-mm-ddThh:mm:ss.[mmmuuunnnZ]';
  irf.log('critical',errStr), error(errStr)
end

utcNew = utc; lUtc = size(utc,2); flagAddZ = true;
if all(all(utc(:,end)=='Z')), lUtc = lUtc - 1; flagAddZ = false; end

if lUtc == MAX_NUM_IDX && ~flagAddZ, return, end
if lUtc == idxDotIDX_DOT-1, utcNew(:,idxDotIDX_DOT) = '.'; lUtc = lUtc + 1; end
if lUtc < MAX_NUM_IDX, utcNew(:,(lUtc+1):MAX_NUM_IDX) = '0'; end % Pad with zeros
utcNew(:,MAX_NUM_IDX+1) = 'Z';

end