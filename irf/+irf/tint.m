function TintTT2000 = tint(start,stop)
%TINT  Factory of time intervals (EpochTT2000)
%
%  TintTT2000 = tint(start,stop)
%  TintTT2000 = tint('utc/utc')
%  TintTT2000 = tint(start, duration_sec_double)
%  TintTT2000 = tint(start, duration_ns_int64)

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,2)

switch nargin
  case 1
    if ischar(start)
      errStr = ['Expecting UTC time range: '...
        'yyyy-mm-ddThh:mm:ss[.mmmuuunnnZ]/yyyy-mm-ddThh:mm:ss[.mmmuuunnnZ]'];
      if ~isvector(start), irf.log('critical',errStr), error(errStr), end
      ii = strfind(start(1,:),'/');
      if isempty(ii), irf.log('critical',errStr), error(errStr), end
      t1 = irf.utc_validate_and_pad(start(:,1:ii(1)-1));
      t2 = irf.utc_validate_and_pad(start(:,ii(1)+1:end));
      TintTT2000 = irf.time_array([t1; t2]);
    else
      if ~isvector(start)
        errStr = 'Expecting vector: [START STOP]';
        irf.log('critical',errStr), error(errStr)
      end
      TintTT2000 = irf.time_array(start([1 end]));
    end
  case 2
    if ischar(start) && ischar(stop)
      if ~isvector(start) || ~isvector(stop)
        errStr = 'Expecting UTC string: yyyy-mm-ddThh:mm:ss[.mmmuuunnnZ]';
        irf.log('critical',errStr), error(errStr)
      end
      TintTT2000 = irf.time_array([start; stop]);
    elseif (ischar(start)||isa(start,'GenericTimeArray')) && isnumeric(stop)
      TintTT2000 = irf.time_array(start, [0 stop]);
    else
      TintTT2000 = irf.time_array([start; stop]);
    end
  otherwise % should not be here
end

if TintTT2000(2) <= TintTT2000(1),
  errStr = 'TINT has negative or zero duration';
  irf.log('critical',errStr), error(errStr)
end