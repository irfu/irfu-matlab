function TintTT = tint(start,stop)
%IRF.TINT  Factory of time intervals (EpochTT)
%
%  TintTT2000 = irf.tint(start,stop)
%  TintTT2000 = irf.tint('utc/utc')
%  TintTT2000 = irf.tint(start, duration_sec_double)
%  TintTT2000 = irf.tint(start, duration_ns_int64)

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
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
      t1 = GenericTimeArray.validate_and_pad_utc(start(:,1:ii(1)-1));
      t2 = GenericTimeArray.validate_and_pad_utc(start(:,ii(1)+1:end));
      TintTT = irf.time_array([t1; t2]);
    else
      if ~isvector(start)
        errStr = 'Expecting vector: [START STOP]';
        irf.log('critical',errStr), error(errStr)
      end
      TintTT = irf.time_array(start([1 end]));
    end
    
  case 2
    if ischar(start) && ischar(stop)
      if ~isvector(start) || ~isvector(stop)
        errStr = 'Expecting UTC string: yyyy-mm-ddThh:mm:ss[.mmmuuunnnZ]';
        irf.log('critical',errStr), error(errStr)
      end
      TintTT = irf.time_array([start; stop]);
    elseif (ischar(start)||isa(start,'GenericTimeArray')) && isnumeric(stop)
      TintTT = irf.time_array(start, [0 stop]);
    elseif isa(start,'GenericTimeArray') ...
        && isa(stop,'GenericTimeArray')
      TintTT = irf.tint([start.utc '/' stop.utc]);
    else
      TintTT = irf.time_array([start; stop]);
    end
  otherwise % should not be here
end

if TintTT(2) < TintTT(1)
  errStr = 'TINT has stop before start';
  irf.log('critical',errStr), error(errStr)
end