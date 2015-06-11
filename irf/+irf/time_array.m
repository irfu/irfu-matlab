function TaTT = time_array(start,dt)
%TIME_ARRAY  Factory of time arrays (EpochTT)
%
% TaTT = irf.time_array(start,dt)
% TaTT = irf.time_array(times)

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,2)

if isa(start,'GenericTimeArray')
  %
elseif isa(start,'char')
  if nargin==2 && ~isvector(start) || ~ismatrix(start)
    errStr = 'START must be a string yyyy-mm-ddThh:mm:ss[.mmmuuunnnZ]';
    irf.log('critical', errStr), error(errStr)
  end
  %
elseif (isa(start,'double') || isa(start,'int64')) && ...
    ((nargin==2 && isscalar(start)) || isvector(start))
  %
else
  errStr = 'Unrecoglized input';
  irf.log('critical', errStr), error(errStr)
end

Epoch0 = EpochTT(start);
if nargin == 1, TaTT = Epoch0; return, end

if ~isvector(dt) || ~isnumeric(dt)
  errStr = 'DT must be a numeric vector';
  irf.log('critical', errStr), error(errStr)
end

if isa(dt,'int64')
  TaTT = EpochTT(Epoch0.epoch + dt);
elseif isa(dt,'double')
  TaTT = EpochTT(Epoch0.epoch + int64(dt*1e9));
end
  