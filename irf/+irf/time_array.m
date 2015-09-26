function TaTT = time_array(tStart,dtArray)
%TIME_ARRAY  Factory of time arrays (EpochTT)
%
% TaTT = irf.time_array(tStart,dtArray)
% TaTT = irf.time_array(tArray)
%
% Input: 
%   tStart is reference time given as either GenericTimeArray, UTC, or TT in
%     seconds (double) or nanoseconds (int64)
%   dtArray is an array of time values measured since tStart. dtArray can be
%     specified in seconds (double) or nanoseconds (int64). 
%   tArray is time array given either as GenericTimeArray, UTC or TT in
%     seconds (doubel or nanoseconds (int64)
% Output: 
%  TaTT time array as EpochTT
%
% Example:
%  TT = irf.time_array('2015-09-26T20:00:00',[0 10 20])
% 

narginchk(1,2)

if isa(tStart,'GenericTimeArray')
  %
elseif isa(tStart,'char')
  if nargin==2 && ~isvector(tStart) || ~ismatrix(tStart)
    errStr = 'START must be a string yyyy-mm-ddThh:mm:ss[.mmmuuunnnZ]';
    irf.log('critical', errStr), error(errStr)
  end
  %
elseif (isa(tStart,'double') || isa(tStart,'int64')) && ...
    ((nargin==2 && isscalar(tStart)) || isvector(tStart))
  %
else
  errStr = 'Unrecoglized input';
  irf.log('critical', errStr), error(errStr)
end

Epoch0 = EpochTT(tStart);
if nargin == 1, TaTT = Epoch0; return, end

if ~isvector(dtArray) || ~isnumeric(dtArray)
  errStr = 'DT must be a numeric vector';
  irf.log('critical', errStr), error(errStr)
end

if isa(dtArray,'int64')
  TaTT = EpochTT(Epoch0.epoch + dtArray);
elseif isa(dtArray,'double')
  TaTT = EpochTT(Epoch0.epoch + int64(dtArray*1e9));
end
  