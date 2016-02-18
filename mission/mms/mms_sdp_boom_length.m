function [Blen] = mms_sdp_boom_length(scId,tInp)
%MMS_SDP_BOOM_LENGTH   get MMS SDP boom lengths
%
% [Blen] = mms_sdp_boom_length(scId)
%

narginchk(1,2)

% Table of boom lengths
switch scId
  case 1
    Blen(1).time=get_tt([2015 03 13 00 00 00]); Blen(1).len=[ 00 00 00 00 ];
    Blen(2).time=get_tt([2015 04 04 06 43 15]); Blen(2).len=[ 00 00 03 03 ];
    Blen(3).time=get_tt([2015 04 04 06 54 51]); Blen(3).len=[ 03 03 03 03 ];
    Blen(4).time=get_tt([2015 04 04 07 25 34]); Blen(4).len=[ 03 03 17 17 ];
    Blen(5).time=get_tt([2015 04 04 07 58 51]); Blen(5).len=[ 17 17 17 17 ];
    Blen(6).time=get_tt([2015 04 09 14 09 50]); Blen(6).len=[ 17 17 22 22 ];
    Blen(7).time=get_tt([2015 04 13 01 11 13]); Blen(7).len=[ 22 22 22 22 ];
    Blen(8).time=get_tt([2015 04 13 14 22 32]); Blen(8).len=[ 22 22 41 41 ];   
    Blen(9).time=get_tt([2015 04 20 23 21 06]); Blen(9).len=[ 22 22 57 57 ];   
    Blen(10).time=get_tt([2015 04 22 22 55 59]); Blen(10).len=[ 31.3 31.3 57 57 ];   
    Blen(11).time=get_tt([2015 04 22 23 17 22]); Blen(11).len=[ 41 41 57 57 ];   
    Blen(12).time=get_tt([2015 04 24 05 14 08]); Blen(12).len=[ 57 57 57 57 ];   
  case 2
    Blen(1).time=get_tt([2015 03 13 00 00 00]); Blen(1).len=[ 00 00 00 00 ];
    Blen(2).time=get_tt([2015 04 06 10 49 15]); Blen(2).len=[ 00 00 03 03 ];
    Blen(3).time=get_tt([2015 04 06 11 00 32]); Blen(3).len=[ 03 03 03 03 ];
    Blen(4).time=get_tt([2015 04 06 11 31 53]); Blen(4).len=[ 03 03 17 17 ];
    Blen(5).time=get_tt([2015 04 06 12 03 28]); Blen(5).len=[ 17 17 17 17 ];
    Blen(6).time=get_tt([2015 04 09 16 23 38]); Blen(6).len=[ 17 17 22 22 ];
    Blen(7).time=get_tt([2015 04 13 03 20 36]); Blen(7).len=[ 22 22 22 22 ];
    Blen(8).time=get_tt([2015 04 13 04 07 12]); Blen(8).len=[ 22 22 41 41 ];
    Blen(9).time=get_tt([2015 04 21 01 05 21]); Blen(9).len=[ 22 22 42.5 42.5 ];
    Blen(10).time=get_tt([2015 04 21 01 39 52]); Blen(10).len=[ 22 22 57 57 ];
    Blen(11).time=get_tt([2015 04 23 01 25 14]); Blen(11).len=[ 41 41 57 57 ];
    Blen(12).time=get_tt([2015 04 24 09 30 06]); Blen(12).len=[ 57 57 57 57 ];
  case 3
    Blen(1).time=get_tt([2015 03 13 00 00 00]); Blen(1).len=[ 00 00 00 00 ];
    Blen(2).time=get_tt([2015 04 07 06 24 26]); Blen(2).len=[ 00 00 03 03 ];
    Blen(3).time=get_tt([2015 04 07 06 34 54]); Blen(3).len=[ 03 03 03 03 ];
    Blen(4).time=get_tt([2015 04 07 06 49 46]); Blen(4).len=[ 03 03 7.3 7.3 ];
    Blen(5).time=get_tt([2015 04 07 07 49 48]); Blen(5).len=[ 03 03 17 17 ];
    Blen(6).time=get_tt([2015 04 07 08 21 36]); Blen(6).len=[ 17 17 17 17 ];
    Blen(7).time=get_tt([2015 04 10 14 13 18]); Blen(7).len=[ 17 17 22 22 ];
    Blen(8).time=get_tt([2015 04 12 14 04 22]); Blen(8).len=[ 22 22 22 22 ];
    Blen(9).time=get_tt([2015 04 14 14 23 18]); Blen(9).len=[ 22 22 41 41 ];
    Blen(10).time=get_tt([2015 04 19 23 37 15]); Blen(10).len=[ 41 41 41 41 ];
    Blen(11).time=get_tt([2015 04 21 22 53 03]); Blen(11).len=[ 41 41 41.2 41.2 ];
    Blen(12).time=get_tt([2015 04 21 23 24 42]); Blen(12).len=[ 41 41 57 57 ];
    Blen(13).time=get_tt([2015 04 24 00 57 12]); Blen(13).len=[ 57 57 57 57 ];
  case 4
    Blen(1).time=get_tt([2015 03 13 00 00 00]); Blen(1).len=[ 00 00 00 00 ];
    Blen(2).time=get_tt([2015 04 05 10 50 56]); Blen(2).len=[ 00 00 03 03 ];
    Blen(3).time=get_tt([2015 04 05 11 01 46]); Blen(3).len=[ 03 03 03 03 ];
    Blen(4).time=get_tt([2015 04 05 11 32 05]); Blen(4).len=[ 03 03 17 17 ];
    Blen(5).time=get_tt([2015 04 05 12 03 35]); Blen(5).len=[ 17 17 17 17 ];
    Blen(6).time=get_tt([2015 04 10 16 22 20]); Blen(6).len=[ 17 17 22 22 ];
    Blen(7).time=get_tt([2015 04 12 16 16 49]); Blen(7).len=[ 22 22 22 22 ];
    Blen(8).time=get_tt([2015 04 14 16 39 13]); Blen(8).len=[ 22 22 41 41 ];
    Blen(9).time=get_tt([2015 04 20 01 48 52]); Blen(9).len=[ 41 41 41 41 ];
    Blen(10).time=get_tt([2015 04 22 01 36 49]); Blen(10).len=[ 41 41 57 57 ];
    Blen(11).time=get_tt([2015 04 24 03 11 49]); Blen(11).len=[ 57 57 57 57 ];
  otherwise
    error('Invalid scId')
end
if nargin == 1, return, end

if isa(tInp,'GenericTimeArray'), t = EpochTT(tInp).epoch;
elseif isa(tInp,'int64'), t = tInp;
else
  errStr = 'T must be tt2000_ns (int64) or GenericTimeArray';
  irf.log('critical',errStr), error(errStr)
end

if isscalar(t), Blen = Blen(find_idx(t));
elseif isvector(t), Blen = Blen(find_idx_range(t));
else
  errStr = 'T(one point) or TINT (2 points) is expected';
  irf.log('critical',errStr), error(errStr)
end 
  function idx = find_idx(t)
    idx = arrayfun(@(x) x.time.epoch<=t,Blen);
    if ~any(idx)
      errS = 'Invalid time input (outside of range)';
      irf.log('critical',errS), error(errS)
    end
    idx = find(idx); idx = idx(end);
  end
  function idx = find_idx_range(t)
    errS = 'Invalid time input (outside of range)';
    idx1 = arrayfun(@(x) x.time.epoch<=t(1),Blen);
    if ~any(idx1), irf.log('critical',errS), error(errS), end
    idx1(find(idx1,1,'last')) = false;
    idx = arrayfun(@(x) x.time.epoch<=t(end),Blen);
    if ~any(idx), irf.log('critical',errS), error(errS), end
    idx(idx1) = false;
  end
  function ns = get_tt(vec)
    % Convert date to tt2000_ns
    ttns = irf_time(vec,'vector6>ttns');
    ns = EpochTT(ttns);
  end
end

