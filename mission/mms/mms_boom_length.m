function [Blen] = mms_boom_length(tInp,sc_id)
%mms_boom_len   get MMS boom lengths
%
% [Blen] = mms_boom_len(t,sc_id)
%

narginchk(2,2)

if isa(tInp,'GenericTimeArray'), t = tInp.toEpochTT2000().epoch;
elseif isa(tInp,'int64'), t = tInp;
else
  errStr = 'T must be tt2000_ns (int64) or GenericTimeArray';
  irf.log('critical',errStr), error(errStr)
end
t = t(1);

% Table of boom lengths
switch sc_id
  case 1
    if t>=get_ns([ 2015 04 13 14 22 32 ]), Blen = [ 22 22 41 41  ];
    elseif t>=get_ns([ 2015 04 13 01 11 13 ]), Blen = [ 22 22 22 22  ];
    elseif t>=get_ns([ 2015 04 09 14 09 50 ]), Blen = [ 17 17 22 22  ];
    elseif t>=get_ns([ 2015 04 04 07 58 51 ]), Blen = [ 17 17 17 17  ];
    elseif t>=get_ns([ 2015 04 04 07 25 34 ]), Blen = [ 03 03 17 17  ];
    elseif t>=get_ns([ 2015 04 04 06 54 51 ]), Blen = [ 03 03 03 03  ];
    elseif t>=get_ns([ 2015 04 04 06 43 15 ]), Blen = [ 00 00 03 03  ];
    else
      Blen = [ 0 0 0 0 ];
    end
  case 2
    if t>=get_ns([ 2015 04 13 04 07 12 ]), Blen = [ 22 22 41 41  ];
    elseif t>=get_ns([ 2015 04 13 03 20 36 ]), Blen = [ 22 22 22 22  ];
    elseif t>=get_ns([ 2015 04 09 16 23 38 ]), Blen = [ 17 17 22 22  ];
    elseif t>=get_ns([ 2015 04 06 12 03 28 ]), Blen = [ 17 17 17 17  ];
    elseif t>=get_ns([ 2015 04 06 11 31 53 ]), Blen = [ 03 03 17 17  ];
    elseif t>=get_ns([ 2015 04 06 11 00 32 ]), Blen = [ 03 03 03 03  ];
    elseif t>=get_ns([ 2015 04 06 10 49 15 ]), Blen = [ 00 00 03 03  ];
    else
      Blen = [ 0 0 0 0 ];
    end
  case 3
    if t>=get_ns([ 2015 04 14 14 23 18 ]), Blen = [ 22 22 41 41  ];
    elseif t>=get_ns([ 2015 04 12 14 04 22 ]), Blen = [ 22 22 22 22  ];
    elseif t>=get_ns([ 2015 04 10 14 13 18 ]), Blen = [ 17 17 22 22  ];
    elseif t>=get_ns([ 2015 04 07 08 21 36 ]), Blen = [ 17 17 17 17  ];
    elseif t>=get_ns([ 2015 04 07 07 49 48 ]), Blen = [ 03 03 17 17  ];
    elseif t>=get_ns([ 2015 04 07 06 49 46 ]), Blen = [ 03 03 7.3 7.3  ];
    elseif t>=get_ns([ 2015 04 07 06 34 54 ]), Blen = [ 03 03 03 03  ];
    elseif t>=get_ns([ 2015 04 07 06 24 26 ]), Blen = [ 00 00 03 03  ];
    else
      Blen = [ 0 0 0 0 ];
    end
  case 4
    if t>=get_ns([ 2015 04 14 16 39 13 ]), Blen = [ 22 22 41 41  ];
    elseif t>=get_ns([ 2015 04 12 16 16 49 ]), Blen = [ 22 22 22 22  ];
    elseif t>=get_ns([ 2015 04 10 16 22 20 ]), Blen = [ 17 17 22 22  ];
    elseif t>=get_ns([ 2015 04 05 12 03 35 ]), Blen = [ 17 17 17 17  ];
    elseif t>=get_ns([ 2015 04 05 11 32 05 ]), Blen = [ 03 03 17 17  ];
    elseif t>=get_ns([ 2015 04 05 11 01 46 ]), Blen = [ 03 03 03 03  ];
    elseif t>=get_ns([ 2015 04 05 10 50 56 ]), Blen = [ 00 00 03 03  ];
    else
      Blen = [ 0 0 0 0 ];
    end
  otherwise
    error('Invalid sc_id')
end

  function ns = get_ns(vec)
    % Convert date to tt2000_ns
    ns = irf_time(vec,'vector6>ttns');
  end
end

