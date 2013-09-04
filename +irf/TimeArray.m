classdef TimeArray < irf.Time
  %TimeArray Class representing a time array
  %   Many time points
  
% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------
  
  properties
    dt % int64 offset in nano-seconds
  end
  
  methods
    function ta = TimeArray(inp)
      ta = ta@irf.Time(inp(1));
      if isa(inp,'irf.TimeArray')
        ta.tt2000 = inp.tt2000;
        ta.dt = inp.dt;
      elseif isa(inp,'double')
        if min(size(inp))>1
          error('irf:TimeArray:TimeArray:badInputs',...
            'ISDAT epoch input (double) must be a columt or raw vector')
        end
        if size(inp,2)~=1, inp = inp'; end % to column
        ta.dt = int64((inp - inp(1)) * 1e9 );
      elseif isa(inp,'int64')
        if min(size(inp))>1
          error('irf:TimeArray:TimeArray:badInputs',...
            'epoch2000 nano-seconds input (int64) must be a column/row vector')
        end
        if size(inp,2)~=1, inp = inp'; end % to column
        ta.dt = inp - inp(1);
      else
        error('irf:TimeArray:TimeArray:badInputs',...
          'unknown input type')
      end
    end
    
    function display(ta)
      n = ta.length();
      fprintf('irf.TimeArray : %d points\n',n);
      if n>5
        for i = 1:2, displayTime, end
        fprintf('... skipped %d points ...\n',n-4);
        for i = n-1:n, displayTime, end
      else
        for i = 1:n, displayTime, end
      end
      function displayTime
        s_ttt = encodett2000(ta.tt2000+ta.dt(i));
        fprintf('%s\n',s_ttt{:})
      end
    end
    
    function r = eq(ta1,ta2)
      r = ( (ta1.tt2000==ta2.tt2000) && length(ta1.dt)==length(ta2.dt) && ...
        all(ta1.dt==ta2.dt) );
    end
    
    function l = length(ta)
      l = length(ta.dt);
    end
    
    function e = toEpoch(ta)
      epoch0 = toEpoch@irf.Time(ta);
      e = double(ta.dt)*1e-9 + epoch0;
    end
    
    function s = toUTC(ta)
      s = encodett2000(ta.tt2000+ta.dt);
    end
  end % Methods
  
end

