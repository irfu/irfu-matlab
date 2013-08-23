classdef Time
  %TIME Time class
  %   This is one time point
  
% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

  properties
    tt2000 % int64 time in nanoseconds since 12:00, Jan 1, 2000
  end
  
  methods
    function t = Time(inp)
      if isa(inp,'irf.Time')
        t.tt2000 = inp.tt2000;
      elseif isa(inp,'double') % ISDAT epoch
        if all(size(1)==1) 
          t.tt2000 = int64( (inp-946728000.0)*1e9 );
        else
          error('ISDAT epoch input (double) must be scalar')
        end
      elseif isa(inp,'int64')
        if all(size(1)==1) % epoch2000 nano-seconds
          t.tt2000 = inp;
        else
          error('epoch2000 nano-seconds input (int64) must be scalar')
        end
      else
        error('unknown input type')
      end
    end
    
    function e = toEpoch(t)
      e = double(t.tt2000)*1e-9 + 946728000.0;
    end
    
    function s = toISO(t)
      s = epoch2iso(t.toEpoch);
    end
    
    function r = eq(t1,t2)
      r = (t1.tt2000 == t2.tt2000);
    end
    function display(t)
      fprintf('irf.Time : %s\n',t.toISO());
    end
  end
end

