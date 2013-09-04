classdef Time
  %TIME Basic class describing time
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
        if all(size(inp)==1) 
          t.tt2000 = parsett2000(epoch2iso(inp));
        else
          error('MATLAB:Time:Time:badInputs',...
            'ISDAT epoch input (double) must be scalar')
        end
      elseif isa(inp,'int64')
        if all(size(inp)==1) % epoch2000 nano-seconds
          t.tt2000 = inp;
        else
          error('MATLAB:Time:Time:badInputs',...
            'epoch2000 nano-seconds input (int64) must be scalar')
        end
      elseif ischar(inp)
        if inp(end)=='z' || inp(end)=='Z', inp(end)=''; end
        if length(inp)<21
          error('MATLAB:Time:Time:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnn')
        end
        t.tt2000 = parsett2000(inp);
      else
        error('MATLAB:Time:Time:badInputs',...
          'Unknown input type')
      end
    end
    
    function e = toEpoch(t)
      e = iso2epoch2(toUTC(t));
    end
    
    function s = toUTC(t)
      s_tmp = encodett2000(t.tt2000);
      s = s_tmp{:};
    end
    
    function r = eq(t1,t2)
      r = (t1.tt2000 == t2.tt2000);
    end
    
    function display(t)
      fprintf('irf.Time : %s\n',toUTC(t));
    end
    
  end % methods
end

