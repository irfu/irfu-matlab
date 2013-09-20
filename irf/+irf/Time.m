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
      %Create a time object
      %
      % t = Time(time_t)
      % t = Time(double_epoch)
      % t = Time(int64_epoch_tt2000)
      % t = Time(string_utc)
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
        if inp(end)=='Z' || inp(end)=='z', inp(end)=''; end
        s='0123456789';
        if length(inp)<19 || length(inp)>29 || ~all(inp([5 8])=='-') || ...
            inp(11)~='T' || ~all(inp([14 17])==':') || ...
            all(inp(1)~='12') || all(inp(2)~=s) || all(inp(3)~=s) || all(inp(4)~=s) || ...
            all(inp(6)~='01') || all(inp(7)~=s) || all(inp(6:7)==0) || ...
            all(inp(9)~='0123') || all(inp(10)~=s) || all(inp(9:10)==0) || ...
            all(inp(12)~='012') || all(inp(13)~=s) || ...
            all(inp(15)~='012345') || all(inp(16)~=s) || ...
            all(inp(18)~='0123456') || all(inp(19)~=s) || ...
            length(inp)>19 && ~all(double(inp(21:end))>47 & double(inp(21:end))<58)
          error('MATLAB:Time:Time:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        else
          yyyy = str2double(inp(1:4));
          if yyyy>2050 || yyyy<1950
            error('MATLAB:Time:Time:badInputs',...
              'YEAR must be between 1950 and 2050')
          end
          mo = str2double(inp(6:7));
          if mo > 12
            error('MATLAB:Time:Time:badInputs',...
              'MONTH must be between 1 and 12')
          end
          dd = str2double(inp(9:10));
          if any(mo==[1,3,5,7,8,10,12]), maxDay = 31;
          elseif any(mo==[4,6,9,11]),  maxDay = 30;
          elseif any(mo==2)
            if fix(yyyy/4)==yyyy/4, maxDay = 29; else maxDay = 28; end
          else
            maxDay = -1;
          end
          if dd>maxDay
            error('MATLAB:Time:Time:badInputs','Bad value for DAY')
          end
          hh = str2double(inp(12:13));
          if hh>23
            error('MATLAB:Time:Time:badInputs','Bad value for HOURS')
          end
          mm = str2double(inp(15:16));
          if mm>59
            error('MATLAB:Time:Time:badInputs','Bad value for MINUTES')
          end
          ss = str2double(inp(18:19));
          if ss>60
            error('MATLAB:Time:Time:badInputs','Bad value for SEC')
          elseif ss==60 && ~(hh==23 && mm==59 && (mo==6&&dd==30 || mo==12&&dd==31))
            error('MATLAB:Time:Time:badInputs',...
              'Bad value for SEC: Leap second can be only in the end of June of December')
          end
          if length(inp)==29
            % skip
          elseif length(inp)==19
            inp = [inp '.000000000'];
          else
            zeroPad = '000000000';
            inp = [inp zeroPad(29-length(inp))];
          end
        end
        t.tt2000 = parsett2000(inp);
      else
        error('MATLAB:Time:Time:badInputs',...
          'Unknown input type')
      end
    end
    
    function e = toEpoch(t)
      % Convert time to ISDAT epoch
      e = iso2epoch2(toUTC(t));
    end
    
    function s = toUTC(t)
      % Convert Time to UTC string
      s_tmp = encodett2000(t.tt2000);
      s = s_tmp{:};
    end
    
    function r = lt(t1,t2)
      r = (t1.tt2000 < t2.tt2000);
    end
    
    function r = gt(t1,t2)
      r = (t1.tt2000 > t2.tt2000);
    end
    
    function r = le(t1,t2)
      r = (t1.tt2000 <= t2.tt2000);
    end
    
    function r = ge(t1,t2)
      r = (t1.tt2000 >= t2.tt2000);
    end
    
    function r = ne(t1,t2)
      r = (t1.tt2000 ~= t2.tt2000);
    end
    
    function r = eq(t1,t2)
      r = (t1.tt2000 == t2.tt2000);
    end
    
    function r = minus(t1,t2)
      %Subtract an offset
      %
      % time_t2 = time_t1 - int64_offset_ns
      % time_t2 = time_t1 - double_offset_sec
      %
      % INT64 offset is treated as nano-seconds
      % DOUBLE offset is treated as seconds
      %
      % offset_sec = time_t2 - time_t1
      %
      % Compute offset in sec betweeb two times
      if isa(t2,'irf.Time')
        r = double((t1.tt2000 - t2.tt2000)*1e-9);
      elseif isa(t2,'int64')
        r = irf.Time(t1.tt2000 - t2);
      elseif isa(t2,'double')
        r = irf.Time(t1.tt2000 - int64(t2*1e9));
      else
        error('MATLAB:Time:minus:badInputs',...
          'Unknown input type')
      end
    end
    
    function r = plus(t1,offset)
      %Add an offset
      %
      % time_t2 = time_t1 + int64_offset_ns
      % time_t2 = time_t1 + double_offset_sec
      %
      % INT64 OFFSET is treated as nano-seconds
      % DOUBLE OFFSET is treated as seconds
      %
      % timearray_t2 = time_t1 + int64_offsetarray_ns
      % timearray_t2 = time_t1 + double_offsetarray_sec
      
      if isa(offset,'int64')
        if numel(offset)==1
          r = irf.Time(t1.tt2000 + offset);
        else
          r = irf.TimeArray(t1.tt2000 + offset);
        end
      elseif isa(offset,'double')
        if numel(offset)==1
          r = irf.Time(t1.tt2000 + int64(offset*1e9));
        else
          r = irf.TimeArray(t1.tt2000 + int64(offset*1e9));
        end
      else
        error('MATLAB:Time:minus:badInputs',...
          'Unknown input type')
      end
    end
    
    function display(t)
      fprintf('irf.Time : %s\n',toUTC(t));
    end
    
  end % methods
end

