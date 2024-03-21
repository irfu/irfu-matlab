classdef (Abstract) GenericTimeArray
  %GenericTimeArray  Generic (Abstract) class describing a time array
  %
  %  Methods:
  %     disp()
  %     end()
  %     isempty()
  %     length()
  %     size()
  %     start()     First point of the time array
  %     stop()      Last point of the time array
  %     tlim()      Returns index and records within specified time interval
  %     epochUnix()
  %     tts()       Number of seconds
  %     ttns()      Number of nanoseconds
  %
  %     Static:
  %     pad_utc()
  %     validate_utc()
  %     LeapSeconds()

  % ----------------------------------------------------------------------------
  % SPDX-License-Identifier: Beerware
  % "THE BEER-WARE LICENSE" (Revision 42):
  % <yuri@irfu.se> wrote this file.  As long as you retain this notice you
  % can do whatever you want with this stuff. If we meet some day, and you think
  % this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
  % ----------------------------------------------------------------------------

  properties
    epoch
  end

  methods
    function disp(obj)
      %DISP display object
      %
      %  disp(obj)
      n = obj.length();
      fprintf('  %s object : %d records\n\n',class(obj),n);
      if n>5
        for i = 1:2, dispI(i), end
        fprintf('  ... skipped %d records ...\n',n-4);
        for i = n-1:n, dispI(i), end
      else
        for i = 1:n, dispI(i), end
      end
      fprintf('\n')
      function dispI(i)
        S.type = '()'; S.subs={i};
        tmp = subsref(obj,S);
        fprintf('  %s\n',tmp.toUtc)
      end
    end

    function res = isempty(obj)
      %ISEMPTY True for empty array.
      res = isempty(obj.epoch);
    end

    function res = length(obj)
      %LENGTH number of records in array.
      res = length(obj.epoch);
    end

    function res = size(obj)
      %SIZE size of the record's array.
      res = size(obj.epoch);
    end

    function res = minus(obj, obj1)
      %MINUS  time difference in seconds
      % T2 - T1 returns time difference between T2 and T1 in seconds
      %         T1 or T2 can be also of length 1
      if ~isa(obj1,'GenericTimeArray')
        error('irf:GenericTimeArray:minus:badInputs',...
          'inpus must be subclasses of GenericTimeArray')
      end
      if length(obj1) ~= 1 && length(obj) ~= 1 && length(obj) ~= length(obj1)
        error('irf:GenericTimeArray:minus:badInputs',...
          'inpus have different length')
      end
      res = obj.tts - obj1.tts;
    end

    function res = le(obj,obj1)
      %LE  less than of equal to
      if ~isa(obj1,'GenericTimeArray')
        error('irf:GenericTimeArray:le:badInputs',...
          'inpus must be subclasses of GenericTimeArray')
      elseif length(obj1)~=1
        error('irf:GenericTimeArray:le:badInputs',...
          'second argument must have length 1')
      end
      objTmp = convert_epoch(obj1,class(obj));
      res = obj.epoch <= objTmp.epoch;
    end

    function res = ge(obj,obj1)
      %GE  greater than of equal to
      if ~isa(obj1,'GenericTimeArray')
        error('irf:GenericTimeArray:ge:badInputs',...
          'inpus must be subclasses of GenericTimeArray')
      elseif length(obj1)~=1
        error('irf:GenericTimeArray:ge:badInputs',...
          'second argument must have length 1')
      end
      objTmp = convert_epoch(obj1,class(obj));
      res = obj.epoch >= objTmp.epoch;
    end

    function res = lt(obj,obj1)
      %LT  less than
      if ~isa(obj1,'GenericTimeArray')
        error('irf:GenericTimeArray:le:badInputs',...
          'inpus must be subclasses of GenericTimeArray')
      elseif length(obj1)~=1
        error('irf:GenericTimeArray:le:badInputs',...
          'second argument must have length 1')
      end
      objTmp = convert_epoch(obj1,class(obj));
      res = obj.epoch < objTmp.epoch;
    end

    function res = gt(obj,obj1)
      %GT  greater than
      if ~isa(obj1,'GenericTimeArray')
        error('irf:GenericTimeArray:gt:badInputs',...
          'inpus must be subclasses of GenericTimeArray')
      elseif length(obj1)~=1
        error('irf:GenericTimeArray:gt:badInputs',...
          'second argument must have length 1')
      end
      objTmp = convert_epoch(obj1,class(obj));
      res = obj.epoch > objTmp.epoch;
    end

    function res = eq(obj,obj1)
      %EQ  equal to
      if ~isa(obj1,'GenericTimeArray')
        error('irf:GenericTimeArray:gt:badInputs',...
          'inpus must be subclasses of GenericTimeArray')
      end
      len = obj.length(); len1 = obj1.length();
      if len==0 && len1==0, res = true; return
      elseif len==0 || len1==0, res = false; return
      end

      objTmp = convert_epoch(obj1,class(obj));
      if len1==1
        res = obj.epoch == objTmp.epoch;
      elseif len1==len
        res = obj.epoch == objTmp.epoch;
      else
        res = false;
      end
    end

    function res = ne(obj,obj1)
      %NE  not equal to
      if ~isa(obj1,'GenericTimeArray')
        error('irf:GenericTimeArray:gt:badInputs',...
          'inpus must be subclasses of GenericTimeArray')
      end
      res = ~(obj.eq(obj1));
    end

    function res = end(obj,k,~)
      %END last index in array
      switch k
        case 1
          res = length(obj);
        otherwise
          res = 1;
      end
    end

    function res = start(obj)
      %START  first point of the time array
      if isempty(obj)
        error('irf:GenericTimeArray:tlim:badInputs',...
          'empty input')
      end
      S.type = '()'; S.subs={1};
      res = subsref(obj,S);
    end

    function res = stop(obj)
      %STOP  last point of the time array
      if isempty(obj)
        error('irf:GenericTimeArray:tlim:badInputs',...
          'empty input')
      end
      S.type = '()'; S.subs={length(obj)};
      res = subsref(obj,S);
    end

    function [res,varargout] = sort(obj)
      %SORT   sorts time array in ascending order
      %  T2 = T1.sort;
      %  [T2, I] = T1.sort;
      %     I - index matrix
      %     T2 = T1(I);

      if isempty(obj)
        error('irf:GenericTimeArray:sort:badInputs',...
          'empty input')
      end
      nout = max(nargout,1) - 1;

      [sorted_obj,sorted_ind] = sort(obj.ttns);
      res = feval(class(obj),sorted_obj);

      if nout == 0
        varargout = {};
      elseif nout == 1
        varargout = {sorted_ind};
      else
        error('irf:GenericTimeArray:tlim:badOutputs',...
          'unknown output format')
      end
    end

    function res = horzcat(varargin)
      % tt = [t1 t2];
      % tt = [t1 t2:t3];

      objs = varargin;
      nobj = numel(objs);
      obj_ttns = int64([]);
      for iobj = 1:nobj % an obj can be a single EpochTT or a list
        nsubobj = numel(objs{iobj});
        for isubobj = 1:nsubobj
          tmp_obj = objs{iobj};
          obj_ttns = [obj_ttns; tmp_obj.ttns]; %#ok<AGROW>
        end
      end
      res = EpochTT(obj_ttns);
    end

    function res = vertcat(varargin)
      res = horzcat(varargin{:});
    end

    function [varargout] = subsref(obj,idx)
      %SUBSREF handle indexing
      switch idx(1).type
        % Use the built-in subsref for dot notation
        case '.'
          [varargout{1:nargout}] = builtin('subsref',obj,idx);
        case '()'
          tmpEpoch = builtin('subsref',obj.epoch,idx(1));
          out = feval(class(obj),tmpEpoch);
          if numel(idx) > 1
            out = builtin('subsref',out,idx(2:end));
          end
          [varargout{1:nargout}] = out;
          % No support for indexing using '{}'
        case '{}'
          error('irf:GenericTimeArray:subsref',...
            'Not a supported subscripted reference')
      end
    end

    function [idxLim,res] = tlim(obj,inp,mode)
      %TLIM   Returns index and records within specified time interval
      %
      % [IDX,Y] = tlim(X,LIM,[MODE])
      %
      % Where MODE can be:
      %        'and', 0 (default)
      %        'xor', 1
      %
      %	IDX contains indexes of rows that were returned
      %
      % Y is part of the X that is within interval
      % LIM.START <= X(:,1) < LIM.STOP for "AND" mode
      %
      % Y is part of X outside the interval for "XOR" mode:
      % X(:,1) < LIM.START & X(:,1) > LIM.STOP

      if isempty(inp)
        error('irf:GenericTimeArray:tlim:badInputs',...
          'empty limiting array')
      end
      if nargin<3, mode = 0; end
      if ischar(mode)
        switch lower(mode)
          case 'and'
            mode = 0;
          case 'xor'
            mode = 1;
          otherwise
            error('irf:GenericTimeArray:tlim:badInputs',...
              'MODE can be ''and'' (default) of ''xor''')
        end
      elseif isnumeric(mode)
        if ~any(mode==[0 1])
          error('irf:GenericTimeArray:tlim:badInputs',...
            'MODE can be 0 (''and'', default) of 1 (''xor''')
        end
      else
        error('irf:GenericTimeArray:tlim:badInputs',...
          'MODE can be 0 (''and'', default) of 1 (''xor''')
      end
      className = class(obj);
      lim = convert_epoch(inp,className);
      if nargout>1, [idxLim,res] = tlimPrivate(obj,lim,mode);
      else, idxLim = tlimPrivate(obj,lim,mode);
      end
    end

    function res = convert_epoch(obj,className)
      if isa(obj,className)
        res = obj;
      else
        res = feval(className,obj.ttns);
      end
    end

    function res = epochUnix(obj)
      if isempty(obj), res = [];
      else, res = EpochUnix.from_ttns(obj.ttns);
      end
    end
    function s = utc(obj,varargin)
      % s = utc(obj,format)
      s = GenericTimeArray.ttns2utc(obj.ttns,varargin{:});
    end
    function s = toUtc(obj,varargin)
      s = obj.utc(varargin{:});
    end

    % Abstract time operation methods
    to_ttns(epoch,index)   % static function to convert epoch to ttns
    from_ttns(obj,index)   % static function to convert ttns to epoch

    function s = ttns(obj,varargin)
      % s = ttns(obj,[index]) return index points
      s = obj.to_ttns(obj.epoch,varargin{:});
    end
    function s = tts(obj,varargin)
      % s = tts(obj,[index]) return index points
      s = double(obj.ttns(varargin{:}))/1e9;
    end
  end

  methods (Abstract)
    plus(obj,arg)
    colon(obj,varargin)
  end

  methods (Access = private)
    function [idxLim,res] = tlimPrivate(obj,inp,mode)
      % Private version of tLim, can be reloaded
      if mode==0
        idxLim = find((obj.epoch >= inp.epoch(1)) & (obj.epoch <= inp.epoch(end)));
      else
        idxLim = find((obj.epoch < inp.epoch(1)) | (obj.epoch > inp.epoch(end)));
      end
      if nargout>1
        S.type = '()'; S.subs={idxLim};
        tmpEpoch = builtin('subsref',obj.epoch,S);
        res = feval(class(obj),tmpEpoch);
      end
    end
  end

  methods (Static)
    function res = validate_utc( utc )
      %VALIDATE_UTC validate UTC string
      %   validate UTC string : yyyy-mm-ddThh:mm:ss.[mmmuuunnnZ]'

      MAX_NUM_IDX = 29;
      res = false;

      if isempty(utc), return, end
      if ~ismatrix(utc), return, end

      lUtc = size(utc,2); if all(all(utc(:,end)=='Z')), lUtc = lUtc-1; end
      idxDash = [5 8]; idxColon = [14 17]; idxT = 11; idxDot = 20;
      idxNonDigit = [idxDash idxT idxColon];
      if lUtc>=idxDot, idxNonDigit = [idxNonDigit idxDot]; end

      if lUtc < 19 || lUtc > MAX_NUM_IDX+1 || lUtc == idxDot || ...
          ~all(all(utc(:,idxT)=='T')) || ...
          ~all(all(utc(:,idxDash)=='-')) || ~all(all(utc(:,idxColon)==':')) || ...
          (lUtc>=idxDot && ~all(all(utc(:,idxDot)=='.')))
        return
      end

      idxTmp = 1:lUtc;
      idxDigit = setxor(idxNonDigit,idxTmp);
      if ~all(all(isstrprop(utc(:,idxDigit),'digit')))
        return
      end

      if  ~all( utc(:,1)=='1' | utc(:,1)=='2' ) || ... % year starts with 1 or 2
          ~all( utc(:,6)=='0' | utc(:,6)=='1' ) || ... % month starts with 0 or 1
          any( utc(:,6)=='0' & utc(:,7)=='0' )         % no month 00
        return
      end
      if any( utc(:,6)=='1' & ~(utc(:,7)=='0' | utc(:,7)=='1' | utc(:,7)=='2')) || ... % month > 12
          ~all( utc(:,9)=='0' | utc(:,9)=='1' | utc(:,9)=='2' | utc(:,9)=='3') || ...
          any( utc(:,9)=='0' & utc(:,10)=='0' ) || ... % no day 00
          any( utc(:,9)=='3' & ~(utc(:,10)=='0' | utc(:,10)=='1' ))
        return
      end
      res = true;
    end
    function res = pad_utc( utc )
      %Add missing Z and zeros to comply with yyyy-mm-ddThh:mm:ss.[mmmuuunnnZ]
      MAX_NUM_IDX = 29; idxDotIDX_DOT = 20;
      res = utc; lUtc = size(utc,2); flagAddZ = true;
      if all(all(utc(:,end)=='Z')), lUtc = lUtc - 1; flagAddZ = false; end
      if lUtc == MAX_NUM_IDX && ~flagAddZ, return, end
      if lUtc == idxDotIDX_DOT-1
        res(:,idxDotIDX_DOT) = '.'; lUtc = lUtc + 1;
      end
      if lUtc < MAX_NUM_IDX, res(:,(lUtc+1):MAX_NUM_IDX) = '0'; end % Pad with zeros
      res(:,MAX_NUM_IDX+1) = 'Z';
    end
    function res = validate_and_pad_utc(utc)
      % UTC_VALIDATE  validate UTC string and pad with zeros
      %
      % utcNew = validate_and_pad_utc(utc)
      %
      %  Validates that string is compliand with the UTC format:
      %  yyyy-mm-ddThh:mm:ss.[mmmuuunnnZ] and pads with missing zeros and Z.
      if ~GenericTimeArray.validate_utc(utc)
        error('GenericTimeArray:validate_and_pad_utc',...
          'UTC must be a string: yyyy-mm-ddThh:mm:ss.[mmmuuunnnZ]');
      end
      res = GenericTimeArray.pad_utc(utc);
    end
    function ttns = utc2ttns(utc,index)
      % Convert UTC to TT nanoseconds
      if nargin == 2
        utc = utc(index,:);
      end
      ttns = spdfparsett2000(GenericTimeArray.validate_and_pad_utc(utc));
      if ttns==int64(-9223372036854775805)
        error('GenericTimeArray:utc2ttns',...
          'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
      end
    end
    function utc = ttns2utc(ttns,format)
      % Convert TT nanoseconds to UTC
      if nargin<2 || isempty(format), format = 2; end
      utc =  char(spdfencodett2000(int64(ttns)));
      if isnumeric(format)
        switch format
          case 0, utc = utc(:,1:26);
          case 1, utc = utc(:,1:23);
          case 2
          otherwise
            error('irf:EpochUnix:toUtc:badInputs',...
              'wrong format value')
        end
        utc(:,end+1)='Z';
      elseif ischar(format)
        fmt = format;
        iyyyy = strfind(fmt,'yyyy');
        imm   = strfind(fmt,'mm');
        idd   = strfind(fmt,'dd');
        iHH   = strfind(fmt,'HH');
        iMM   = strfind(fmt,'MM');
        iSS   = strfind(fmt,'SS');
        immm  = strfind(fmt,'mmm');
        iuuu  = strfind(fmt,'uuu');
        innn  = strfind(fmt,'nnn');
        tVec9 = irf_time(ttns,'ttns>vector9');
        utc = repmat(fmt,size(ttns,1),1);
        if iyyyy, utc(:,iyyyy:iyyyy+3)= num2str(tVec9(:,1),'%04u'); end
        if imm,   utc(:,imm:imm+1)    = num2str(tVec9(:,2),'%02u'); end
        if idd,   utc(:,idd:idd+1)    = num2str(tVec9(:,3),'%02u'); end
        if iHH,   utc(:,iHH:iHH+1)    = num2str(tVec9(:,4),'%02u'); end
        if iMM,   utc(:,iMM:iMM+1)    = num2str(tVec9(:,5),'%02u'); end
        if iSS,   utc(:,iSS:iSS+1)    = num2str(tVec9(:,6),'%02u'); end
        if immm,  utc(:,immm:immm+2)  = num2str(tVec9(:,7),'%03u'); end
        if iuuu,  utc(:,iuuu:iuuu+2)  = num2str(tVec9(:,8),'%03u'); end
        if innn,  utc(:,innn:innn+2)  = num2str(tVec9(:,9),'%03u'); end
      end
    end
    function res = leap_seconds()
      % Try to read CDFLeapSeconds.txt (from NASA_cdf_patch). If not found
      % revert to hard coded values.
      % Source of both: <http://maia.usno.navy.mil/ser7/tai-utc.dat>
      % Updated: 20120401
      % Leap Seconds Table - used by CDF
      % Update it when a leap second(s) is added.
      %Year Month Day Leap Seconds      Drift
      if(exist('CDFLeapSeconds.txt', 'file'))
        n = importdata('CDFLeapSeconds.txt', ' ', 6);
        res = n.data;
        % Info about CDFLeapSeconds and time of last update.
        %info_str=['CDFLeapSeconds.txt was used.', n.textdata{2,1}(2:end)];
      else
        % File not found?
        %info_str=['No CDFLeapSeconds.txt found. Using hardcoded values,
        %updated 20160801.'];
        res = [...
          1960   1    1    1.4178180  37300.0  0.001296
          1961   1    1    1.4228180  37300.0  0.001296
          1961   8    1    1.3728180  37300.0  0.001296
          1962   1    1    1.8458580  37665.0  0.0011232
          1963  11    1    1.9458580  37665.0  0.0011232
          1964   1    1    3.2401300  38761.0  0.001296
          1964   4    1    3.3401300  38761.0  0.001296
          1964   9    1    3.4401300  38761.0  0.001296
          1965   1    1    3.5401300  38761.0  0.001296
          1965   3    1    3.6401300  38761.0  0.001296
          1965   7    1    3.7401300  38761.0  0.001296
          1965   9    1    3.8401300  38761.0  0.001296
          1966   1    1    4.3131700  39126.0  0.002592
          1968   2    1    4.2131700  39126.0  0.002592
          1972   1    1   10.0            0.0  0.0
          1972   7    1   11.0            0.0  0.0
          1973   1    1   12.0            0.0  0.0
          1974   1    1   13.0            0.0  0.0
          1975   1    1   14.0            0.0  0.0
          1976   1    1   15.0            0.0  0.0
          1977   1    1   16.0            0.0  0.0
          1978   1    1   17.0            0.0  0.0
          1979   1    1   18.0            0.0  0.0
          1980   1    1   19.0            0.0  0.0
          1981   7    1   20.0            0.0  0.0
          1982   7    1   21.0            0.0  0.0
          1983   7    1   22.0            0.0  0.0
          1985   7    1   23.0            0.0  0.0
          1988   1    1   24.0            0.0  0.0
          1990   1    1   25.0            0.0  0.0
          1991   1    1   26.0            0.0  0.0
          1992   7    1   27.0            0.0  0.0
          1993   7    1   28.0            0.0  0.0
          1994   7    1   29.0            0.0  0.0
          1996   1    1   30.0            0.0  0.0
          1997   7    1   31.0            0.0  0.0
          1999   1    1   32.0            0.0  0.0
          2006   1    1   33.0            0.0  0.0
          2009   1    1   34.0            0.0  0.0
          2012   7    1   35.0            0.0  0.0
          2015   7    1   36.0            0.0  0.0
          2017   1    1   37.0            0.0  0.0
          ];
      end
    end
  end
end
