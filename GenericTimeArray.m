classdef (Abstract) GenericTimeArray
  %GenericTimeArray Generic (Abstract) class describing a time array
  
% ----------------------------------------------------------------------------
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
    
    function res = end(obj,k,n)
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
    
    function sref = subsref(obj,idx)
      %SUBSREF handle indexing
        switch idx(1).type
          % Use the built-in subsref for dot notation
          case '.'
            sref = builtin('subsref',obj,idx);
          case '()'
            tmpEpoch = builtin('subsref',obj.epoch,idx);
            sref = feval(class(obj),tmpEpoch);
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
      lim = inp.(['to' className]);
      [idxLim,res] = tlimPrivate(obj,lim,mode);
    end
    
    % Anstract methods
    toUtc(obj)
    %toUtc  convert to UTC time string
  end
  
  methods (Access = private)
    tlimPrivate(obj,inp,mode)
  end
  
  methods (Static)
    function [ output_args ] = validate_iso_time_str( input_args )
      %verify_iso_time_str Summary of this function goes here
      %   Detailed explanation goes here
      
      output_args = true;
    end
  end
end
