classdef EpochTT < GenericTimeArray
  %EpochTT Class representing TT epoch (CDF TT2000), nanoseconds since 2000.
  %
  % EpochTT(t) - initialize class, where t can be:
  %             - vector of number (double) of seconds since TT2000 epoch
  %             - vector of integer number (int64) of nanoseconds as TT2000
  %             - UTC string array

  % ----------------------------------------------------------------------------
  % SPDX-License-Identifier: Beerware
  % "THE BEER-WARE LICENSE" (Revision 42):
  % <yuri@irfu.se> wrote this file.  As long as you retain this notice you
  % can do whatever you want with this stuff. If we meet some day, and you think
  % this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
  % ----------------------------------------------------------------------------

  methods
    function obj = EpochTT(inp)
      if nargin==0, return, end
      if isa(inp,'double')
        if min(size(inp))>1
          error('irf:EpochTT:EpochTT:badInputs',...
            'input must be a column or row vector')
        end
        obj.epoch = int64(inp(:)*1e9); % column vector
      elseif isa(inp,'int64')
        if min(size(inp))>1
          error('irf:EpochTT:EpochTT:badInputs',...
            'int64 input (nanoseconds since 2000) must be a column or row vector')
        end
        obj.epoch = inp(:); % column vector
      elseif isa(inp,'char')
        if ~GenericTimeArray.validate_utc(inp)
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        end
        obj.epoch = GenericTimeArray.utc2ttns(inp);
      elseif isa(inp,'string')
        if ~GenericTimeArray.validate_utc(inp)
          error('irf:EpochUnix:EpochUnix:badInputs',...
            'UTC string input (string) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        end
        obj.epoch = GenericTimeArray.utc2ttns(inp);  
      elseif isa(inp,'GenericTimeArray')
        if isa(inp,'EpochTT')
          obj = inp;
        else
          obj = EpochTT(inp.ttns);
        end
      else
        error('irf:EpochUnix:EpochUnix:badInputs',...
          'Expected inputs: int64 (nanoseconds since 2000), double (seconds since 1970) or char (yyyy-mm-ddThh:mm:ss.mmmuuunnnZ)')
      end
    end
    function objOut = plus(obj,arg)
      if isnumeric(arg)
        if isa(arg,'double')
          inp = int64(arg*1e9);
        elseif isa(arg,'int64')
          inp = arg;
        else
          error('Input type not defined');
        end
        objOut = obj;
        objOut.epoch = obj.epoch + inp(:);
      end
    end
    function outObj = colon(obj,varargin)
      if nargin == 2 && isa(varargin{1},'EpochTT')
        tns = obj.start.ttns:int64(1e9):varargin{1}.stop.ttns;
        outObj = EpochTT(tns);
      elseif nargin == 3 && isa(varargin{2},'EpochTT') && isnumeric(varargin{1})
        tns = obj.start.ttns:int64(varargin{1}*1e9):varargin{2}.stop.ttns;
        outObj = EpochTT(tns);
      end
    end
  end

  methods (Static)
    function output = from_ttns(input,index) % for consistency with other GenericTimeArray routines
      if nargin == 1
        output = input;
      else
        output = input(index);
      end
    end
    function output = to_ttns(input,index)
      if nargin == 1
        output = input;
      else
        output = input(index);
      end
    end
  end

end