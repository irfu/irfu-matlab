classdef EpochCdf16 < GenericTimeArray
  %EpochCdf16 Class representing CDF epoch16, seconds and picoseconds.
  %   Detailed explanation goes here

  % ----------------------------------------------------------------------------
  % SPDX-License-Identifier: Beerware
  % "THE BEER-WARE LICENSE" (Revision 42):
  % <yuri@irfu.se> wrote this file.  As long as you retain this notice you
  % can do whatever you want with this stuff. If we meet some day, and you think
  % this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
  % ----------------------------------------------------------------------------

  properties
    ps
  end

  methods
    function obj = EpochCdf16(inp,inp2)
      if nargin==0, return, end
      if nargin==2, obj.ps = inp2; end
      if isa(inp,'double')
        if min(size(inp))>1
          error('irf:EpochCdf:EpochCdf:badInputs',...
            'double input (CDF epoch) must be a column or row vector')
        end
        obj.epoch = inp(:); % column vector
      elseif isa(inp,'char')
        if GenericTimeArray.validate_utc(inp)
          epochTmp = iso2epoch(inp);
          obj.epoch = EpochCdf.epoch2cdfepoch(epochTmp);
        else
          error('irf:EpochCdf:EpochCdf:badInputs',...
            'UTC string input (char) must be in the form yyyy-mm-ddThh:mm:ss.mmmuuunnnZ')
        end
      else
        error('irf:EpochCdf:EpochCdf:badInputs',...
          'Expected inputs: double (CDF epoch) or char (yyyy-mm-ddThh:mm:ss.mmmuuunnnZ)')
      end
    end
    function s = toUtc(obj,format)
      % s = toUtc(obj,format)
      if nargin<2, format = 3; end
      d = cdflib.epoch16Breakdown([obj.epoch; obj.ps]);
      switch format
        case 0
          s = num2str(d,'%04d-%02d-%02dT%02d:%02d:%02d.%03d%03dZ');
        case 1
          s = num2str(d,'%04d-%02d-%02dT%02d:%02d:%02d.%03dZ');
        case 2
          s = num2str(d,'%04d-%02d-%02dT%02d:%02d:%02d.%03d%03d%03dZ');
        case 3
          s = num2str(d,'%04d-%02d-%02dT%02d:%02d:%02d.%03d%03d%03d%03dZ');
        otherwise
          error('wrong format')
      end
    end
    function res = toEpochUnix(obj)
      res = EpochUnix(EpochCdf.cdfepoch2epoch(obj.epoch));
    end
  end
  methods (Static)
    function t_out = cdfepoch162epoch(s,ps)
      t_out = (s-62167219200)+ps*1e-12;
    end
    function [s_out,ps_out] = epoch2cdfepoch16(t_in)
      s_out = fix(t_in)+62167219200; ps_out = mod(t_in,1)*1e12;
    end
  end
end