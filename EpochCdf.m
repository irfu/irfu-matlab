classdef EpochCdf < GenericTimeArray
  %EpochCdf Class representing CDF epoch
  %   Detailed explanation goes here

  % ----------------------------------------------------------------------------
  % SPDX-License-Identifier: Beerware
  % "THE BEER-WARE LICENSE" (Revision 42):
  % <yuri@irfu.se> wrote this file.  As long as you retain this notice you
  % can do whatever you want with this stuff. If we meet some day, and you think
  % this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
  % ----------------------------------------------------------------------------

  methods
    function obj = EpochCdf(inp)
      if nargin==0, return, end
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
      if nargin<2, format = 0; end
      epochTmp = EpochCdf.cdfepoch2epoch(obj.epoch);
      s = epoch2iso(epochTmp,format);
    end
    function res = toEpochUnix(obj)
      res = EpochUnix(EpochCdf.cdfepoch2epoch(obj.epoch));
    end
  end
  methods (Static)
    function t_out = cdfepoch2epoch(t_in)
      t_out = (t_in-62167219200000)/1000;
    end
    function t_out = epoch2cdfepoch(t_in)
      t_out = t_in*1000 + 62167219200000;
    end
  end
end