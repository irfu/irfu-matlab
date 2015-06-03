classdef mms_edp_Sweep
  %MMS_EDP_SWEEP MMS EDP Sweep
  %   
  % To construct new object use:
  %
  %   obj = mms_edp_Sweep(fileName)
  %
  %   where fileName is a L1b sweep CDF file name

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

  
  properties (SetAccess = immutable)
    sweep
    varPref
    nSweeps
  end
  
  methods
    function obj = mms_edp_Sweep(fileName)
      % obj = mms_edp_Sweep(fileName)
      if nargin==0, 
        obj.nSweeps = []; obj.sweep = []; obj.varPref = ''; return
      end
      if isempty(regexp(fileName,'mms[1-4]_edp_srvy_l1b_sweeps_', 'once'))
        msg = 'File name must be mms[1-4]_edp_srvy_l1b_sweeps_*.cdf';
        irf.log('critical',msg),error(msg)
      end
      obj.sweep = dataobj(fileName);
      if isempty(obj.sweep),
        msg = 'No data loaded'; irf.log('critical',msg),error(msg)
      end
      obj.varPref = lower(obj.sweep.GlobalAttributes.Source_name{:}(1:4));
      obj.nSweeps = obj.sweep.data.([obj.varPref '_sweep_start']).nrec;
    end % CONSTRUCTOR
    
    function hout = plot(obj,h,iSweep)
      % Plot sweep
      %
      % hout = plot(obj,[h],iSweep)
      % 
      % Input: 
      %    iSweep - sweep number in the file
      %    h - axes handle [optional]
      if nargin==2, iSweep = h; h = []; 
      elseif ~isgraphics(h,'axes')
        msg = 'H bist be an axes handle'; irf.log('critical',msg),error(msg)
      end  
      if 1>iSweep || iSweep>obj.nSweeps || int8(iSweep)~=iSweep
        msg = sprintf('ISWEEP must be 1..%d',obj.nSweeps);
        irf.log('critical',msg),error(msg) %#ok<SPERR>
      end
      sweepTime = EpochTT2000([...
        obj.sweep.data.([obj.varPref '_sweep_start']).data(iSweep)...
        obj.sweep.data.([obj.varPref '_sweep_stop']).data(iSweep)]);
      [idx, epoch] = EpochTT2000(obj.sweep.data.Epoch.data).tlim(sweepTime);
      p1 = obj.sweep.data.([obj.varPref '_sweep_swept']).data(iSweep);
      % The "other probe" is the other probe in the pair 1-2, 3-4, 5-6
      if fix(p1/2)*2==p1, p2 = p1 - 1; else p2 = p1 + 1; end
      voltage1 =  obj.sweep.data.([obj.varPref '_edp_sweeps']).data(idx,p1);
      voltage2 =  obj.sweep.data.([obj.varPref '_edp_sweeps']).data(idx,p2);
      [idxBias,eBias] = ...
        EpochTT2000(obj.sweep.data.epoch_sweepsamp.data).tlim(sweepTime);
      bias1 =  obj.sweep.data.([obj.varPref '_sweep_bias1']).data(idxBias);
      bias2 =  obj.sweep.data.([obj.varPref '_sweep_bias2']).data(idxBias);
      % Find current values (biasRes) corresponding to voltages
      biasRes1 = zeros(size(voltage1))*NaN; biasRes2 = biasRes1;
      for i=1:length(idxBias)
        if i == length(idxBias)
          ii = epoch.tlim(EpochTT2000(...
            [eBias.stop().epoch sweepTime.stop().epoch]));
        else ii = epoch.tlim(eBias(i+[0 1]));
        end
        biasRes1(ii) = bias1(i);  biasRes2(ii) = bias2(i);
      end
      % Plot
      c = 'kkrrbb';
      lineStyleP1 = [c(p1) 'x-']; lineStyleP2 = [c(p2) 'o-'];
      if isempty(h),
        plot(voltage1,biasRes1,lineStyleP1,voltage2,biasRes2,lineStyleP2);
        h = gca;
      else
        plot(h,voltage1,biasRes1,lineStyleP1,voltage2,biasRes2,lineStyleP2);
      end
      set(h,'YLim',[min([bias1;  bias2])-10 max([bias1;  bias2])+10]);
      t = title(h,[obj.varPref ' ' toUtc(sweepTime(1),1)]);
      if isa(h,'handle'), set(t,'FontSize',12); end % Needed for HG2
      ylabel(h,...
        ['Bias [' getunits(obj.sweep,[obj.varPref '_sweep_bias1']) ']'])
      xlabel(h,...
        ['Voltage [' getunits(obj.sweep,[obj.varPref '_edp_sweeps']) ']'])
      legend(h,sprintf('V%d',p1),sprintf('V%d',p2))
      grid(h,'on')
      if nargout, hout = h; end
    end % PLOT
  end
  
end

