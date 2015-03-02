classdef mms_edp_Sweep
  %MMS_EDP_SWEEP MMS EDP Sweep
  %   Detailed explanation goes here
  
  properties (SetAccess = immutable)
    sweep
    mmsId
    varPref
    nSweeps
  end
  
  methods
    function obj = mms_edp_Sweep(fileName)
      if nargin==0, 
        obj.mmsId = []; obj.nSweeps = []; obj.sweep = []; obj.varPref = '';
        return
      end
      if ~regexp(fileName,'^mms[1-4]_edp_srvy_l1b_sweeps_')
        msg = 'File name must be mms[1-4]_edp_srvy_l1b_sweeps_*.cdf';
        irf.log('critical',msg),error(msg)
      end
      obj.mmsId = str2double(fileName(4));
      obj.varPref = ['mms' fileName(4)];
      obj.sweep = dataobj(fileName,'KeepTT2000');
      if isempty(obj.sweep),
        msg = 'No data loaded'; irf.log('critical',msg),error(msg)
      end
      obj.nSweeps = obj.sweep.data.([obj.varPref '_sweep_start']).nrec;
    end
    
    function hout = plot(obj,h,iSweep)
      if nargin==2, iSweep = h; h = []; 
      elseif ~isgraphics(h,'axes')
        msg = 'H bist be an axes handle'; irf.log('critical',msg),error(msg)
      end  
      if 1>iSweep || iSweep>obj.nSweeps || int8(iSweep)~=iSweep
        msg = sprintf('ISWEEP must be 1..%s',obj.nSweeps);
        irf.log('critical',msg),error(msg) %#ok<SPERR>
      end
      lim = EpochTT2000([...
        obj.sweep.data.([obj.varPref '_sweep_start']).data(iSweep)...
        obj.sweep.data.([obj.varPref '_sweep_stop']).data(iSweep)]);
      [idx, epoch] = EpochTT2000(obj.sweep.data.Epoch.data).tlim(lim);
      p1 = obj.sweep.data.([obj.varPref '_sweep_swept']).data(iSweep);
      p2 = p1 + 1;
      s1 =  obj.sweep.data.([obj.varPref '_edp_sweeps']).data(idx,p1);
      s2 =  obj.sweep.data.([obj.varPref '_edp_sweeps']).data(idx,p2);
      [idxBias,eBias] = EpochTT2000(obj.sweep.data.epoch_sweepsamp.data).tlim(lim);
      bias1 =  obj.sweep.data.([obj.varPref '_sweep_bias1']).data(idxBias);
      bias2 =  obj.sweep.data.([obj.varPref '_sweep_bias1']).data(idxBias);
      
      biasRes1 = zeros(size(s1))*NaN; biasRes2 = biasRes1;
      for i=1:length(idxBias)
        if i == length(idxBias)
        else
          ii = epoch.tlim(eBias(i+[0 1]));
          biasRes1(ii) = bias1(i);  biasRes2(ii) = bias2(i); 
        end
      end
      
      c = 'kkrrbb';
      if isempty(h),
        plot(s1,biasRes1,c(p1),s2,biasRes2,[c(p2) '--']);
        h = gca;
      else
        plot(h,s1,biasRes1,c(p1),s2,biasRes2,[c(p2) '--']);
      end
      yl = h.YLim; h.YLim = yl*1.05;
      t=title(h,[toUtc(lim(1),1) ' -- ' toUtc(lim(2),1)]); t.FontSize=8;
      ylabel(h,...
        ['Bias [' getunits(obj.sweep,[obj.varPref '_sweep_bias1']) ']'])
      xlabel(h,...
        ['Voltage [' getunits(obj.sweep,[obj.varPref '_edp_sweeps']) ']'])
      legend(h,sprintf('V%d',p1),sprintf('V%d',p2))
      grid(h,'on')
      if nargout, hout = h; end
    end
  end
  
end

