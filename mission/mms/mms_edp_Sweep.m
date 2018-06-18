classdef mms_edp_Sweep < handle
  %MMS_EDP_SWEEP MMS EDP Sweep
  %   
  % To construct new object use:
  %
  %   obj = mms_edp_Sweep(fileName[,printStatus])
  %
  %   where fileName is a L1b sweep CDF file name
  %   and printStatus indicates whether to print (1) or not print (0) status for each sweep

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

  
  properties (SetAccess = immutable)
    sweep
    scId
    nSweeps
    pTable
  end
  properties (SetAccess = private)
    p1
    p2
    p3
    p4
    p5
    p6
  end
  
  methods
    
    function obj = mms_edp_Sweep(fileName,printStatus)
      if nargin==0
        obj.nSweeps = []; obj.sweep = []; obj.scId = ''; return
      end
      if nargin==1
        printStatus = 0;
      end
      if isempty(regexp(fileName,'mms[1-4]_edp_srvy_l1b_sweeps_', 'once'))
        msg = 'File name must be mms[1-4]_edp_srvy_l1b_sweeps_*.cdf';
        irf.log('critical',msg),error(msg)
      end
      obj.sweep = dataobj(fileName);
      if isempty(obj.sweep)
        msg = 'No data loaded'; irf.log('critical',msg),error(msg)
      end
      obj.scId = lower(obj.sweep.GlobalAttributes.Source_name{:}(1:4));
      obj.nSweeps = obj.sweep.data.([obj.scId '_sweep_start']).nrec;
      obj.pTable = [];
      obj.p1 = []; obj.p2 = []; obj.p3 = []; obj.p4 = []; obj.p5 = []; obj.p6 = [];
      for iSweep = 1:obj.nSweeps
      %for iSweep = 1:5
        [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
          eVolt, eBias, v01, v02]...
          = obj.getSweep(iSweep);
        % Determine type of sweep:
        % ++ = up only, -- = down only, +- = up-down, -+ = down-up, 00 = wiggly
        bias_no_NaN = biasRes1(isfinite(biasRes1));
        n = length(bias_no_NaN); n2 = round(n/2);
        if n>1
          p = polyfit((1:n)',bias_no_NaN(1:n),1); slope = p(1);
          p = polyfit((1:n2)',bias_no_NaN(1:n2),1); slope1 = p(1);
          p = polyfit((n2+1:n)',bias_no_NaN(n2+1:n),1); slope2 = p(1);
          if slope > 0.1
            type = '++';
          elseif slope < -0.1
            type = '--';
          elseif slope1 > 0.1 && slope2 < -0.1
            type = '+-';
          elseif slope1 < -0.1 && slope2 > 0.1
            type = '-+';
          else
            type = '00';
          end
        else
          type = '**';
        end
        if printStatus
          if ~isempty(biasRes1(isnan(biasRes1)))
            nanstr = [num2str(length(biasRes1(isnan(biasRes1)))) ' NaNs '];
          else
            nanstr = '';
          end
          if length(voltage1) > 1
            millis = [
              num2str(1e3*(eVolt(1)-sweepTime(1)),'%6.3f') ' ' ...
              num2str(1e3*(eVolt(2)-sweepTime(1)),'%6.3f') ' ' ...
              num2str(1e3*(eBias(1)-sweepTime(1)),'%6.3f') ' ' ...
              num2str(1e3*(eBias(2)-sweepTime(1)),'%6.3f') ' ' ...
              num2str(round(1/(eVolt(2)-eVolt(1)))) ' ' ...
              num2str(round(1/(eBias(2)-eBias(1)))) ' ' ];
          else
            millis = '';
          end
          disp([num2str(iSweep,'%3.2d') ': ' ...
            sweepTime(1).toUtc ' ' ...
            num2str((sweepTime(2)-sweepTime(1)),'%3.2f') 's ' ...
            num2str(prb1) ' ' num2str(prb2) ' ' ...
            'tbl=' num2str(obj.sweep.data.([obj.scId '_sweep_table']).data(iSweep)) ' ' ...
            'dwl=' num2str(obj.sweep.data.([obj.scId '_sweep_dwell']).data(iSweep)) ' ' ...
            'stp=' num2str(obj.sweep.data.([obj.scId '_sweep_steps']).data(iSweep)) ' ' ...
            'len=' num2str(length(voltage1),'%4.4d') ' ' ...
            'rt=' num2str(length(voltage1)/...
            (double(obj.sweep.data.([obj.scId '_sweep_dwell']).data(iSweep))*...
            double(obj.sweep.data.([obj.scId '_sweep_steps']).data(iSweep)))) ' ' ...
            'typ=' type ' ' ...
            '[' num2str(round(min(biasRes1(isfinite(biasRes1))))) ','...
            num2str(round(max(biasRes1(isfinite(biasRes1))))) '] ' ...
            nanstr ...
            millis ...
            ]);
          %if isnan(biasRes1(1))
          %  [biasRes1(1:16), voltage1(1:16)], [biasRes2(1:16), voltage2(1:16)]
          %end
        end
        if isempty(voltage1)
          disp(['*** Warning, empty sweep ',num2str(iSweep),' of ',num2str(obj.nSweeps)]);
        end
        tmp.time = sweepTime;
        tmp.iPh = NaN;
        tmp.impedance = NaN;
        tmp.phase = NaN;
        tmp.type = type;
        switch prb1
          case 1, p_1='p1'; p_2='p2';
          case 2, p_1='p2'; p_2='p1';
          case 3, p_1='p3'; p_2='p4';
          case 4, p_1='p4'; p_2='p3';
          case 5, p_1='p5'; p_2='p6';
          case 6, p_1='p6'; p_2='p5';
          otherwise, error('Bad probe number')
        end
        obj.pTable = [ obj.pTable [ prb1; length(obj.(p_1)) + 1 ] ];
        obj.(p_1) = [ obj.(p_1) tmp ];
        obj.(p_2) = [ obj.(p_2) tmp ];
      end
    end % CONSTRUCTOR
    
    function analyze(obj,iSweep,sps)
      % Analyze sweep
      %
      % This function analyzes a sweep and determines iPh, impedance,
      % and phase for one sweep in a sweep object obj
      %
      % analyze(obj,iSweep[,sps])
      %
      % Input: iSweep - sweep number in the file
      %        sps is a sunpulse structure obtained from sunpulse_from_hk101()
      %
      % If sps is absent no phase is computed
      narginchk(2,3);
      if(nargin==2), doPhase=false; else, doPhase=true; end
      [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
          eVolt, eBias, v01, v02]...
          = getSweep(obj,iSweep);
      if isempty(voltage1)
        disp(['*** Warning, not analyzed empty sweep ',num2str(iSweep),...
          ' of ',num2str(obj.nSweeps)]);
        return
      end
      if size(voltage1,1) <2
        disp(['*** Warning, not analyzed almost empty sweep ',num2str(iSweep),...
          ' of ',num2str(obj.nSweeps)]);
        return
      end
      % Determine type of sweep:
      % ++ = up only, -- = down only, +- = up-down, -+ = down-up, 00 = wiggly
      type = obj.(['p' num2str(obj.pTable(1,iSweep))])(obj.pTable(2,iSweep)).type;
      if strcmp(type,'00')
        disp(['*** Warning, not analyzed type ',type,' sweep ',num2str(iSweep),...
          ' of ',num2str(obj.nSweeps)]);
        return
      end
%       tmp1.iPh = NaN; tmp2.iPh = NaN;
%       tmp1.impedance = NaN; tmp2.impedance = NaN;
%       tmp1.phase = NaN; tmp2.phase = NaN;
      movedNan = 0;
%      while isnan(biasRes1(1))
%        biasRes1(1:end-1) = biasRes1(2:end);
%        biasRes2(1:end-1) = biasRes2(2:end);
%        movedNan = movedNan + 1;
%      end
      if movedNan > 0
        disp(['*** Warning, ' num2str(movedNan) ' NaNs removed by downshift of biasRes in sweep '...
          num2str(iSweep)])
      end
      if (prb1 == 1 || prb1 == 3) && ...
          round(min(biasRes1)) == -250 && round(max(biasRes1)) == 100
        biasRes1 = (biasRes1+250)/(100+250)*(50+460) - 460;
        biasRes2 = (biasRes2+250)/(100+250)*(50+460) - 460;
        disp('*** Warning, biasRes tables changed from [-250,100] to [-460,50]');
      end
      % Compute photoemission iPh and impedance dV/dI
      tmp1 = compute_IPh_Impedance(voltage1, v01, biasRes1);
      % do the same for the other probe
      tmp2 = compute_IPh_Impedance(voltage2, v02, biasRes2);
      if doPhase
        if strcmp(type,'+-') || strcmp(type,'-+') || strcmp(type,'00')
          disp(['*** Warning, phase not computed for type ',type,' sweep ',num2str(iSweep),...
            ' of ',num2str(obj.nSweeps)]);
          ph_tmp1 = NaN;
          ph_tmp2 = NaN;
        else
          t_iph = sweepTime.start.epoch ...
            + (-tmp1.iPh - biasRes1(1)) / (biasRes1(end)-biasRes1(1)) ...
            * (sweepTime.stop.epoch-sweepTime.start.epoch);
          if isfield(sps,'zphase')
            phase = mms_defatt_phase(sps,t_iph);
            ph_tmp1 = phase.data;
          else
            [ph_tmp1, ~] = mms_sdp_phase_2(sps,t_iph);
          end
          t_iph = sweepTime.start.epoch ...
            + (-tmp2.iPh - biasRes2(1)) / (biasRes2(end)-biasRes2(1)) ...
            * (sweepTime.stop.epoch-sweepTime.start.epoch);
          if isfield(sps,'zphase')
            phase = mms_defatt_phase(sps,t_iph);
            ph_tmp2 = phase.data;
          else
            [ph_tmp2, ~] = mms_sdp_phase_2(sps,t_iph);
          end
        end
      else
        ph_tmp1 = NaN;
        ph_tmp2 = NaN;
      end
      switch obj.pTable(1,iSweep)
        case 1, p_1='p1'; p_2='p2'; angle1=30; angle2=210;
        case 2, p_1='p2'; p_2='p1'; angle1=210; angle2=30;
        case 3, p_1='p3'; p_2='p4'; angle1=120; angle2=300;
        case 4, p_1='p4'; p_2='p3'; angle1=300; angle2=120;
        case 5, p_1='p5'; p_2='p6'; angle1=0; angle2=0;
        case 6, p_1='p6'; p_2='p5'; angle1=0; angle2=0;
        otherwise, error('Bad probe number')
      end
      obj.(p_1)(obj.pTable(2,iSweep)).iPh = tmp1.iPh;
      obj.(p_2)(obj.pTable(2,iSweep)).iPh = tmp2.iPh;
      obj.(p_1)(obj.pTable(2,iSweep)).impedance = tmp1.impedance;
      obj.(p_2)(obj.pTable(2,iSweep)).impedance = tmp2.impedance;
      obj.(p_1)(obj.pTable(2,iSweep)).phase = mod(ph_tmp1+angle1+90,360)-90;
      obj.(p_2)(obj.pTable(2,iSweep)).phase = mod(ph_tmp2+angle2+90,360)-90;
      return

      % Help function
      function tmp = compute_IPh_Impedance(voltage, v0, biasRes)
        % Compute photoemission iPh and impedance dV/dI
        narginchk(3,3);
         % Default output
        tmp = struct('iPh', NaN, 'impedance', NaN, 'phase', NaN);
        % Second version, to be refined...
        % Find maximum negative voltage
        vmin = min(voltage);
        % Find v0 for first point
        %v0 = voltage1(1);
        %v0 = v01;
        % Find -Iph for the average voltage (2*V0+Vmin)/3
        % skipping the first 8 bias values (settling time)...
        [~,indIph] = min(abs(voltage(9:end)-(2*v0+vmin)/3));
        tmp.iPh = -biasRes(8+indIph);
        % Find dVdI = slope of curve between -0.8*Iph and -0.7*Iph
        [~,ind80] = min(abs(biasRes-(-0.8*tmp.iPh)));
        [~,ind70] = min(abs(biasRes-(-0.7*tmp.iPh)));
        %disp(['1 ind80 ',num2str(ind80),' ind70 ',num2str(ind70),' n ',num2str(ind70-ind80+1)]);
        if ind70-ind80 > 1
          p = polyfit(biasRes(ind80:ind70),voltage(ind80:ind70),1);
          tmp.impedance = 1000*p(1);
        end
        if ind70-ind80 < 1
          p = polyfit(biasRes(ind70:ind80),voltage(ind70:ind80),1);
          tmp.impedance = 1000*p(1);
        end
      end

    end % analyze

    function analyze_all(obj,sps,printStatus)
      % Analyze sweeps
      %
      % This function analyzes sweeps and determines iPh, impedance,
      % and phase for each sweep in a sweep object obj
      %
      % analyze_all(obj[,sps[,printStatus]])
      %
      % Input: sps is a sunpulse structure obtained from sunpulse_from_hk101()
      % printStatus indicates whether to print (1) or not print (0) status for each sweep
      %
      % If sps is absent no phase is computed
      % If printStatus is absent, no status is printed for each sweep
      %
      if nargin == 1
        for iSweep = 1:obj.nSweeps
          analyze(obj,iSweep)
        end
      elseif nargin == 2
        for iSweep = 1:obj.nSweeps
          analyze(obj,iSweep,sps)
        end
      elseif nargin == 3
        for iSweep = 1:obj.nSweeps
          if printStatus, disp(['Analyzing sweep ',num2str(iSweep)]), end
          analyze(obj,iSweep,sps)
        end
      end
      return
    end % analyze_all
    
    function [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        eVolt, eBias, v01, v02]...
        = debug(obj,iSweep)
      % Debug sweep
      %
      % [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
      %   eVolt, eBias, v01, v02]...
      %   = debug(obj,iSweep)
      %
      % Input: iSweep - sweep number in the file
      %
      % Output:
      %  sweepTime
      %  prb1 - swept probe
      %  prb2 - other probe
      %  voltage1 - voltage values of swept probe
      %  biasRes1 - bias values of swept probe
      %  voltage2 - voltage values of other probe
      %  biasRes2 - bias values of other probe
      %  eVolt - epoch of voltage values
      %  eBias - epoch of bias values
      %  v01 - last voltage value of swept probe before sweep
      %  v02 - last voltage value of other probe before sweep
      %
      [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        eVolt, eBias, v01, v02]...
        = getSweep(obj,iSweep);
      return
    end
      
    function list(obj,iSweep)
      % List sweep
      %
      % list(obj,iSweep)
      %
      % Input: iSweep - sweep number in the file
      %
      [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        eVolt, eBias, v01, v02]...
        = getSweep(obj,iSweep);
      disp(['Sweep start: ' irf_time(sweepTime(1).epoch,'ttns>utc')])
      disp(['Sweep end:   ' irf_time(sweepTime(2).epoch,'ttns>utc')])
      time = sweepTime(1).epochUnix + (0:length(biasRes1)-1)' * ...
        (sweepTime(2).epochUnix-sweepTime(1).epochUnix) / ...
        (length(biasRes1)-0);
      time_1 = sweepTime(1).epochUnix - ...
        (sweepTime(2).epochUnix-sweepTime(1).epochUnix) / ...
        (length(biasRes1)-0);
      disp([ irf_time(time_1,'utc') ...
        ' ' '---' ' ' '---' ...
        ' ' num2str(v01) ' ' num2str(v02) ])
      for i = 1:length(biasRes1)
        disp([ irf_time(time(i),'utc') ...
          ' ' num2str(biasRes1(i)) ' ' num2str(biasRes2(i)) ...
          ' ' num2str(voltage1(i)) ' ' num2str(voltage2(i)) ...
          ' ' num2str((eBias(fix((i-1)*length(eBias)/length(biasRes1))+1)-sweepTime(1))*1e3) ...
          ' ' num2str((eVolt(i)-sweepTime(1))*1e3) ])
      end
% for debugging 2015-06-02
      %eBias
      %eVolt
      return
    end % list

    function hout = plot(obj,h,iSweep)
      % Plot sweep
      %
      % Plots a sweep as I(V)
      %
      % hout = plot(obj,[h],iSweep)
      % 
      % Input: 
      %    iSweep - sweep number in the file
      %    h - axes handle [optional]
      %
      if nargin==2, iSweep = h; h = []; 
      elseif ~isgraphics(h,'axes')
        msg = 'H bist be an axes handle'; irf.log('critical',msg),error(msg)
      end  
      if 1>iSweep || iSweep>obj.nSweeps || round(iSweep)~=iSweep
        msg = sprintf('ISWEEP must be 1..%d',obj.nSweeps);
        irf.log('critical',msg),error(msg) %#ok<SPERR>
      end
      [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        eVolt, eBias, v01, v02]...
        = getSweep(obj,iSweep);
      if isempty(voltage1)
        disp(['*** Warning, not plotted empty sweep ',num2str(iSweep),...
          ' of ',num2str(obj.nSweeps)]);
        return
      end
      dwl = obj.sweep.data.([obj.scId '_sweep_dwell']).data(iSweep);
      stp = obj.sweep.data.([obj.scId '_sweep_steps']).data(iSweep);
      len = length(voltage1);
      % Determine type of sweep:
      % ++ = up only, -- = down only, +- = up-down, -+ = down-up, 00 = wiggly
      if iSweep <= length(obj.pTable)
        type = obj.(['p' num2str(obj.pTable(1,iSweep))])(obj.pTable(2,iSweep)).type;
      else
        type = '  ';
      end
      movedNan = 0;
%      while isnan(biasRes1(1))
%        biasRes1(1:end-1) = biasRes1(2:end);
%        biasRes2(1:end-1) = biasRes2(2:end);
%        movedNan = movedNan + 1;
%      end
      if movedNan > 0
        disp(['*** Warning, ' num2str(movedNan) ' NaNs removed by downshift of biasRes in sweep '...
          num2str(iSweep)])
      end
      if (prb1 == 1 || prb1 == 3) && ...
          round(min(biasRes1)) == -250 && round(max(biasRes1)) == 100
        biasRes1 = (biasRes1+250)/(100+250)*(50+460) - 460;
        biasRes2 = (biasRes2+250)/(100+250)*(50+460) - 460;
        disp('*** Warning, biasRes tables changed from [-250,100] to [-460,50]')
      end
% Plot
      if isempty(h), clf, else, clf(h), end
      c = 'krgbmc';
      lineStyleP1 = [c(prb1) '.-']; lineStyleP2 = [c(prb2) '.-'];
      if strcmp(type,'++') || strcmp(type,'--') || strcmp(type,'00')
        if isempty(h)
          plot(voltage1,biasRes1,lineStyleP1,...
            voltage2,biasRes2,lineStyleP2,voltage1(1),NaN);
          h = gca;
        else
          plot(h,voltage1,biasRes1,lineStyleP1,...
            voltage2,biasRes2,lineStyleP2,voltage1(1),NaN);
        end
      else
        n=length(biasRes1(isfinite(biasRes1))); n2=round(n/2);
        lineStyleP3 = [c(prb1) '.--']; lineStyleP4 = [c(prb2) '.--'];
        if isempty(h)
          plot(voltage1(n2:n),biasRes1(n2:n),lineStyleP3,...
            voltage2(n2:n),biasRes2(n2:n),lineStyleP4);
          hold(gca,'on')
          plot(voltage1(1:n2),biasRes1(1:n2),lineStyleP1,...
            voltage2(1:n2),biasRes2(1:n2),lineStyleP2,voltage1(1),NaN);
          hold(gca,'off')
          h = gca;
        else
          plot(h,voltage1(n2:n),biasRes1(n2:n),lineStyleP3,...
            voltage2(n2:n),biasRes2(n2:n),lineStyleP4);
          hold(gca,'on')
          plot(h,voltage1(1:n2),biasRes1(1:n2),lineStyleP1,...
            voltage2(1:n2),biasRes2(1:n2),lineStyleP2,voltage1(1),NaN);
          hold(gca,'off')
        end
      end
      ylim = [min([biasRes1; biasRes2])-10 max([biasRes1; biasRes2])+10];
      set(h,'YLim',ylim);
      t = title(h,[obj.scId ' ' toUtc(sweepTime(1),1) ' sweep ' num2str(iSweep)...
        ', len=' num2str(len) ' stp=' num2str(stp) ' dwl=' num2str(dwl)...
        ' rate=' num2str(128*len/double(stp)/double(dwl))]);
      %t = title(h,[obj.scId ' ' toUtc(sweepTime(1),1) ' sweep ' num2str(iSweep)]);
      if isa(h,'handle'), set(t,'FontSize',12); end % Needed for HG2
      ylabel(h,...
        ['Bias [' getunits(obj.sweep,[obj.scId '_sweep_bias1']) ']'])
      xlabel(h,...
        ['Voltage [' getunits(obj.sweep,[obj.scId '_edp_sweeps']) ']'])
% get values for Iph and dVdI for the relevant probes
      if iSweep <= length(obj.pTable)
        switch obj.pTable(1,iSweep)
          case 1, p_1='p1'; p_2='p2';
          case 2, p_1='p2'; p_2='p1';
          case 3, p_1='p3'; p_2='p4';
          case 4, p_1='p4'; p_2='p3';
          case 5, p_1='p5'; p_2='p6';
          case 6, p_1='p6'; p_2='p5';
        end
        Iph1 = obj.(p_1)(obj.pTable(2,iSweep)).iPh;
        Iph2 = obj.(p_2)(obj.pTable(2,iSweep)).iPh;
        dVdI1 = obj.(p_1)(obj.pTable(2,iSweep)).impedance;
        dVdI2 = obj.(p_2)(obj.pTable(2,iSweep)).impedance;
        pha1 = obj.(p_1)(obj.pTable(2,iSweep)).phase;
        pha2 = obj.(p_2)(obj.pTable(2,iSweep)).phase;
      else
        Iph1 = NaN; Iph2 = NaN;
        dVdI1 = NaN; dVdI2 = NaN;
        pha1 = NaN; pha2 = NaN;
      end
      % Output for debugging
      %[double(prb1) double(prb2) Iph1 Iph2 dVdI1 dVdI2]
      legend(h,...
        sprintf('V%d  Iph %0.1f nA,  %0.2f MOhm @ %0.1f nA',prb1,Iph1,dVdI1,0.75*Iph1),...
        sprintf('V%d  Iph %0.1f nA,  %0.2f MOhm @ %0.1f nA',prb2,Iph2,dVdI2,0.75*Iph2),...
        sprintf('Phase%d  %0.1f,  Phase%d  %0.1f degrees Type %s',prb1,pha1,prb2,pha2,type),...
        'Location','NorthWest')
      grid(h,'on')
      if nargout, hout = h; end
    end % PLOT
    
    function hout = plot_time(obj,h,iSweep,sps,evfile)
      % Plot sweep
      %
      % Plots a sweep as I(t), V(t) and dVdI(t) [ and Phase(t) ]
      %
      % hout = plot2(obj,[h],iSweep[,sps])
      % 
      % Input: 
      %    iSweep - sweep number in the file
      %    h - axes handle [optional]
      %    sps - a sunpulse structure obtained from sunpulse_from_hk101()
      %    if sps is not given, no phase is plotted
      %
      if nargin==4, evfile = sps; plotev = true; else, plotev = false; end
      if nargin>=3, sps = iSweep; iSweep = h; h = []; doPhase = true;
      elseif nargin == 2, iSweep = h; h = []; doPhase = false;
      elseif ~isgraphics(h,'axes')
        msg = 'H must be an axes handle'; irf.log('critical',msg),error(msg)
      end  
      if 1>iSweep || iSweep>obj.nSweeps || round(iSweep)~=iSweep
        msg = sprintf('ISWEEP must be 1..%d',obj.nSweeps);
        irf.log('critical',msg),error(msg) %#ok<SPERR>
      end
      [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        eVolt, eBias, v01, v02]...
        = getSweep(obj,iSweep);
      if isempty(voltage1)
        disp(['*** Warning, not plotted empty sweep ',num2str(iSweep),...
          ' of ',num2str(obj.nSweeps)]);
        return
      end
      dwl = obj.sweep.data.([obj.scId '_sweep_dwell']).data(iSweep);
      stp = obj.sweep.data.([obj.scId '_sweep_steps']).data(iSweep);
      len = length(voltage1);
      ndwl = len/stp;
      impedance1 = NaN(size(biasRes1)); impedance2 = impedance1;
      ind_impedance = ndwl/2:ndwl:length(biasRes1)-ndwl/2;
      for i=ind_impedance
        i1=i+0*ndwl-1-(mod(i-1,ndwl));
        i2=i+1*ndwl-1-(mod(i-1,ndwl));
        i3=i+2*ndwl-1-(mod(i-1,ndwl));
        if i == ndwl/2
          impedance1(i) = (voltage1(i2)-voltage1(i3))/(biasRes1(i2)-biasRes1(i3)) * 1000;
          impedance2(i) = (voltage2(i2)-voltage2(i3))/(biasRes2(i2)-biasRes2(i3)) * 1000;
        elseif i == length(biasRes1)-ndwl/2
          impedance1(i) = (voltage1(i1)-voltage1(i2))/(biasRes1(i1)-biasRes1(i2)) * 1000;
          impedance2(i) = (voltage2(i1)-voltage2(i2))/(biasRes2(i1)-biasRes2(i2)) * 1000;
        else
          impedance1(i) = ((voltage1(i1)-voltage1(i2))/(biasRes1(i1)-biasRes1(i2))+...
            (voltage1(i2)-voltage1(i3))/(biasRes1(i2)-biasRes1(i3))) / 2 * 1000;
          impedance2(i) = ((voltage2(i1)-voltage2(i2))/(biasRes2(i1)-biasRes2(i2))+...
            (voltage2(i2)-voltage2(i3))/(biasRes2(i2)-biasRes2(i3))) / 2 * 1000;
        end
      end
      time = sweepTime(1).epochUnix + (0:length(biasRes1)-1)' * ...
        (sweepTime(2).epochUnix-sweepTime(1).epochUnix) / ...
        (length(biasRes1)-0);
      if doPhase
        time_ph = sweepTime.epoch(1) + int64(0:length(biasRes1)-1)' * ...
          (sweepTime.epoch(2)-sweepTime.epoch(1)) / ...
          int64(length(biasRes1)-0);
          if isfield(sps,'zphase')
            phase = mms_defatt_phase(sps,time_ph);
            ph_tmp = phase.data;
          else
            [ph_tmp, ~] = mms_sdp_phase_2(sps,time_ph);
          end
        switch prb1
          case 1, p_1='p1'; p_2='p2'; angle1=30; angle2=210; %#ok<NASGU>
          case 2, p_1='p2'; p_2='p1'; angle1=210; angle2=30; %#ok<NASGU>
          case 3, p_1='p3'; p_2='p4'; angle1=120; angle2=300; %#ok<NASGU>
          case 4, p_1='p4'; p_2='p3'; angle1=300; angle2=120; %#ok<NASGU>
          case 5, p_1='p5'; p_2='p6'; angle1=0; angle2=0; %#ok<NASGU>
          case 6, p_1='p6'; p_2='p5'; angle1=0; angle2=0; %#ok<NASGU>
          otherwise, error('Bad probe number')
        end
        phase1 = mod(ph_tmp+angle1+90,360)-90;
        %phase2 = mod(ph_tmp+angle2+90,360)-90;
      end
      clf
      if doPhase
        if plotev
          evobj=dataobj(evfile);
          efield = getmat(evobj,[obj.scId,'_edp_dce_sensor']);
          vfield = getmat(evobj,[obj.scId,'_edp_dcv_sensor']);
          ind = find(vfield(:,1)>=2*time(1)-time(2) & ...
            vfield(:,1)<= 2*time(end)-time(length(time)-1));
        switch prb1
            case 1
                v12 = [efield(ind,1) efield(ind,2)];
                v1 = [efield(ind,1) vfield(ind,2) vfield(ind,2)-efield(ind,2)*0.12];
            case 3
                v12 = [efield(ind,1) efield(ind,3)];
                v1 = [efield(ind,1) vfield(ind,3) vfield(ind,3)-efield(ind,3)*0.12];
            case 5
                v12 = [efield(ind,1) efield(ind,4)];
                v1 = [efield(ind,1) vfield(ind,4) vfield(ind,4)-efield(ind,4)*0.12];
            otherwise, error('Bad probe number')
        end
          h1=irf_plot({[[2*time(1)-time(2); time] [NaN; biasRes1]],...
            [[2*time(1)-time(2); time] double([v01; voltage1]) double([v02; voltage2])], ...
            [v1(:,1) v1(:,2) v1(:,3)],...
            [v12(:,1) v12(:,2)],...
            [time(ind_impedance)...
            impedance1(ind_impedance)...
            impedance2(ind_impedance)],...
            [time phase1]},'.-');
        else
          h1=irf_plot({[[2*time(1)-time(2); time] [NaN; biasRes1]],...
            [[2*time(1)-time(2); time] double([v01; voltage1]) double([v02; voltage2])], ...
            [time(ind_impedance)...
            impedance1(ind_impedance)...
            impedance2(ind_impedance)],...
            [time phase1]},'.-');
        end
      else
        h1=irf_plot({[[2*time(1)-time(2); time] [NaN; biasRes1]],...
          [[2*time(1)-time(2); time] double([v01; voltage1]) double([v02; voltage2])], ...
          [time(ind_impedance)...
          impedance1(ind_impedance)...
          impedance2(ind_impedance)]...
          },'.-');
      end
      title(h1(1),sprintf('%s %s sweep %i, len=%d stp=%d dwl=%d rate=%d', ...
        obj.scId, toUtc(sweepTime(1),1), iSweep, len, stp, dwl, ...
        128*len/double(stp)/double(dwl)));
      ylabel(h1(1),{'Bias',['[' getunits(obj.sweep,[obj.scId '_sweep_bias1']) ']']});
      ylabel(h1(2),{'Voltage',['[' getunits(obj.sweep,[obj.scId '_edp_sweeps']) ']']});
      legend(h1(1),sprintf('p%d',prb1))
      legend(h1(2),['p',num2str(prb1)],['p',num2str(prb2)])
      imp1 = impedance1(ind_impedance); medi1 = median(imp1(isfinite(imp1)));
      imp2 = impedance2(ind_impedance); medi2 = median(imp2(isfinite(imp2)));
      n1 = obj.pTable(1,iSweep);
      n2 = n1 - 1 + 2 * mod(n1,2);
      obj.(['p' num2str(n1)])(obj.pTable(2,iSweep)).impedance = medi1;
      obj.(['p' num2str(n2)])(obj.pTable(2,iSweep)).impedance = medi2;
      if doPhase
        if plotev
          ylabel(h1(3),{'V',['[' getunits(obj.sweep,[obj.scId '_edp_sweeps']) ']']});
          ylabel(h1(4),{'E',['[' 'mV/m' ']']});
          ylabel(h1(5),{'dVdI','[MOhm]'})
          ylabel(h1(6),{'Phase','[deg]'})
          legend(h1(3),['p',num2str(prb1)],['p',num2str(prb2)])
          legend(h1(4),['p',num2str(prb1),num2str(prb2)])
          legend(h1(5),['p',num2str(prb1),' median=',num2str(medi1,'%4.1f')],...
              ['p',num2str(prb2),' median=',num2str(medi2,'%4.1f')])
          legend(h1(6),sprintf('p%d',prb1))
        else
          ylabel(h1(3),{'dVdI','[MOhm]'})
          ylabel(h1(4),{'Phase','[deg]'})
          legend(h1(3),['p',num2str(prb1),' median=',num2str(medi1,'%4.1f')],...
              ['p',num2str(prb2),' median=',num2str(medi2,'%4.1f')])
          legend(h1(4),sprintf('p%d',prb1))
        end
      else
        ylabel(h1(3),{'dVdI','[MOhm]'})
        legend(h1(3),['p',num2str(prb1),' median=',num2str(medi1,'%4.1f')],...
            ['p',num2str(prb2),' median=',num2str(medi2,'%4.1f')])
      end
      if nargout, hout = h1; end
    end % PLOT_TIME
    
  end % methods
  
  methods(Access=private)
    
    function [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        Epoch, eBias, v01, v02]...
        = getSweep(obj,iSweep)
      % get probe, voltage and current values for one sweep
      sweepTime = EpochTT([...
        obj.sweep.data.([obj.scId '_sweep_start']).data(iSweep)+0e9/128 ...
        obj.sweep.data.([obj.scId '_sweep_stop']).data(iSweep)+0e9/128]);
      [idx, Epoch] = tlim(EpochTT(obj.sweep.data.Epoch.data), sweepTime);
      prb1 = obj.sweep.data.([obj.scId '_sweep_swept']).data(iSweep);
      % The "other probe" is the other probe in the pair 1-2, 3-4, 5-6
      if fix(prb1/2)*2==prb1, prb2 = prb1 - 1; else, prb2 = prb1 + 1; end
      voltage1 =  obj.sweep.data.([obj.scId '_edp_sweeps']).data(idx,prb1);
      voltage2 =  obj.sweep.data.([obj.scId '_edp_sweeps']).data(idx,prb2);
      if ~isempty(voltage1)
        v01 = obj.sweep.data.([obj.scId '_edp_sweeps']).data(idx(1)-1,prb1);
        v02 = obj.sweep.data.([obj.scId '_edp_sweeps']).data(idx(1)-1,prb2);
      else
        v01 = NaN; v02 = NaN;
      end
% for debugging 2015-06-02
%      sweepTime = EpochTT([...
%        obj.sweep.data.([obj.scId '_sweep_start']).data(iSweep)+0e9/128-1e5 ...
%        obj.sweep.data.([obj.scId '_sweep_stop']).data(iSweep)+0e9/128]);
      [idxBias,eBias] = ...
        tlim(EpochTT(obj.sweep.data.epoch_sweepsamp.data+0e9/128), sweepTime);
      %[idxBias,eBias] = ...
      %  tlim(EpochTT(obj.sweep.data.epoch_sweepsamp.data+0e9/128), sweepTime);
% for debugging 2015-06-02
%      sweepTime = EpochTT([...
%        obj.sweep.data.([obj.scId '_sweep_start']).data(iSweep)+0e9/128-0e5 ...
%        obj.sweep.data.([obj.scId '_sweep_stop']).data(iSweep)+0e9/128]);
      bias1 =  obj.sweep.data.([obj.scId '_sweep_bias1']).data(idxBias);
      bias2 =  obj.sweep.data.([obj.scId '_sweep_bias2']).data(idxBias);
      % Find current values (biasRes) corresponding to voltages
      biasRes1 = NaN(size(voltage1)); biasRes2 = biasRes1;
      for i=1:length(idxBias)
        if i == length(idxBias)
          ii = tlim(Epoch, irf.tint([eBias.stop.toUtc,'/',sweepTime.stop.toUtc]));
        else, ii = tlim(Epoch, eBias(i+[0 1]));
        end
        biasRes1(ii) = bias1(i);  biasRes2(ii) = bias2(i);
      end
    end % getSweep
    
  end % methods(Access=private)
  
end

