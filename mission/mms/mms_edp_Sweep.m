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
  % SPDX-License-Identifier: Beerware
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
        irf.log('critical', msg); error(msg);
      end
      obj.sweep = dataobj(fileName);
      if isempty(obj.sweep)
        msg = 'No data loaded'; irf.log('critical', msg); error(msg);
      end
      obj.scId = lower(obj.sweep.GlobalAttributes.Source_name{:}(1:4));
      obj.nSweeps = obj.sweep.data.([obj.scId, '_sweep_start']).nrec;
      obj.pTable = [];
      obj.p1 = []; obj.p2 = []; obj.p3 = []; obj.p4 = []; obj.p5 = []; obj.p6 = [];
      for iSweep = 1:obj.nSweeps
        [sweepTime, prb1, prb2, voltage1, biasRes1, ~, ~, eVolt, eBias, ...
          ~, ~] = obj.getSweep(iSweep);
        % Determine type of sweep:
        % ++ = up only, -- = down only, +- = up-down, -+ = down-up, 00 = wiggly
        bias_no_NaN = biasRes1(isfinite(biasRes1));
        n = length(bias_no_NaN); n2 = round(n/2);
        if n>1
          p = polyfit((1:n)',    bias_no_NaN(1:n),    1); slope = p(1);
          p = polyfit((1:n2)',   bias_no_NaN(1:n2),   1); slope1 = p(1);
          p = polyfit((n2+1:n)', bias_no_NaN(n2+1:n), 1); slope2 = p(1);
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
          if any(isnan(biasRes1))
            nanstr = [num2str(sum(isnan(biasRes1))), ' NaNs '];
          else
            nanstr = '';
          end
          if length(voltage1) > 1
            millis = [
              num2str(1e3*(eVolt(1)-sweepTime.start), '%6.3f'), ' ', ...
              num2str(1e3*(eVolt(2)-sweepTime.start), '%6.3f'), ' ', ...
              num2str(1e3*(eBias(1)-sweepTime.start), '%6.3f'), ' ', ...
              num2str(1e3*(eBias(2)-sweepTime.start), '%6.3f'), ' ', ...
              num2str(round(1/(eVolt(2)-eVolt(1)))), ' ', ...
              num2str(round(1/(eBias(2)-eBias(1)))), ' ' ];
          else
            millis = '';
          end
          disp([num2str(iSweep, '%3.2d'), ': ', ...
            sweepTime.start.toUtc, ' ', ...
            num2str((sweepTime.stop-sweepTime.start), '%3.2f'), 's ', ...
            num2str(prb1), ' ', num2str(prb2), ' ', ...
            'tbl=', num2str(obj.sweep.data.([obj.scId, '_sweep_table']).data(iSweep)), ' ', ...
            'dwl=', num2str(obj.sweep.data.([obj.scId, '_sweep_dwell']).data(iSweep)), ' ', ...
            'stp=', num2str(obj.sweep.data.([obj.scId, '_sweep_steps']).data(iSweep)), ' ', ...
            'len=', num2str(length(voltage1), '%4.4d'), ' ', ...
            'rt=', num2str(length(voltage1)/...
            (double(obj.sweep.data.([obj.scId, '_sweep_dwell']).data(iSweep))*...
            double(obj.sweep.data.([obj.scId, '_sweep_steps']).data(iSweep)))), ' ', ...
            'typ=', type, ' ', ...
            '[', num2str(round(min(biasRes1(isfinite(biasRes1))))), ',', ...
            num2str(round(max(biasRes1(isfinite(biasRes1))))), '] ', ...
            nanstr, ...
            millis, ...
            ]);
        end
        if isempty(voltage1)
          disp(['*** Warning, empty sweep ', num2str(iSweep), ' of ', num2str(obj.nSweeps)]);
          irf.log('warning', ['Empty sweep ', num2str(iSweep), ' of ', num2str(obj.nSweeps)]);
        end
        tmp.time = sweepTime;
        tmp.iPh = NaN; % Old original "iPh" value
        tmp.iPh_knee = NaN; % New knee current value, new "2021 iPh"
        tmp.impedance = NaN;
        tmp.optimal_bias = NaN;
        tmp.optimal_gradient = NaN;
        tmp.phase = NaN;
        tmp.phase_knee = NaN;
        tmp.type = type;
        tmp.impedance_src = ''; % keep track source used for ".impedance" (is it from "plot_time()" or "analyse()")
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
      % "iPh_knee" (current value of the typical sweep "knee"),
      % and spin phase of the "knee" for one sweep in a sweep object obj
      %
      % analyze(obj,iSweep[,sps])
      %
      % Input: iSweep - sweep number in the file
      %        sps is a sunpulse structure obtained from
      %        sunpulse_from_hk101() or spin phase information from
      %        ancillary defatt files.
      % If using ancillary defatt, it can be loaded first using:
      %  sps = mms_load_ancillary(defattFileName, 'defatt');
      %  idxBad = diff(sps.time)==0; % Identify duplicates to be removed
      %  sps.time(idxBad) = []; sps.zphase(idxBad) = [];
      %
      % If sps is absent no phase is computed
      narginchk(2,3);
      if(nargin==2), doPhase=false; else, doPhase=true; end
      [sweepTime, prb1, ~, voltage1, biasRes1, voltage2, biasRes2,...
        ~, ~, v01, v02] = getSweep(obj,iSweep);
      if isempty(voltage1)
        disp(['*** Warning, not analyzed empty sweep ', num2str(iSweep), ...
          ' of ', num2str(obj.nSweeps)]);
        return
      end
      if size(voltage1,1) <2
        disp(['*** Warning, not analyzed almost empty sweep ', num2str(iSweep), ...
          ' of ', num2str(obj.nSweeps)]);
        return
      end
      % Determine type of sweep:
      % ++ = up only, -- = down only, +- = up-down, -+ = down-up, 00 = wiggly
      type = obj.(['p' num2str(obj.pTable(1,iSweep))])(obj.pTable(2,iSweep)).type;
      if (strcmp(type,'00') || strcmp(type, '**'))
        disp(['*** Warning, not analyzed type ', type, ' sweep ', num2str(iSweep), ...
          ' of ', num2str(obj.nSweeps)]);
        return
      end
      if (prb1 == 1 || prb1 == 3) && ...
          round(min(biasRes1)) == -250 && round(max(biasRes1)) == 100
        %% ThoNi 2021/11: When is/was this shift ever used?
        biasRes1 = (biasRes1+250)/(100+250)*(50+460) - 460;
        biasRes2 = (biasRes2+250)/(100+250)*(50+460) - 460;
        disp('*** Warning, biasRes tables changed from [-250,100] to [-460,50]');
      end
      % Compute photoemission iPh (and its timestamp) and impedance dV/dI
      tmp1 = compute_IPh_Impedance(voltage1, v01, biasRes1, type, sweepTime);
      % do the same for the other probe
      tmp2 = compute_IPh_Impedance(voltage2, v02, biasRes2, type, sweepTime);
      % Default NaN-filled phase
      if ismember(type, {'+-', '-+'})
        ph_knee_tmp1 = {NaN, NaN};
        ph_knee_tmp2 = {NaN, NaN};
        ph_tmp1 = {NaN, NaN};
        ph_tmp2 = {NaN, NaN};
      else
        ph_knee_tmp1 = NaN;
        ph_knee_tmp2 = NaN;
        ph_tmp1 = NaN;
        ph_tmp2 = NaN;
      end
      if doPhase
        switch type
          case '00'
            disp(['*** Warning, phase not computed for type ', type, ...
              ' sweep ', num2str(iSweep), ' of ', num2str(obj.nSweeps)]);
            phaseTimesIn = [];
          case {'+-', '-+'}
            % sweep has been split
            phaseTimesIn = {tmp1.iPh_knee_time{1}, tmp1.iPh_knee_time{2}, ...
              tmp2.iPh_knee_time{1}, tmp2.iPh_knee_time{2}, ...
              tmp1.iPh_time{1}, tmp1.iPh_time{2}, ...
              tmp2.iPh_time{1}, tmp2.iPh_time{2}};
            phasesOut = {'ph_knee_tmp1{1}', 'ph_knee_tmp1{2}', ...
              'ph_knee_tmp2{1}', 'ph_knee_tmp2{2}', ...
              'ph_tmp1{1}', 'ph_tmp1{2}', ...
              'ph_tmp2{1}', 'ph_tmp2{2}'};
          case {'++', '--'}
            % Plain "++" or "--" sweeps
            phaseTimesIn = {tmp1.iPh_knee_time, tmp2.iPh_knee_time, ...
              tmp1.iPh_time, tmp2.iPh_time};
            phasesOut = {'ph_knee_tmp1', 'ph_knee_tmp2', ...
              'ph_tmp1', 'ph_tmp2'};
          otherwise
            irf.log('warning', ['Unexpected type', type,' to compute z-phase for...']);
            phaseTimesIn = [];
        end
        if isfield(sps, 'zphase')
          % phase from defatt
          for iPhase = 1:length(phaseTimesIn)
            if ~isnan(phaseTimesIn{iPhase})
              phase = mms_defatt_phase(sps, phaseTimesIn{iPhase}); %#ok<NASGU>
              eval([phasesOut{iPhase}, '=phase.data;']);
            end
          end
        else
          % phase from sunpulses
          for iPhase = 1:length(phaseTimesIn)
            if ~isnan(phaseTimesIn{iPhase})
              [phase, ~] = mms_sdp_phase_2(sps, phaseTimesIn{iPhase}); %#ok<ASGLU>
              eval([phasesOut{iPhase}, '=phase;']);
            end
          end
        end
      end % doPhase
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
      if strcmp(obj.(p_1)(obj.pTable(2,iSweep)).impedance_src, 'time')
        % we have impedance computed from "plot_time", give warning before
        % replacing it.
        irf.log('warning', 'Replacing "impedance" computed in plot_time() with computed value from analyse().');
      end
      obj.(p_1)(obj.pTable(2,iSweep)).impedance_src = 'analyse';
      obj.(p_2)(obj.pTable(2,iSweep)).impedance_src = 'analyse';
      obj.(p_1)(obj.pTable(2,iSweep)).impedance = tmp1.impedance;
      obj.(p_2)(obj.pTable(2,iSweep)).impedance = tmp2.impedance;
      obj.(p_1)(obj.pTable(2,iSweep)).optimal_gradient = tmp1.optimal_gradient;
      obj.(p_2)(obj.pTable(2,iSweep)).optimal_gradient = tmp2.optimal_gradient;
      obj.(p_1)(obj.pTable(2,iSweep)).optimal_bias = tmp1.optimal_bias;
      obj.(p_2)(obj.pTable(2,iSweep)).optimal_bias = tmp2.optimal_bias;

      if ismember(type, {'+-', '-+'})
        obj.(p_1)(obj.pTable(2,iSweep)).phase_knee = {mod(ph_knee_tmp1{1}+angle1+90,360)-90, mod(ph_knee_tmp1{2}+angle1+90,360)-90};
        obj.(p_2)(obj.pTable(2,iSweep)).phase_knee = {mod(ph_knee_tmp2{1}+angle2+90,360)-90, mod(ph_knee_tmp2{2}+angle2+90,360)-90};
        obj.(p_1)(obj.pTable(2,iSweep)).phase = {mod(ph_tmp1{1}+angle1+90,360)-90, mod(ph_tmp1{2}+angle1+90,360)-90};
        obj.(p_2)(obj.pTable(2,iSweep)).phase = {mod(ph_tmp2{1}+angle2+90,360)-90, mod(ph_tmp2{2}+angle2+90,360)-90};
      else
        obj.(p_1)(obj.pTable(2,iSweep)).phase_knee = mod(ph_knee_tmp1+angle1+90,360)-90;
        obj.(p_2)(obj.pTable(2,iSweep)).phase_knee = mod(ph_knee_tmp2+angle2+90,360)-90;
        obj.(p_1)(obj.pTable(2,iSweep)).phase = mod(ph_tmp1+angle1+90,360)-90;
        obj.(p_2)(obj.pTable(2,iSweep)).phase = mod(ph_tmp2+angle2+90,360)-90;
      end
      % New "iPh knee" values
      obj.(p_1)(obj.pTable(2,iSweep)).iPh_knee = tmp1.iPh_knee;
      obj.(p_2)(obj.pTable(2,iSweep)).iPh_knee = tmp2.iPh_knee;
      return

      % Help function
      function tmp = compute_IPh_Impedance(voltage, v0, biasRes, type, sweepTime)
        % Compute photoemission iPh, iPh_knee (and its timestamp) and impedance dV/dI
        narginchk(5,5);
        %% Default output
        tmp = struct('iPh', NaN);
        if ismember(type, {'+-', '-+'})
          % iPh knee and impedance on "up" and "down" part separately
          tmp.iPh_knee = {NaN, NaN};
          tmp.impedance = {NaN, NaN};
          tmp.optimal_gradient = {NaN, NaN};
          tmp.optimal_bias = {NaN, NaN};
          tmp.iPh_knee_time = {NaN, NaN}; % Time of derived iPh_knee (first and second part)
          tmp.iPh = {NaN, NaN};
          tmp.iPh_time = {NaN, NaN};  % Time of derived iPh (first & second part)
        else
          tmp.iPh_knee = NaN;
          tmp.impedance = NaN;
          tmp.optimal_gradient = NaN;
          tmp.optimal_bias = NaN;
          tmp.iPh_knee_time = NaN;
          tmp.iPh_time = NaN;
        end

        %% Find settling time (with fallback max 10% of sweep) based on
        % voltage response ("--" sweep should have negative gradient "++"
        % sweep should have positive gradient of voltage, so discard the
        % first part of the sweep not as aligned with this expectation).
        % Also exclude any NaN's returned from getSweep
        switch type
          case {'++', '+-'}
            % Initially increasing bias current => implies increasing probe
            % voltage response
            settleTime = find(diff(voltage)>0 & ~isnan(biasRes(1:end-1)), 1, 'first');
          case {'--', '-+'}
            % Initially decreasing bias current => implies decreasing probe
            % voltage resonse
            settleTime = find(diff(voltage)<0 & ~isnan(biasRes(1:end-1)), 1, 'first');
          otherwise
            % type="00"?. For now fallback to 8.
            settleTime = 8;
        end
        if settleTime > 0.1*length(voltage)
          % Fallback to 8 or first non-NaN, if computed settleTime is more
          % than 10% of the entire sweep. (some bad I/V sweep responses)
          settleTime = max(8, find(~isnan(biasRes), 1, 'first'));
        end
        voltageSegm = voltage(settleTime:end);
        biasResSegm = biasRes(settleTime:end);
        oldFallback = false; % <- Is set to true if new method fails.

        %% Split up and down part of sweep
        switch type
          case {'+-', '-+'}
            % Split upward and downward part
            splitInd = length(biasRes)/2-settleTime;
            % first half
            biasFirst = biasResSegm(1:splitInd);
            voltFirst = voltageSegm(1:splitInd);
            [tmp.iPh_knee{1}, tmp.optimal_gradient{1}, tmp.optimal_bias{1}] = compute_IPh_knee(voltFirst, biasFirst, type);
            % second half
            biasSecond = biasResSegm(splitInd+1:end);
            voltSecond = voltageSegm(splitInd+1:end);
            [tmp.iPh_knee{2}, tmp.optimal_gradient{2}, tmp.optimal_bias{2}] = compute_IPh_knee(voltSecond, biasSecond, type);
            if isnan(tmp.iPh_knee{1}) || isnan(tmp.iPh_knee{2})
              % Any of the new values not ok. Fallback to old method
              oldFallback = true;
            else
              % iPh_knee_time{1} time of iPh_knee for the first part of the sweep,
              % iPh_knee_time{2} time of iPh_knee for the second part of the sweep
              % (computed similarly but segment start at half time of the
              % sweep for second half).
              tmp.iPh_knee_time{1} = sweepTime.start.epoch ...
                + (-tmp.iPh_knee{1} - biasRes(1)) / (biasRes(splitInd)-biasRes(1)) ...
                * ((sweepTime.stop.epoch-sweepTime.start.epoch)/2);
              tmp.iPh_knee_time{2} = sweepTime.start.epoch ...
                + (-tmp.iPh_knee{2} - biasRes(splitInd+1)) / (biasRes(end)-biasRes(splitInd+1)) ...
                * ((sweepTime.stop.epoch-sweepTime.start.epoch)/2) ...
                + (sweepTime.stop.epoch-sweepTime.start.epoch)/2;
            end
          case {'++', '--'}
            % Monotone upward or downward, no need to split
            [tmp.iPh_knee, tmp.optimal_gradient, tmp.optimal_bias] = compute_IPh_knee(voltageSegm, biasResSegm, type);
            if isnan(tmp.iPh_knee)
              logStr = 'Failed new "knee" method, fallback to old method.';
              irf.log('debug', logStr);
              oldFallback = true;
            else
              % Compute time of iPh_knee
              tmp.iPh_knee_time = sweepTime.start.epoch ...
                + (-tmp.iPh_knee - biasRes(1)) / (biasRes(end)-biasRes(1)) ...
                * (sweepTime.stop.epoch-sweepTime.start.epoch);
            end
          otherwise
            % "type" is not supported by compute_IPh_knee yet, use old
            % method as a fallback
            oldFallback = true;
        end

        %% Old method
        if ismember(type, {'+-', '-+'})
          % Split upward and downward part (was not done in the original
          % old method)
          splitInd = length(biasRes)/2-8;
          % first half
          biasFirst = biasResSegm(1:splitInd);
          voltFirst = voltageSegm(1:splitInd);
          % Find maximum negative voltage
          vmin = min(voltFirst);
          [~, old_indIph] = min(abs(voltFirst-(2*v0+vmin)/3));
          tmp.iPh{1} = -biasFirst(old_indIph);
          % Compute time of iPh
          timeArray = int64(linspace(0, double(sweepTime.stop.ttns-sweepTime.start.ttns), length(biasRes))) + sweepTime.start.ttns;
          tmp.iPh_time{1} = timeArray(8+old_indIph);
          % second half
          biasSecond = biasResSegm(splitInd+1:end);
          voltSecond = voltageSegm(splitInd+1:end);
          [~, old_indIph] = min(abs(voltSecond-(2*v0+vmin)/3));
          tmp.iPh{2} = -biasSecond(old_indIph);
          tmp.iPh_time{2} = timeArray(splitInd + old_indIph);
        else
          % basic one way no need to split
          % Find maximum negative voltage
          vmin = min(voltage(9:end));
          [~, old_indIph] = min(abs(voltage(9:end)-(2*v0+vmin)/3));
          tmp.iPh = -biasRes(8+old_indIph);
          % Compute time of iPh
          timeArray = int64(linspace(0, double(sweepTime.stop.ttns-sweepTime.start.ttns), length(biasRes))) + sweepTime.start.ttns;
          tmp.iPh_time = timeArray(8+old_indIph);
        end
        if oldFallback
          if ismember(type, {'+-', '-+'})
            tmp.iPh_knee{1} = tmp.iPh;
            tmp.iPh_knee{2} = tmp.iPh;
            % tmp.iPh_knee_time remains default "{NaN, NaN}"
          else
            tmp.iPh_knee = tmp.iPh;
            tmp.iPh_knee_time = sweepTime.start.epoch ...
              + (-tmp.iPh - biasRes(1)) / (biasRes(end)-biasRes(1)) ...
              * (sweepTime.stop.epoch-sweepTime.start.epoch);
          end
        end

        % TODO:
        % Future improvements could be to possibly try to compute "optimal"
        % dV/dI (max value) per sweep instead value at 70-80% of iPh..
        % Or possibly try to make use of iPh_knee.
        if ismember(type, {'+-', '-+'})
          % first half
          [~,ind80] = min(abs(biasFirst-(-0.8*tmp.iPh{1})));
          [~,ind70] = min(abs(biasFirst-(-0.7*tmp.iPh{1})));
          if ind70-ind80 > 1
            p = polyfit(biasFirst(ind80:ind70), voltFirst(ind80:ind70), 1);
            tmp.impedance{1} = 1000*p(1);
          elseif ind70-ind80 < 1
            p = polyfit(biasFirst(ind70:ind80), voltFirst(ind70:ind80), 1);
            tmp.impedance{1} = 1000*p(1);
          end
          % second half
          [~,ind80] = min(abs(biasSecond-(-0.8*tmp.iPh{2})));
          [~,ind70] = min(abs(biasSecond-(-0.7*tmp.iPh{2})));
          if ind70-ind80 > 1
            p = polyfit(biasSecond(ind80:ind70), voltSecond(ind80:ind70), 1);
            tmp.impedance{2} = 1000*p(1);
          elseif ind70-ind80 < 1
            p = polyfit(biasSecond(ind70:ind80), voltSecond(ind70:ind80), 1);
            tmp.impedance{2} = 1000*p(1);
          end
        else
          % Find dVdI = slope of curve between -0.8*Iph and -0.7*Iph
          [~,ind80] = min(abs(biasRes-(-0.8*tmp.iPh)));
          [~,ind70] = min(abs(biasRes-(-0.7*tmp.iPh)));
          if ind70-ind80 > 1
            p = polyfit(biasRes(ind80:ind70),voltage(ind80:ind70),1);
            tmp.impedance = 1000*p(1);
          elseif ind70-ind80 < 1
            p = polyfit(biasRes(ind70:ind80),voltage(ind70:ind80),1);
            tmp.impedance = 1000*p(1);
          end
        end

      end

      function [iPhKnee, optimalGradient, optimalBias] = compute_IPh_knee(voltage, bias, type)
        % Help function to compute iPh_knee, a value which nominal
        % operating point wants to avoid..
        narginchk(3,3);
        % Default output
        iPhKnee = NaN; optimalGradient=NaN; optimalBias=NaN;
        % MMS1 2021-04-18T22:59:26 sweep 5 (p3 & p4) have some strange
        % behaviour after the saturation point making the voltage go back
        % towards zero (just a few volts), so instead of a strict
        % instrument saturation voltage (-107V) use a limit of <= -100 V for MMS.
        saturLimit = -100; % V, voltages below this value is considered saturated
        % mms1 20210402 sweep1 p2 visually have iPh knee first "spread" at 0.28 V
        % mms1 20210411 sweep1 p1 visually have iPh knee first "spread" at 0.29 V
        % mms3 20191201 sweep1 p1 visually have iPh knee first "spread" at 0.23 V
        % mms3 20191213 sweep4 p1 visually have iPh knee first "spread" at 0.26 V
        biasSpreadLimit = 0.3; % V, voltage spread limit (per bias current) for detecting iPh
        % Locate any saturated voltages (<= saturLimit V)
        voltageSaturated = voltage <= saturLimit;
        if ~all(voltageSaturated) && ismember(type, {'++', '--', '+-', '-+'})
          % Not all data was saturated
          voltageNonSatur = voltage(~voltageSaturated);
          biasNonSatur = bias(~voltageSaturated);
          biasUniq = unique(biasNonSatur);
          biasUniqLen = length(biasUniq);
          if biasUniqLen == length(biasNonSatur)
            logStr = 'Only one voltage datapoint per bias current. Sweep voltage spread-response method does NOT work!';
            irf.log('warning', logStr);
            return % return, with default NaN values (=> use fallback)
          end
          biasSpread = NaN(biasUniqLen, 1);
          for iBias=1:biasUniqLen
            biasSpread(iBias) = max(voltageNonSatur(biasNonSatur == biasUniq(iBias))) - min(voltageNonSatur(biasNonSatur == biasUniq(iBias)));
          end
          % "iPh" identification based on voltage spread should be less
          % than "biasSpreadLimit" V per bias current. And ensure it is not
          % the first datapoint (as it can be just after voltage saturation
          % and we later want to move down one step)
          indIphSpread = find(biasSpread < biasSpreadLimit & [false; true(biasUniqLen-1,1)], 1, 'first');
          if isempty(indIphSpread)
            logStr=  'Failed to detect iPh based on voltage spread-response, fallback to old method';
            irf.log('notice', logStr);
            return % return, with default NaN values (=> use fallback)
          else
            iPhKnee = -biasUniq(indIphSpread-1); % bias current just before linear part
            biasOverKnee=biasNonSatur(biasNonSatur>=(-iPhKnee)); % keeping only biases which are to be considered for gradient calculation
            biasForGradient=unique(biasOverKnee);
            voltageOverKnee=voltageNonSatur(1:length(biasOverKnee)); % getting the corresponding voltage values to calculate the average
            biasGroups=findgroups(biasOverKnee);
            averageVoltages=splitapply(@mean, voltageOverKnee, biasGroups);
            [optimalGradient, optIndex] = max(gradient(biasForGradient, averageVoltages));
            optimalBias = -biasForGradient(optIndex); % Keep currents direction convention
          end
        else
          % Unsupported type "00" or entirely saturated (e.g. MMS2 sdp2
          % after its probe failure)
          logStr = ['Either all data is saturated with voltage below ', ...
            num2str(saturLimit), 'V or "type" is not supported'];
          irf.log('notice', logStr);
          return % return, with default NaN values (=> use fallback)
        end % all data saturated? & supported type?
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
          analyze(obj, iSweep)
        end
      elseif nargin == 2
        for iSweep = 1:obj.nSweeps
          analyze(obj, iSweep, sps)
        end
      elseif nargin == 3
        for iSweep = 1:obj.nSweeps
          if printStatus, disp(['Analyzing sweep ', num2str(iSweep)]), end
          analyze(obj, iSweep, sps)
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
        eVolt, eBias, v01, v02] = getSweep(obj, iSweep);
      return
    end

    function list(obj, iSweep)
      % List sweep
      %
      % list(obj, iSweep)
      %
      % Input: iSweep - sweep number in the file
      %
      [sweepTime, ~, ~, voltage1, biasRes1, voltage2, biasRes2,...
        eVolt, eBias, v01, v02]...
        = getSweep(obj,iSweep);
      disp(['Sweep start: ', sweepTime.start.utc]);
      disp(['Sweep end:   ', sweepTime.stop.utc]);
      time = sweepTime.start + (0:length(biasRes1)-1)' * ...
        (sweepTime.stop - sweepTime.start) / ...
        length(biasRes1);
      time_1 = EpochTT(sweepTime.start.ttns - ...
        (sweepTime.stop.ttns - sweepTime.start.ttns) / length(biasRes1));
      disp([ time_1.utc, ...
        ' ', '---', ' ', '---', ...
        ' ', num2str(v01), ' ', num2str(v02) ]);
      for i = 1:length(biasRes1)
        disp([time(i).utc, ...
          ' ', num2str(biasRes1(i)), ' ', num2str(biasRes2(i)), ...
          ' ', num2str(voltage1(i)), ' ', num2str(voltage2(i)), ...
          ' ', num2str((eBias(fix((i-1)*length(eBias)/length(biasRes1))+1)-sweepTime.start)*1e3), ...
          ' ', num2str((eVolt(i)-sweepTime.start)*1e3) ]);
      end
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
        msg = 'H must be an axes handle'; irf.log('critical',msg); error(msg);
      end
      if 1>iSweep || iSweep>obj.nSweeps || round(iSweep)~=iSweep
        msg = sprintf('ISWEEP must be 1..%d',obj.nSweeps);
        irf.log('critical',msg); error(msg); %#ok<SPERR>
      end
      [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2, ...
        ~, ~, ~, ~] = getSweep(obj,iSweep);
      if isempty(voltage1)
        disp(['*** Warning, not plotted empty sweep ', num2str(iSweep), ...
          ' of ', num2str(obj.nSweeps)]);
        return
      end
      dwl = obj.sweep.data.([obj.scId, '_sweep_dwell']).data(iSweep);
      stp = obj.sweep.data.([obj.scId, '_sweep_steps']).data(iSweep);
      len = length(voltage1);
      % Determine type of sweep:
      % ++ = up only, -- = down only, +- = up-down, -+ = down-up, 00 = wiggly
      if iSweep <= length(obj.pTable)
        type = obj.(['p' num2str(obj.pTable(1,iSweep))])(obj.pTable(2,iSweep)).type;
      else
        type = '  ';
      end
      if (prb1 == 1 || prb1 == 3) && ...
          round(min(biasRes1)) == -250 && round(max(biasRes1)) == 100
        %% ThoNi 2021/11: When is/was this shift ever used?
        biasRes1 = (biasRes1+250)/(100+250)*(50+460) - 460;
        biasRes2 = (biasRes2+250)/(100+250)*(50+460) - 460;
        disp('*** Warning, biasRes tables changed from [-250,100] to [-460,50]')
      end
      % Plot
      if isempty(h), clf; else, clf(h); end
      c = 'krgbmc';
      lineStyleP1 = [c(prb1), '.-'];
      lineStyleP2 = [c(prb2), '.-'];
      if strcmp(type,'++') || strcmp(type,'--') || strcmp(type,'00')
        if isempty(h)
          plot(voltage1, biasRes1, lineStyleP1, ...
            voltage2, biasRes2, lineStyleP2, ...
            voltage1(1), NaN);
          h = gca;
        else
          plot(h, voltage1, biasRes1, lineStyleP1, ...
            voltage2, biasRes2, lineStyleP2, ...
            voltage1(1), NaN);
        end
      else
        n = length(biasRes1(isfinite(biasRes1)));
        n2 = round(n/2);
        lineStyleP3 = [c(prb1), '.--'];
        lineStyleP4 = [c(prb2), '.--'];
        if isempty(h)
          plot(voltage1(n2:n), biasRes1(n2:n), lineStyleP3, ...
            voltage2(n2:n), biasRes2(n2:n), lineStyleP4);
          hold(gca, 'on')
          plot(voltage1(1:n2), biasRes1(1:n2), lineStyleP1, ...
            voltage2(1:n2), biasRes2(1:n2), lineStyleP2, ...
            voltage1(1), NaN);
          hold(gca, 'off')
          h = gca;
        else
          plot(h, voltage1(n2:n), biasRes1(n2:n), lineStyleP3, ...
            voltage2(n2:n), biasRes2(n2:n), lineStyleP4);
          hold(gca, 'on')
          plot(h, voltage1(1:n2), biasRes1(1:n2), lineStyleP1, ...
            voltage2(1:n2), biasRes2(1:n2), lineStyleP2, ...
            voltage1(1), NaN);
          hold(gca, 'off')
        end
      end
      ylim = [min([biasRes1; biasRes2])-10 max([biasRes1; biasRes2])+10];
      set(h, 'YLim', ylim);
      t = title(h, [obj.scId, ' ', toUtc(sweepTime(1),1), ' sweep ', num2str(iSweep), ...
        ', len=', num2str(len), ' stp=', num2str(stp), ' dwl=', num2str(dwl), ...
        ' rate=', num2str(128*len/double(stp)/double(dwl))]);
      %t = title(h,[obj.scId ' ' toUtc(sweepTime(1),1) ' sweep ' num2str(iSweep)]);
      if isa(h, 'handle'), set(t, 'FontSize', 12); end % Needed for HG2
      ylabel(h, ...
        ['Bias [', getunits(obj.sweep,[obj.scId, '_sweep_bias1']), ']'])
      xlabel(h, ...
        ['Voltage [', getunits(obj.sweep,[obj.scId, '_edp_sweeps']), ']'])
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
        dVdI1_optimal = obj.(p_1)(obj.pTable(2,iSweep)).optimal_gradient;
        dVdI2_optimal = obj.(p_2)(obj.pTable(2,iSweep)).optimal_gradient;
        dVdI1bias_optimal = obj.(p_1)(obj.pTable(2,iSweep)).optimal_bias;
        dVdI2bias_optimal = obj.(p_2)(obj.pTable(2,iSweep)).optimal_bias;
        pha1 = obj.(p_1)(obj.pTable(2,iSweep)).phase;
        pha2 = obj.(p_2)(obj.pTable(2,iSweep)).phase;
        Iph1_knee = obj.(p_1)(obj.pTable(2,iSweep)).iPh_knee;
        Iph2_knee = obj.(p_2)(obj.pTable(2,iSweep)).iPh_knee;
        pha1_knee = obj.(p_1)(obj.pTable(2,iSweep)).phase_knee;
        pha2_knee = obj.(p_2)(obj.pTable(2,iSweep)).phase_knee;
      else
        Iph1 = NaN; Iph2 = NaN;
        dVdI1 = NaN; dVdI2 = NaN;
        dVdI1_optimal = NaN; dVdI2_optimal = NaN;
        dVdI1bias_optimal = NaN; dVdI2bias_optimal = NaN;
        pha1 = NaN; pha2 = NaN;
        Iph1_knee = NaN; Iph2_knee = NaN;
        pha1_knee = NaN; pha2_knee = NaN;
      end
      % Output for debugging
      % [double(prb1) double(prb2) Iph1 Iph2 dVdI1 dVdI2]
      if ismember(type, {'+-', '-+'})
        legend(h, ...
          sprintf('V%d  Iph %0.1f & %0.1f (knee: %0.1f & %0.1f) nA,  Imp: %0.2f (& %0.2f) MOhm @ %0.1f (& %0.1f) nA, Optimal Imp: %0.2f (& %0.2f) MOhm @ %0.1f (& %0.1f) nA', prb1, Iph1{1}, Iph1{2}, Iph1_knee{1}, Iph1_knee{2}, dVdI1{1}, dVdI1{2}, 0.75*Iph1{1}, 0.75*Iph1{2}, dVdI1_optimal{1}, dVdI1bias_optimal{2}, dVdI1bias_optimal{1}, dVdI1bias_optimal{2}),...
          sprintf('V%d  Iph %0.1f & %0.1f (knee: %0.1f & %0.1f) nA,  Imp: %0.2f (& %0.2f) MOhm @ %0.1f (& %0.1f) nA, Optimal Imp: %0.2f (& %0.2f) MOhm @ %0.1f (& %0.1f) nA', prb2, Iph2{1}, Iph2{2}, Iph2_knee{1}, Iph2_knee{2}, dVdI2{1}, dVdI2{2}, 0.75*Iph2{1}, 0.75*Iph2{2}, dVdI2_optimal{1}, dVdI2bias_optimal{2}, dVdI2bias_optimal{1}, dVdI2bias_optimal{2}),...
          sprintf('Phase%d  %0.1f & %0.1f (knee: %0.1f & %0.1f),  Phase%d  %0.1f & %0.1f (knee: %0.1f & %0.1f) degrees Type %s', prb1, pha1{1}, pha1{2}, pha1_knee{1}, pha1_knee{2}, prb2, pha2{1}, pha2{2}, pha2_knee{1}, pha2_knee{2}, type),...
          'Location', 'NorthWest', ...
          'AutoUpdate', 'off');
      else
        legend(h, ...
          sprintf('V%d  Iph %0.1f (knee: %0.1f) nA,  %0.2f MOhm @ %0.1f nA, Optimal Imp: %0.2f MOhm @ %0.1f nA', prb1, Iph1, Iph1_knee, dVdI1, 0.75*Iph1, dVdI1_optimal, dVdI1bias_optimal), ...
          sprintf('V%d  Iph %0.1f (knee: %0.1f) nA,  %0.2f MOhm @ %0.1f nA, Optimal Imp: %0.2f MOhm @ %0.1f nA', prb2, Iph2, Iph2_knee, dVdI2, 0.75*Iph2, dVdI2_optimal, dVdI2bias_optimal),...
          sprintf('Phase%d  %0.1f (knee: %0.1f),  Phase%d  %0.1f (knee: %0.1f) degrees Type %s', prb1, pha1, pha1_knee, prb2, pha2, pha2_knee, type),...
          'Location', 'NorthWest', ...
          'AutoUpdate', 'off');
      end
      grid(h, 'on')
      if nargout, hout = h; end
    end % PLOT

    function hout = plot_time(obj, h, iSweep, sps, evfile)
      % Plot sweep
      %
      % Plots a sweep as I(t), V(t) and dVdI(t) [ and Phase(t) ]
      %
      % hout = plot_time(obj, [h], iSweep [, sps, evfile])
      %
      % Input:
      %    iSweep - sweep number in the file
      %    h - axes handle [optional]
      %    sps - a sunpulse structure obtained from sunpulse_from_hk101()
      %    if sps is not given, no phase is plotted
      %    evfile - corresponding dce file
      %
      if nargin==4, evfile = sps; plotev = true; else, plotev = false; end
      if nargin>=3, sps = iSweep; iSweep = h; h = []; doPhase = true;
      elseif nargin == 2, iSweep = h; h = []; doPhase = false;
      elseif ~isgraphics(h, 'axes')
        msg = 'H must be an axes handle';
        irf.log('critical', msg); error(msg);
      end
      if 1>iSweep || iSweep>obj.nSweeps || round(iSweep)~=iSweep
        msg = sprintf('ISWEEP must be 1..%d', obj.nSweeps);
        irf.log('critical', msg); error(msg); %#ok<SPERR>
      end
      [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        ~, ~, v01, v02] = getSweep(obj,iSweep);
      if isempty(voltage1)
        disp(['*** Warning, not plotted empty sweep ', num2str(iSweep), ...
          ' of ', num2str(obj.nSweeps)]);
        irf.log('warning', ['Not plotting empty sweep ', num2str(iSweep), ...
          ' of ', num2str(obj.nSweeps)]);
        return
      end
      type = obj.(['p' num2str(obj.pTable(1,iSweep))])(obj.pTable(2,iSweep)).type;
      dwl = obj.sweep.data.([obj.scId, '_sweep_dwell']).data(iSweep);
      stp = obj.sweep.data.([obj.scId, '_sweep_steps']).data(iSweep);
      len = length(voltage1);
      ndwl = len/stp;
      impedance1 = NaN(size(biasRes1)); impedance2 = impedance1;
      ind_impedance = ndwl/2:ndwl:length(biasRes1)-ndwl/2;
      for i=ind_impedance
        i1=i+0*ndwl-1-(mod(i-1, ndwl));
        i2=i+1*ndwl-1-(mod(i-1, ndwl));
        i3=i+2*ndwl-1-(mod(i-1, ndwl));
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
      time = sweepTime.start + (0:length(biasRes1)-1)' * ...
        (sweepTime.stop - sweepTime.start) / length(biasRes1);
      clf
      switch prb1
        case 1
          p_1='p1'; p_2='p2'; angle1=30; angle2=210; %#ok<NASGU>
        case 2
          p_1='p2'; p_2='p1'; angle1=210; angle2=30; %#ok<NASGU>
        case 3
          p_1='p3'; p_2='p4'; angle1=120; angle2=300; %#ok<NASGU>
        case 4
          p_1='p4'; p_2='p3'; angle1=300; angle2=120; %#ok<NASGU>
        case 5
          p_1='p5'; p_2='p6'; angle1=0; angle2=0; %#ok<NASGU>
        case 6
          p_1='p6'; p_2='p5'; angle1=0; angle2=0; %#ok<NASGU>
        otherwise
          error('Bad probe number');
      end
      if doPhase
        if isfield(sps, 'zphase')
          phase = mms_defatt_phase(sps, time.ttns);
          ph_tmp = phase.data;
        else
          [ph_tmp, ~] = mms_sdp_phase_2(sps, time.ttns);
        end
        phase1 = mod(ph_tmp+angle1+90, 360)-90;
        % phase2 = mod(ph_tmp+angle2+90, 360)-90;
        if plotev
          evobj = dataobj(evfile);
          efield = get_ts(evobj, [obj.scId, '_edp_dce_sensor']);
          vfield = get_ts(evobj, [obj.scId, '_edp_dcv_sensor']);
          ind = find(vfield.time.ttns >= 2*time.start.ttns - time(2).ttns & ...
            vfield.time.ttns <= 2*time.stop.ttns - time(end-1).ttns);
          switch prb1
            case 1
              v12 = irf.ts_scalar(efield.time(ind), efield.data(ind, 1));
              v1 = irf.ts_vec_xy(efield.time(ind), [vfield.data(ind, 1), vfield.data(ind, 1)-efield.data(ind, 1)*0.12]); % SDP distance 120 m (*10^-3 since mV->V)
            case 3
              v12 = irf.ts_scalar(efield.time(ind), efield.data(ind, 2));
              v1 = irf.ts_vec_xy(efield.time(ind), [vfield.data(ind, 2), vfield.data(ind, 2)-efield.data(ind, 2)*0.12]);
            case 5
              v12 = irf.ts_scalar(efield.time(ind), efield.data(ind, 3));
              v1 = irf.ts_vec_xy(efield.time(ind), [vfield.data(ind, 3), vfield.data(ind, 3)-efield.data(ind, 3)*0.0292]); % ADP distance 29.2 m
            otherwise, error('Bad probe number')
          end
          h1=irf_plot({...
            irf.ts_scalar(time, biasRes1), ...
            irf.ts_vec_xy(EpochTT([2*time(1).ttns-time(2).ttns; time.ttns]), [double([v01; voltage1]), double([v02; voltage2])]), ...
            v1, ...
            v12, ...
            irf.ts_vec_xy(time(ind_impedance), [impedance1(ind_impedance), impedance2(ind_impedance)]), ...
            irf.ts_scalar(time, phase1)}, ...
            '.-');
        else
          h1 = irf_plot({...
            irf.ts_scalar(time, biasRes1), ...
            irf.ts_vec_xy(EpochTT([2*time(1).ttns-time(2).ttns; time.ttns]), [double([v01; voltage1]), double([v02; voltage2])]), ...
            irf.ts_vec_xy(time(ind_impedance), [impedance1(ind_impedance), impedance2(ind_impedance)]), ...
            irf.ts_scalar(time, phase1)}, ...
            '.-');
        end % if plotev
      else
        h1 = irf_plot({...
          irf.ts_scalar(time, biasRes1), ...
          irf.ts_vec_xy(EpochTT([2*time(1).ttns-time(2).ttns; time.ttns]), [double([v01; voltage1]), double([v02; voltage2])]), ...
          irf.ts_vec_xy(time(ind_impedance), [impedance1(ind_impedance), impedance2(ind_impedance)])}, ...
          '.-');
      end % if doPhase
      title(h1(1), sprintf('%s %s sweep %i, len=%d stp=%d dwl=%d rate=%d', ...
        obj.scId, sweepTime.start.utc(1), iSweep, len, stp, dwl, ...
        128*len/double(stp)/double(dwl)));
      ylabel(h1(1), {'Bias', ['[', getunits(obj.sweep, [obj.scId, '_sweep_bias1']), ']']});
      ylabel(h1(2), {'Voltage', ['[', getunits(obj.sweep, [obj.scId, '_edp_sweeps']), ']']});
      legend(h1(1), p_1);
      legend(h1(2), p_1, p_2);
      switch type
        case {'+-', '-+'}
          % separate upward and downward parts
          ind_imp_part1 = ind_impedance(1:length(ind_impedance)/2);
          ind_imp_part2 = ind_impedance(length(ind_impedance)/2+1:end);
          imp1 = impedance1(ind_imp_part1);
          medi1 = {median(imp1(isfinite(imp1)))};
          imp1 = impedance1(ind_imp_part2);
          medi1{2} = median(imp1(isfinite(imp1)));
          imp2 = impedance2(ind_imp_part1);
          medi2 = {median(imp2(isfinite(imp2)))};
          imp2 = impedance2(ind_imp_part2);
          medi2{2} = median(imp2(isfinite(imp2)));
        case {'--', '++', '00'}
          % strictly increasing or decreasing or (wiggle "00")
          imp1 = impedance1(ind_impedance);
          medi1 = median(imp1(isfinite(imp1)));
          imp2 = impedance2(ind_impedance);
          medi2 = median(imp2(isfinite(imp2)));
        otherwise
          % Unknown type
          irf.log('warning', ['Type:', type, ' is not yet supported. Fallback to NaN.']);
          medi1 = NaN;
          medi2 = NaN;
      end
      if strcmp(obj.(p_1)(obj.pTable(2,iSweep)).impedance_src, 'analyse')
        % We have impedance computed from "analyse", give warning before
        % replacing it.
        irf.log('warning', 'Replacing "impedance" computed in analyse() with computed value from plot_time().');
      end
      obj.(p_1)(obj.pTable(2,iSweep)).impedance_src = 'time';
      obj.(p_2)(obj.pTable(2,iSweep)).impedance_src = 'time';
      obj.(p_1)(obj.pTable(2, iSweep)).impedance = medi1;
      obj.(p_2)(obj.pTable(2, iSweep)).impedance = medi2;
      if doPhase
        if plotev
          ylabel(h1(3), {'V', ['[', getunits(obj.sweep, [obj.scId, '_edp_sweeps']), ']']});
          ylabel(h1(4), {'E', '[mV/m]'});
          ylabel(h1(5), {'dVdI', '[MOhm]'})
          ylabel(h1(6), {'Phase', '[deg]'})
          legend(h1(3), p_1, p_2)
          legend(h1(4), ['p', num2str(prb1), num2str(prb2)])
          if ismember(type, {'+-', '-+'})
            legend(h1(5), [p_1, ' median=', num2str(medi1{1}, '%4.1f'), ...
              ' & ', num2str(medi1{2}, '%4.1f')], ...
              [p_2, ' median=', num2str(medi2{1}, '%4.1f'), ' & ', ...
              num2str(medi2{2}, '%4.1f')])
          else
            legend(h1(5), [p_1, ' median=', num2str(medi1, '%4.1f')], ...
              [p_2, ' median=', num2str(medi2, '%4.1f')])
          end
          legend(h1(6), p_1)
        else
          ylabel(h1(3), {'dVdI', '[MOhm]'})
          ylabel(h1(4), {'Phase', '[deg]'})
          if ismember(type, {'-+', '+-'})
            legend(h1(3), [p_1, ' median=', num2str(medi1{1}, '%4.1f'), ...
              ' & ', num2str(medi1{2}, '%4.1f')], ...
              [p_2, ' median=', num2str(medi2{1}, '%4.1f'), ' & ', ...
              num2str(medi2{2}, '%4.1f')]);
          else
            legend(h1(3), [p_1, ' median=', num2str(medi1, '%4.1f')], ...
              [p_2, ' median=', num2str(medi2, '%4.1f')]);
          end
          legend(h1(4), p_1)
        end
      else
        ylabel(h1(3), {'dVdI', '[MOhm]'})
        if ismember(type, {'+-', '-+'})
          legend(h1(3), [p_1, ' median=', num2str(medi1{1}, '%4.1f'), ...
            ' & ', num2str(medi1{2}, '%4.1f')], ...
            [p_2, ' median=', num2str(medi2{1}, '%4.1f'), ' & ', ...
            num2str(medi2{2}, '%4.1f')]);
        else
          legend(h1(3), [p_1, ' median=', num2str(medi1, '%4.1f')], ...
            [p_2, ' median=', num2str(medi2, '%4.1f')]);
        end
      end
      if nargout, hout = h1; end
    end % PLOT_TIME

  end % methods

  methods(Access=private)

    function [sweepTime, prb1, prb2, voltage1, biasRes1, voltage2, biasRes2,...
        Epoch, eBias, v01, v02]...
        = getSweep(obj, iSweep)
      % get probe, voltage and current values for one sweep
      sweepTime = EpochTT([...
        obj.sweep.data.([obj.scId, '_sweep_start']).data(iSweep), ...
        obj.sweep.data.([obj.scId, '_sweep_stop']).data(iSweep)]);
      % Internally use "sweepTimeSeg" for extracting the data, this has a
      % margin of 42200 ns applied to the "sweep_start" time.
      % This margin is dervied from mms1 sweep file v0.6.0 for 20150427
      % which had a time jitter of <= 42198 ns between the official
      % "sweep_start" time and the "Epoch" and "epoch_sweepsamp".
      % Other files with time jitter include mms4 file v0.7.2 for 20160512
      % (with <= 100 ns) and mms3 file v0.7.2 for 20170102 (with <= 286 ns).
      % No other files have had jitter was observed in any sweep file as
      % per 2021-09-01.
      sweepTimeSeg = EpochTT(sweepTime.ttns - int64([42200; 0]));
      % Extract corresponding voltage response
      [idx, Epoch] = tlim(EpochTT(obj.sweep.data.Epoch.data), sweepTimeSeg);
      % Sanity check
      if isempty(idx)
        irf.log('warning', ['Problem extracting voltage response for sweep number ', num2str(iSweep)]);
      end
      prb1 = obj.sweep.data.([obj.scId, '_sweep_swept']).data(iSweep);
      % The "other probe" is the other probe in the pair 1-2, 3-4, 5-6
      if fix(prb1/2)*2==prb1, prb2 = prb1 - 1; else, prb2 = prb1 + 1; end
      voltage1 = obj.sweep.data.([obj.scId, '_edp_sweeps']).data(idx, prb1);
      voltage2 = obj.sweep.data.([obj.scId, '_edp_sweeps']).data(idx, prb2);
      if ~isempty(voltage1)
        v01 = obj.sweep.data.([obj.scId, '_edp_sweeps']).data(idx(1)-1, prb1);
        v02 = obj.sweep.data.([obj.scId, '_edp_sweeps']).data(idx(1)-1, prb2);
      else
        v01 = NaN; v02 = NaN;
      end
      % Extract corresponding bias currents
      % (could potentially use cummulative "sweep_steps" instead of time)
      [idxBias, eBias] = ...
        tlim(EpochTT(obj.sweep.data.epoch_sweepsamp.data), sweepTimeSeg);
      % Sanity check
      if isempty(idxBias)
        irf.log('critical', ['Problem extracting bias currents for sweep number ', num2str(iSweep)]);
      else
        if idxBias(1) ~= sum(obj.sweep.data.([obj.scId, '_sweep_steps']).data(1:iSweep-1))+1
          % Unexpected start point compared with cummulative "sweep_steps"
          % Print warning (to notify user that perhaps the margins must be
          % updated).
          logStr = ['Sweep ', num2str(iSweep), ' start appears to be off by ', ...
            num2str(obj.sweep.data.([obj.scId, '_sweep_start']).data(iSweep) - obj.sweep.data.epoch_sweepsamp.data(idxBias(1)-1)), ...
            ' ns.'];
          irf.log('warning', logStr);
        elseif length(idxBias) ~= obj.sweep.data.([obj.scId, '_sweep_steps']).data(iSweep)
          % Unexpected length
          logStr = ['Sweep ', num2str(iSweep), ...
            ' extracted bias currents (len: ', num2str(length(idxBias)), ...
            ') does not align with expected number of steps ', ...
            num2str(obj.sweep.data.([obj.scId, '_sweep_steps']).data(iSweep))];
          irf.log('warning', logStr);
        end
      end
      bias1 = obj.sweep.data.([obj.scId, '_sweep_bias1']).data(idxBias);
      bias2 = obj.sweep.data.([obj.scId, '_sweep_bias2']).data(idxBias);
      % Find current values (biasRes) corresponding to voltages
      biasRes1 = NaN(size(voltage1)); biasRes2 = biasRes1;
      for i=1:length(idxBias)
        if i == length(idxBias)
          ii = tlim(Epoch, irf.tint([eBias.stop.toUtc, '/', sweepTime.stop.toUtc]));
        else, ii = tlim(Epoch, eBias(i+[0 1]));
        end
        biasRes1(ii) = bias1(i);  biasRes2(ii) = bias2(i);
      end
    end % getSweep

  end % methods(Access=private)

end
