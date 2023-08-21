function mms_sweep_plot_cron(dayToRun)
% dayToRun - date to process (format "YYYYMMDD")
year = dayToRun(1:4);
tint = irf.tint([year,'-',dayToRun(5:6),'-',dayToRun(7:8), 'T00:00:00.000000000Z'], ...
  [year,'-',dayToRun(5:6),'-',dayToRun(7:8), 'T23:59:59.999999999Z']);
dataRoot='/data/mms/';
sweepFolder = fullfile(dataRoot, 'irfu', 'plots', 'edp', 'sweeps', filesep);
mms.db_init('local_file_db', dataRoot);
lineStyle = {'black.', 'red.', 'green.', 'blue.', 'magenta.', 'cyan.'};
nowStr = char(datetime("now", "TimeZone","UTC", "Format","uuuu-MM-dd"));
%   % F10.7 From Swarm database (this contains all published data from 1947
%   % and forward and is updated by a separate daily cronjob).
%   % Source:
%   % ftp://ftp.geolab.nrcan.gc.ca/data/solar_flux/daily_flux_values/fluxtable.txt
%   %
%   % Exclude 'NULL' values and only reqest julian dates after 2015/01/01
%   (pre-launch of MMS), and let the database ensure it is monotonically
%   increasing with time.
%   sqlquery = 'SELECT julianDate,f FROM f107org WHERE NOT f ISNULL AND julianDate>=2457023.5 ORDER BY julianDate ASC;';
%   conn = sqlite('/data/swarm/sqlitedb/indices/f107_sqlite.db');
%   res = fetch(conn, sqlquery);
%   close(conn);
%   % Convert from Juliandate to EpochTT (note, julian date does not
%   % properly handle leap seconds but the F10.7 measurements are only done
%   % three times per day and our sweeps are typically done two times per
%   % day so it is perhaps not critical that we get seconds and fractions
%   % of second correct here.
%   a = datetime([res{:,1}]', 'ConvertFrom', 'juliandate');
%   f107_time = EpochTT(irf_time([a.Year, a.Month, a.Day, a.Hour, a.Minute, a.Second],'vector6>ttns'));
%   f107_ts = irf.ts_scalar(f107_time, [res{:,2}]');
%   figure('units','normalized','outerposition',[0 0 1 1]);
%   irf_plot(f107_ts);
%   ylabel({'F10.7','[s.f.u.]'});
%   title(['Plot created: ',nowStr,'. F10.7 vs time.']);
%   print(gcf, '-dpng', [sweepFolder,'summary_plots/F107vsTime.png']);
%   close all % close plots

for iSc = 1:4
  scStr = num2str(iSc);
  listSweep = mms.db_list_files(['mms', scStr, '_edp_srvy_l1b_sweeps'], tint);
  if size(listSweep)>1, error('Unexpected sweep file'); end
  if isempty(listSweep), irf.log('warning', 'No sweep file found'); continue; end
  %sObj = mms_edp_Sweep([listSweep(1).path, filesep, listSweep(1).name], 1);
  sInfo = evalc('sObj=mms_edp_Sweep([listSweep(1).path,filesep,listSweep(1).name],1);');
  fid = fopen(fullfile(sweepFolder, 'logs', ['all_sweep', scStr, '.txt']),'a');
  fprintf(fid, [listSweep(1).name, ':\n']);
  fprintf(fid, '%s \n', string(sInfo));
  fclose(fid);
  % Locate defatt file(-s) for date and s/c
  listDefatt = mms.db_list_files(['mms', scStr, '_ancillary_defatt'], tint);
  for ii=1:length(listDefatt)
    irf.log('debug', ['Z-Phase using: ', listDefatt(ii).name]);
    dataTmp = mms_load_ancillary([listDefatt(ii).path, filesep, ...
      listDefatt(ii).name], 'defatt');
    idxBad = diff(dataTmp.time)==0; % Identify first duplicate
    dataTmp.time(idxBad) = []; dataTmp.zphase(idxBad) = [];
    if ii==1
      defatt.time = dataTmp.time; defatt.zphase = dataTmp.zphase;
    else
      % multiple defatt, combine each field and run sort & unique on time
      combined = [defatt.time; dataTmp.time];
      [~, srt] = sort(combined);
      [defatt.time, usrt] = unique(combined(srt));
      combined = [defatt.zphase; dataTmp.zphase];
      defatt.zphase = combined(srt(usrt));
    end % if ii==1
  end % for ii=1:length(listDefatt)
  if isempty(listDefatt)
    irf.log('warning', 'Did not find any corresponding defatt, fall back to HK101.');
    tint = tint + [int64(-5*60*10^9) int64(5*60*10^9)]; % Add 5 minutes margin (to have some hk101 before start of day and some after end of day)
    defatt = mms.sunpulse_from_hk101(iSc, tint); % Save as "defatt" (meaning it can be used directly below instead of what was defatt..)
    if isempty(defatt)
      error('Unknown zphase');
    end
  end
  % Maneuvers or Eclipses
  [maneuvers, timeline, eclipse] = mms_maneuvers(tint, iSc);
  eclipseOrManeuver = [eclipse; maneuvers];
  if(iSc==4 && tint.stop.ttns > EpochTT('2016-06-12T05:28:48.200Z').ttns)
    warning('off','MATLAB:polyfit:PolyNotUnique');  % MMS4 p4 failed
  elseif(iSc==2 && tint.stop.ttns > EpochTT('2018-09-21T06:04:45.810Z').ttns)
    warning('off','MATLAB:polyfit:PolyNotUnique');  % MMS2 p2 failed
  elseif(iSc==3 && tint.stop.ttns > EpochTT('2022-05-08T23:18:44.300Z').ttns)
    warning('off','MATLAB:polyfit:PolyNotUnique');  % MMS3 p2 failed
  else
    warning('on','MATLAB:polyfit:PolyNotUnique');  % Other s/c (with working probes) should warn
  end
  analyze_all(sObj, defatt, 0, eclipseOrManeuver); % Process the sweep data, with zphase from defatt if available, no extra debug messages, and eclipseOrManeuver information (if available).
  figure('units', 'normalized', 'outerposition', [0 0 1 1]); % maximized window
  for iSweep=1:sObj.nSweeps
    % Plot individual sweeps
    plot(sObj, iSweep);
    % Create output folders (in case of new year)
    outFolder = sprintf('%sindividual/MMS%i_probe%i-%i/%s', sweepFolder, ...
      iSc, sObj.pTable(1,iSweep), sObj.pTable(1,iSweep)+1, year);
    if ~exist(outFolder,'dir'), mkdir(outFolder); end
    fName = sprintf('%s/mms%i_%i_%s_%03i.png', outFolder, iSc, ...
      sObj.pTable(1,iSweep), dayToRun, iSweep);
    print(gcf, '-dpng', fName);
    clf('reset')
  end % for iSweep=1:sObj.nSweeps

  %% CREATE TSeries
  c_eval('p?_time=[sObj.p?.time];p?_time=p?_time(1:2:end);', 1:6); % Create time array keeping only the start time
  c_eval('S.p?_impedance_ts=irf.ts_scalar(p?_time,[sObj.p?.impedance]);', 1:6); % TSeries of impedance
  c_eval('S.p?_iPh_ts=irf.ts_scalar(p?_time,[sObj.p?.iPh]);', 1:6); % TSeries of iPh
  c_eval('S.p?_iPh_knee_ts=irf.ts_scalar(p?_time,[sObj.p?.iPh_knee]);', 1:6); % TSeries of iPh_knee
  c_eval('S.p?_optimal_bias_ts=irf.ts_scalar(p?_time,[sObj.p?.optimal_bias]);', 1:6); % TSeries of optimal_bias
  c_eval('S.p?_optimal_gradient_ts=irf.ts_scalar(p?_time,[sObj.p?.optimal_gradient]);', 1:6); % TSeries of optimal_gradient
  c_eval('S.p?_phase_ts=irf.ts_scalar(p?_time,mod([sObj.p?.phase],360));', 1:6); % TSeries of phase
  c_eval('S.p?_phase_knee_ts=irf.ts_scalar(p?_time,mod([sObj.p?.phase_knee],360));', 1:6); % TSeries of phase
  c_eval('S.p?_type_ts=irf.ts_vec_xy(p?_time, double(reshape([sObj.p?.type],2,length(sObj.p?))''));',1:6); % TSeries of type ('--' or '++' etc converted from ASCII char to double)
  c_eval('S.p?_eclipseManeuver_ts=irf.ts_scalar(p?_time,[sObj.p?.eclipseManeuver]);', 1:6); % TSeries of eclipse and/or Maneuver time overlap with sweep
  %save(sprintf('%sobj/mms%i_%s_sweepTsObj.mat', sweepFolder, iSc, dayToRun), 'S');

  %% Load the combined data ("Sw") and combine with the data in "S".
  savedSweepTS = fullfile(sweepFolder, 'obj', ['mms', scStr, '_SweepTsCombined.mat']);
  if exist(savedSweepTS ,'file')
    load(savedSweepTS, 'Sw');
    fName = fieldnames(S);
    for iField=1:length(fName)
      Sw.(fName{iField}) = combine(Sw.(fName{iField}), S.(fName{iField}));
    end
  else
    Sw = S;
  end
  save(savedSweepTS, 'Sw');
  close all % close plots..

  %% Load DAC settings
  savedDac = fullfile(dataRoot, 'irfu', 'plots', 'edp', 'DAC', 'obj', ...
    ['mms', scStr, '_dacTsCombined.mat']);
  if exist(savedDac, 'file')
    load(savedDac, 'p1_dac', 'p2_dac', 'p3_dac', 'p4_dac', 'p5_dac', 'p6_dac'); %#ok<NASGU>
    % iPh & iPh_knee have opposite current direction as compared with DAC
    % so simply flip DAC values.
    c_eval('p?_dac.data = -p?_dac.data;', 1:6);
    plotDAC=true;
  else
    plotDAC=false;
  end

  %% Plot all combined data
  % Moving median, Matlab built in uneven sampling time taken into account...
  c_eval('t?=datetime(Sw.p?_iPh_ts.time.utc(''yyyy-mm-ddTHH:MM:SS.mmmZ''),''InputFormat'',''uuuu-MM-dd''''T''''HH:mm:ss.SSS''''Z'''''',''TimeZone'',''UTCLeapSeconds'');', 1:6);
  c_eval('indEclipse?=Sw.p?_eclipseManeuver_ts.data==1;', 1:6);

  % iPh vs time
  %c_eval('iPh_movm?=movmean(Sw.p?_iPh_ts.data, days(15), ''omitnan'', ''SamplePoints'', t?);')
  c_eval('iPh_movm?=movmedian(Sw.p?_iPh_ts.data(~indEclipse?), days(15), ''omitnan'', ''SamplePoints'', t?(~indEclipse?));', 1:6);
  c_eval('figure(''units'', ''normalized'', ''outerposition'', [0 0 1 1]); h?=irf_plot({Sw.p?_iPh_ts(~indEclipse?), irf.ts_scalar(Sw.p?_iPh_ts.time(~indEclipse?), iPh_movm?), Sw.p?_iPh_ts(indEclipse?)},''comp'',''linestyle'',{lineStyle{?},''-black'', ''O''});', 1:6);
  c_eval('title(h?,[''Plot created: '',nowStr,''. MMS'',scStr,'' I_{ph} vs time from sweep on probe ?.'']);', 1:6);
  c_eval('ylabel(h?,{''I_{ph}'',''[nA]''});', 1:6);
  c_eval('ylim(h?,[-55 555]);', 1:6);
  if (plotDAC)
    c_eval('hold(h?, ''on'');irf_plot(h?, p?_dac);', 1:6);
    c_eval('legend(h?,''I_{ph}, sweep p? (excl. eclipse and/or maneuvers)'', ''15 days moving median (excl. eclipse and/or maneuvers)'', ''I_{ph}, sweep p? (only eclipse and/or maneuvers)'', ''Commanded DAC current'');', 1:6);
    c_eval('set(h?.Children(2), ''LineWidth'', 2);', 1:6);
  else
    c_eval('legend(h?,''I_{ph}, sweep p? (excl. eclipse and/or maneuvers)'', ''15 days moving median (excl. eclipse and/or maneuvers)'', ''I_{ph}, sweep p? (only eclipse and/or maneuvers)'');', 1:6);
  end
  c_eval('set(h?.Children(1), ''LineWidth'', 2);', 1:6);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/SDP/iPhVsTime_mms'',scStr,''_p?.png'']);', 1:4);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/ADP/iPhVsTime_mms'',scStr,''_p?.png'']);', 5:6);
  close all % close plots

  % Iph_knee vs time
  c_eval('iPh_knee_movm?=movmedian(Sw.p?_iPh_knee_ts.data(~indEclipse?), days(15), ''omitnan'', ''SamplePoints'', t?(~indEclipse?));', 1:6);
  c_eval('figure(''units'', ''normalized'', ''outerposition'', [0 0 1 1]); h?=irf_plot({Sw.p?_iPh_knee_ts(~indEclipse?), irf.ts_scalar(Sw.p?_iPh_knee_ts.time(~indEclipse?), iPh_knee_movm?), Sw.p?_iPh_knee_ts(indEclipse?)},''comp'',''linestyle'',{lineStyle{?},''-black'', ''O''});', 1:6);
  c_eval('title(h?,[''Plot created: '',nowStr,''. MMS'',scStr,'' I_{ph}_{,}_{knee} vs time from sweep on probe ?.'']);', 1:6);
  c_eval('ylabel(h?,{''I_{ph}_{,}_{knee}'',''[nA]''});', 1:6);
  c_eval('ylim(h?,[-55 555]);', 1:6);
  if (plotDAC)
    c_eval('hold(h?, ''on'');irf_plot(h?, p?_dac);', 1:6);
    c_eval('legend(h?,''I_{ph}_{,}_{knee}, sweep p? (excl. eclipse and/or maneuvers)'', ''15 days moving median (excl. eclipse and/or maneuvers)'', ''I_{ph}_{,}_{knee}, sweep p? (only eclipse and/or maneuvers)'', ''Commanded DAC currents'');', 1:6);
    c_eval('set(h?.Children(2), ''LineWidth'', 2);', 1:6);
  else
    c_eval('legend(h?,''I_{ph}_{,}_{knee}, sweep p? (excl. eclipse and/or maneuvers)'', ''15 days moving median (excl. eclipse and/or maneuvers)'', ''I_{ph}_{,}_{knee}, sweep p? (only eclipse and/or maneuvers)'');', 1:6);
  end
  c_eval('set(h?.Children(1), ''LineWidth'', 2);', 1:6);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/SDP/iPhKneeVsTime_mms'',scStr,''_p?.png'']);', 1:4);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/ADP/iPhKneeVsTime_mms'',scStr,''_p?.png'']);', 5:6);
  close all % close plots

  % "Optimal_bias" vs time
  c_eval('iPh_optimalBias_movm?=movmedian(Sw.p?_optimal_bias_ts.data(~indEclipse?), days(15), ''omitnan'', ''SamplePoints'', t?(~indEclipse?));', 1:6);
  c_eval('figure(''units'', ''normalized'', ''outerposition'', [0 0 1 1]); h?=irf_plot({Sw.p?_optimal_bias_ts(~indEclipse?), irf.ts_scalar(Sw.p?_optimal_bias_ts.time(~indEclipse?), iPh_optimalBias_movm?), Sw.p?_optimal_bias_ts(indEclipse?)},''comp'',''linestyle'',{lineStyle{?},''-black'', ''O''});', 1:6);
  c_eval('title(h?,[''Plot created: '',nowStr,''. MMS'',scStr,'' I_{optimal bias} vs time from sweep on probe ?.'']);', 1:6);
  c_eval('ylabel(h?,{''I_optimal_bias'',''[nA]''});', 1:6);
  c_eval('ylim(h?,[-55 555]);', 1:6);
  if (plotDAC)
    c_eval('hold(h?, ''on'');irf_plot(h?, p?_dac);', 1:6);
    c_eval('legend(h?,''I_{optimal bias}, sweep p? (excl. eclipse and/or maneuvers)'', ''15 days moving median (excl. eclipse and/or maneuvers)'', ''I_{optimal bias}, sweep p? (only eclipse and/or maneuvers)'', ''Commanded DAC currents'');', 1:6);
    c_eval('set(h?.Children(2), ''LineWidth'', 2);', 1:6);
  else
    c_eval('legend(h?,''I_{optimal bias}, sweep p? (excl. eclipse and/or maneuvers)'', ''15 days moving median (excl. eclipse and/or maneuvers)'', ''I_{optimal bias}, sweep p? (only eclipse and/or maneuvers)'');', 1:6);
  end
  c_eval('set(h?.Children(1), ''LineWidth'', 2);', 1:6);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/SDP/iOptimalVsTime_mms'',scStr,''_p?.png'']);', 1:4);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/ADP/iOptimalVsTime_mms'',scStr,''_p?.png'']);', 5:6);
  close all % close plots

  % impedance vs time
  c_eval('imp_movm?=movmedian(Sw.p?_impedance_ts.data(~indEclipse?), days(15),''omitnan'',''SamplePoints'',t?(~indEclipse?));', 1:6);
  c_eval('figure(''units'', ''normalized'', ''outerposition'', [0 0 1 1]); h?=irf_plot({Sw.p?_impedance_ts(~indEclipse?), irf.ts_scalar(Sw.p?_iPh_ts.time(~indEclipse?), imp_movm?), Sw.p?_impedance_ts(indEclipse?)},''comp'',''linestyle'',{lineStyle{?},''-black'', ''O''});', 1:6);
  c_eval('ylabel(h?,{''Impedance'',''[MOhm]''});', 1:6);
  c_eval('title(h?,[''Plot created: '',nowStr,''. MMS'',scStr,'' impedance vs time from sweep on probe ?.'']);', 1:6);
  c_eval('legend(h?,''impedance from P? Sweep (excl. eclipse and/or maneuvers)'',''15 days moving median (excl. eclipse and/or maneuvers)'',''impedance from P? Sweep (only eclipse and/or maneuvers)'');', 1:6);
  c_eval('ylim(h?,[-5 55]);', 1:6);
  c_eval('set(h?.Children(1),''LineWidth'',2);', 1:6);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/SDP/ImpedanceVsTime_mms'',scStr,''_p?.png'']);', 1:4);
  c_eval('print(h?.Parent, ''-dpng'', [sweepFolder,''summary_plots/ADP/ImpedanceVsTime_mms'',scStr,''_p?.png'']);', 5:6);
  close all % close plots

  % phase vs time
  figure('units','normalized','outerposition', [0 0 1 1]);
  h=irf_plot({Sw.p1_phase_ts, Sw.p2_phase_ts, Sw.p3_phase_ts, Sw.p4_phase_ts}, ...
    'comp', 'linestyle', lineStyle(1:4));
  ylabel(h,{'Phase, computed from DefAtt','[deg]'});
  title(h, ['Plot created: ',nowStr,'. MMS',scStr,' phase vs time.']);
  legend(h, 'p1', 'p2', 'p3', 'p4');
  print(h.Parent, '-dpng', [sweepFolder,'summary_plots/SDP/PhaseVsTime_mms',scStr,'_p1234.png']);
  close all % close plots
  figure('units','normalized','outerposition', [0 0 1 1]);
  h=irf_plot({Sw.p5_phase_ts, Sw.p6_phase_ts}, ...
    'comp', 'linestyle', lineStyle(5:6));
  ylabel(h,{'Phase, computed from DefAtt','[deg]'});
  title(h, ['Plot created: ',nowStr,'. MMS',scStr,' phase vs time.']);
  legend(h, 'p5', 'p6');
  print(h.Parent, '-dpng', [sweepFolder,'summary_plots/ADP/PhaseVsTime_mms',scStr,'_p56.png']);
  close all % close plots
  % phase of knee vs time
  figure('units','normalized','outerposition', [0 0 1 1]);
  h=irf_plot({Sw.p1_phase_knee_ts, Sw.p2_phase_knee_ts, Sw.p3_phase_knee_ts, Sw.p4_phase_knee_ts}, ...
    'comp', 'linestyle', lineStyle(1:4));
  ylabel(h,{'Phase at knee, computed from DefAtt','[deg]'});
  title(h, ['Plot created: ',nowStr,'. MMS',scStr,' phase_{knee} vs time.']);
  legend(h, 'p1', 'p2', 'p3', 'p4');
  print(h.Parent, '-dpng', [sweepFolder,'summary_plots/SDP/PhaseKneeVsTime_mms',scStr,'_p1234.png']);
  close all % close plots
  figure('units','normalized','outerposition', [0 0 1 1]);
  h=irf_plot({Sw.p5_phase_knee_ts, Sw.p6_phase_knee_ts}, ...
    'comp', 'linestyle', lineStyle(5:6));
  ylabel(h,{'Phase at knee, computed from DefAtt','[deg]'});
  title(h, ['Plot created: ',nowStr,'. MMS',scStr,' phase_{knee} vs time.']);
  legend(h, 'p5', 'p6');
  print(h.Parent, '-dpng', [sweepFolder,'summary_plots/ADP/PhaseKneeVsTime_mms',scStr,'_p56.png']);
  close all % close plots

  % Do not use probes after probe failure when computing how the photo
  % current or impedance depends on phase (or how it did depend on phase
  % before the probe failure). Keeping the values up to this point helps
  % to highlight the failure in the plots of iPh and impedance vs time.
  if(iSc == 4) % MMS4 p4 failed
    Sw.p4_iPh_ts.data(Sw.p4_iPh_ts.time >= EpochTT('2016-06-12T05:28:48.200Z')) = NaN;
    Sw.p4_iPh_knee_ts.data(Sw.p4_iPh_ts.time >= EpochTT('2016-06-12T05:28:48.200Z')) = NaN;
    Sw.p4_impedance_ts.data(Sw.p4_impedance_ts.time >= EpochTT('2016-06-12T05:28:48.200Z')) = NaN;
  elseif(iSc == 2) % MMS2 p2 failed
    Sw.p2_iPh_ts.data(Sw.p2_iPh_ts.time >= EpochTT('2018-09-21T06:04:45.810Z')) = NaN;
    Sw.p2_iPh_knee_ts.data(Sw.p2_iPh_ts.time >= EpochTT('2018-09-21T06:04:45.810Z')) = NaN;
    Sw.p2_impedance_ts.data(Sw.p2_impedance_ts.time >= EpochTT('2018-09-21T06:04:45.810Z')) = NaN;
  elseif(iSc == 3) % MMS2 p2 failed
    Sw.p2_iPh_ts.data(Sw.p2_iPh_ts.time >= EpochTT('2022-05-08T23:18:44.300Z')) = NaN;
    Sw.p2_iPh_knee_ts.data(Sw.p2_iPh_ts.time >= EpochTT('2022-05-08T23:18:44.300Z')) = NaN;
    Sw.p2_impedance_ts.data(Sw.p2_impedance_ts.time >= EpochTT('2022-05-08T23:18:44.300Z')) = NaN;
  end

  % iPh vs phase
  %c_eval('ind?=Sw.p?_iPh_ts.data>50&Sw.p?_iPh_ts.data<500;', 1:4);
  c_eval('ind?=Sw.p?_iPh_ts.data>0;', 1:6);
  nPhaseSegm = 360 / 5; % Number of phase segments
  c_eval('phaseInd?=discretize(Sw.p?_phase_ts.data,nPhaseSegm);', 1:6); % Index of which phase segment data was measured in
  c_eval('iPhvsPhaseSegm?=zeros(nPhaseSegm,1);',1:6);
  c_eval('iPhvsPhaseSegm!(?)=median(Sw.p!_iPh_ts.data(phaseInd!==? & ind!));', 1:nPhaseSegm, 1:6);
  c_eval('fig?=figure(''units'',''normalized'',''outerposition'',[0 0 1 1]);subplot(1,2,1);polarplot(deg2rad(Sw.p?_phase_ts.data(ind?)),Sw.p?_iPh_ts.data(ind?),lineStyle{?});rlim([0 600]);title(''All times, where I_{ph}>0.'');subplot(1,2,2);polarplot(deg2rad(360/(2*nPhaseSegm)+(0:360/nPhaseSegm:360)),[iPhvsPhaseSegm?; iPhvsPhaseSegm?(1)],''-black'',''LineWidth'',2);rlim([0 600]);title(''Median I_{ph}, where I_{ph}>0, over 5 degrees of phase.'');rticks([0 100 200 300 400 500 600]);suptitle([''Plot created: '',nowStr,''. MMS'',scStr,'' I_{ph} vs phase, p?.'']);', 1:6);
  %c_eval('fig?=figure(''units'',''normalized'',''outerposition'',[0 0 1 1]);subplot(1,2,1);polarplot(deg2rad(Sw.p?_phase_ts.data(ind?)),Sw.p?_iPh_ts.data(ind?),lineStyle{?});rlim([0 600]);title(''All times, where 0<I_{ph}.'');subplot(1,2,2);polarplot(deg2rad(2.5+0:360/nPhaseSegm:(360-1)),iPhvsPhaseSegm?,''-black'',''LineWidth'',2);rlim([0 600]);title(''Median I_{ph}, where 0<I_{ph}, over 5 degrees of phase.'');suptitle([''MMS'',scStr,'' I_{ph} vs phase, p?.'']);', 5:6);
  c_eval('print(fig?, ''-dpng'', [sweepFolder,''summary_plots/SDP/iPhVsPhasePolar_mms'',scStr,''_p?.png'']);', 1:4);
  c_eval('print(fig?, ''-dpng'', [sweepFolder,''summary_plots/ADP/iPhVsPhasePolar_mms'',scStr,''_p?.png'']);', 5:6);
  close all % close plots

  % iPh_knee vs phase_knee
  %c_eval('ind?=Sw.p?_iPh_ts.data>50&Sw.p?_iPh_ts.data<500;', 1:4);
  c_eval('ind?=Sw.p?_iPh_knee_ts.data>0;', 1:6);
  nPhaseKneeSegm = 360 / 5; % Number of phase segments
  c_eval('phaseKneeInd?=discretize(Sw.p?_phase_knee_ts.data,nPhaseKneeSegm);', 1:6); % Index of which phase segment data was measured in
  c_eval('iPhKneevsPhaseKneeSegm?=zeros(nPhaseKneeSegm,1);',1:6);
  c_eval('iPhKneevsPhaseKneeSegm!(?)=median(Sw.p!_iPh_knee_ts.data(phaseKneeInd!==? & ind!));', 1:nPhaseKneeSegm, 1:6);
  c_eval('fig?=figure(''units'',''normalized'',''outerposition'',[0 0 1 1]);subplot(1,2,1);polarplot(deg2rad(Sw.p?_phase_knee_ts.data(ind?)),Sw.p?_iPh_knee_ts.data(ind?),lineStyle{?});rlim([0 600]);title(''All times, where I_{ph}_{,}_{knee}>0.'');subplot(1,2,2);polarplot(deg2rad(360/(2*nPhaseKneeSegm)+(0:360/nPhaseKneeSegm:360)),[iPhKneevsPhaseKneeSegm?; iPhKneevsPhaseKneeSegm?(1)],''-black'',''LineWidth'',2);rlim([0 600]);title(''Median I_{ph}_{,}_{knee}, where I_{ph}_{,}_{knee}>0, over 5 degrees of phase.'');rticks([0 100 200 300 400 500 600]);suptitle([''Plot created: '',nowStr,''. MMS'',scStr,'' I_{ph}_{,}_{knee} vs phase of knee, p?.'']);', 1:6);
  %c_eval('fig?=figure(''units'',''normalized'',''outerposition'',[0 0 1 1]);subplot(1,2,1);polarplot(deg2rad(Sw.p?_phase_ts.data(ind?)),Sw.p?_iPh_ts.data(ind?),lineStyle{?});rlim([0 600]);title(''All times, where 0<I_{ph}.'');subplot(1,2,2);polarplot(deg2rad(2.5+0:360/nPhaseSegm:(360-1)),iPhvsPhaseSegm?,''-black'',''LineWidth'',2);rlim([0 600]);title(''Median I_{ph}, where 0<I_{ph}, over 5 degrees of phase.'');suptitle([''MMS'',scStr,'' I_{ph} vs phase, p?.'']);', 5:6);
  c_eval('print(fig?, ''-dpng'', [sweepFolder,''summary_plots/SDP/iPhKneeVsPhasePolar_mms'',scStr,''_p?.png'']);', 1:4);
  c_eval('print(fig?, ''-dpng'', [sweepFolder,''summary_plots/ADP/iPhKneeVsPhasePolar_mms'',scStr,''_p?.png'']);', 5:6);
  close all % close plots

  % impedance vs phase
  c_eval('ind?=Sw.p?_impedance_ts.data>0;', 1:6);
  nPhaseSegm = 360 / 5; % Number of phase segments
  c_eval('phaseInd?=discretize(Sw.p?_phase_ts.data,nPhaseSegm);', 1:6); % Index of which phase segment data was measured in
  c_eval('ImpvsPhaseSegm?=zeros(nPhaseSegm,1);',1:6);
  c_eval('ImpvsPhaseSegm!(?)=median(Sw.p!_impedance_ts.data(phaseInd!==? & ind!));', 1:nPhaseSegm, 1:6);
  c_eval('fig?=figure(''units'',''normalized'',''outerposition'',[0 0 1 1]);subplot(1,2,1);polarplot(deg2rad(Sw.p?_phase_ts.data(ind?)),Sw.p?_impedance_ts.data(ind?),lineStyle{?});rlim([0 55]);title(''All times, where Impedance>0.'');subplot(1,2,2);polarplot(deg2rad(360/(2*nPhaseSegm)+(0:360/nPhaseSegm:360)),[ImpvsPhaseSegm?; ImpvsPhaseSegm?(1)],''-black'',''LineWidth'',2);rlim([0 55]);title(''Median Impedance, where Impedance>0, over 5 degrees of phase.'');rticks([0 10 20 30 40 50]);suptitle([''Plot created: '',nowStr,''. MMS'',scStr,'' Impedance vs phase, p?.'']);', 1:6);
  c_eval('print(fig?, ''-dpng'', [sweepFolder,''summary_plots/SDP/ImpedanceVsPhasePolar_mms'',scStr,''_p?.png'']);', 1:4);
  c_eval('print(fig?, ''-dpng'', [sweepFolder,''summary_plots/ADP/ImpedanceVsPhasePolar_mms'',scStr,''_p?.png'']);', 5:6);
  close all % close plots

  %    % iPh vs F10.7
  %    % Resample F10.7 cm flux to time of the sweeps
  %    c_eval('f107_resamp?_ts=resample(f107_ts,Sw.p?_impedance_ts,''linear'');', 1:4);
  %    fig = figure('units','normalized','outerposition',[0 0 1 1]);
  %    h=irf_plot({[f107_resamp1_ts.data, Sw.p1_iPh_ts.data], ...
  %      [f107_resamp2_ts.data, Sw.p2_iPh_ts.data], ...
  %      [f107_resamp3_ts.data, Sw.p3_iPh_ts.data], ...
  %      [f107_resamp4_ts.data, Sw.p4_iPh_ts.data]}, ...
  %      'comp','linestyle', lineStyle(1:4));
  %    xlabel(h,{'F10.7','[s.f.u.]'});
  %    ylabel(h,{'Iph','[nA]'});
  %    legend(h,'p1','p2','p3','p4');
  %    title(h,['Plot created: ',nowStr,'. MMS',scStr,' I_{ph} vs F10.7 flux.']);
  %    print(fig, '-dpng', [sweepFolder,'summary_plots/iPhVsF107_mms',scStr,'.png']);
  %    close all % close plots

end % for iSc
end
