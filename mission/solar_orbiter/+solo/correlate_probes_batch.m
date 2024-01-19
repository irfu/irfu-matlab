% solo.correlate_probes_batch
%
% Script to get basic calibration coefficients for BIAS E and SCPOT


% Times for probe potentials discontinuities (solar panel fuses)
discontTimes=EpochTT(solo.ProbePotDiscontinuities);

% Specify output calibration file name
strname = 'd23K123_20230707_test';

% Select time interval. We take +/- 3 days as a margin to better match with the
% old calibration file if we only do a partial update.
calTint = irf.tint('2022-12-01T00:00:00Z/2023-05-28T23:59:59.99Z');
margin = 3*24*60*60; %seconds.
Tint=calTint+[-1,1]*margin;

% If there is a discontinuity in the data, e.g. potential jumps due to the
% solar panels, as in late 2020, early 2021, generate subintervals and
% apply the analysis on them separately.
sub_int_times = EpochTT(solo.split_tint(Tint,discontTimes));
NrOfSubints = length(sub_int_times)-1;

% If you want to plot the calibrated E and PSP, change to =1.
% Not recommended if you are calibrating a large amount of data.
validate_cal_plot=0;

% If you want to plot all the calibration parameters set to = 1.
% (recommended). Good for validation/debugging purposes.
% (Not a problem for large amounts of data)
cal_param_plot = 1;

% If you want to plot the different outputs during the calibration stage,
% set to = 1; (Not a problem for large amounts of data)
stage_plots = 0;

%Pre-define output variables:
DCE_SRF=irf.ts_vec_xyz(EpochTT([]),double.empty(0,3));
DCE_SRF_10s=irf.ts_vec_xyz(EpochTT([]),double.empty(0,3));
PSP=irf.ts_scalar(EpochTT([]),[]);
PSP_10s=irf.ts_scalar(EpochTT([]),[]);

% Generate calibration file for each subinterval separately
for isub=1:length(sub_int_times)-1
  % Load data for the current subinterval
  subTint=sub_int_times(isub:isub+1);
  VDC = solo.db_get_ts('solo_L2_rpw-lfr-surv-cwf-e-1-second', 'VDC', subTint);
  EDC = solo.db_get_ts('solo_L2_rpw-lfr-surv-cwf-e-1-second', 'EDC', subTint);
  QUAL = solo.db_get_ts('solo_L2_rpw-lfr-surv-cwf-e-1-second', 'QUALITY_FLAG', subTint);

  % Remover thrusters, etc
  iiBad = QUAL.data <2;
  EDC.data(iiBad,:) = NaN;
  VDC.data(iiBad,:) = NaN;

  %% Step 1 - get d23
  OutTS_step1 = solo.correlate_probes(VDC, EDC);

  if stage_plots
    h=irf_plot(1,'newfigure');
    irf_plot(h,OutTS_step1.d23,'.-') %Offset between 2,3
    hold on
    irf_plot(h,OutTS_step1.k23,'.-') %Slope between 2,3
    irf_plot(h,OutTS_step1.d123,'.-');
    irf_plot(h,OutTS_step1.k123,'.-');
    irf_plot(h,OutTS_step1.del123,'.-');

    irf_zoom(h,'y');
    legend(h,'d23 (Volts)', 'k23', 'd123 (Volts)','k123, V1=k123*V23+del123','del123')
    title(h,'Step 1')
  end
  % Smooth the data
  d23 = irf.ts_scalar(OutTS_step1.d23.time,movmedian(OutTS_step1.d23.data,5,'omitnan','Endpoints','fill'));
  k23 = irf.ts_scalar(OutTS_step1.k23.time,movmedian(OutTS_step1.k23.data,5,'omitnan','Endpoints','fill'));


  %% Step 2 - get k123
  OutTS_step2 = solo.correlate_probes(VDC,EDC,d23,k23);

  if stage_plots
    h=irf_plot(1,'newfigure');
    irf_plot(h(1),OutTS_step2.d23,'.-')
    hold on
    irf_plot(h(1),OutTS_step2.k23,'.-')
    irf_plot(h(1),OutTS_step2.d123,'.-');
    irf_plot(h(1),OutTS_step2.k123,'.-');
    irf_plot(h(1),OutTS_step2.del123,'.-');
    irf_zoom(h(1),'y');
    legend(h(1),'d23 (Volts)','k23', 'd123 (Volts)','k123, V1=k123*V23+del123','del123')
    title(h(1),'Step 2')
  end

  % Smooth the data
  k123 = irf.ts_scalar(OutTS_step1.k123.time,movmedian(OutTS_step2.k123.data,5,'omitnan','Endpoints','fill'));

  %% Step 3 - get d123
  OutTS_step3 = solo.correlate_probes(VDC,EDC,d23,k23,k123);

  if stage_plots
    h=irf_plot(1,'newfigure');
    irf_plot(h(1),OutTS_step3.d23,'.-')
    hold on
    irf_plot(h(1),OutTS_step3.k23,'.-')
    irf_plot(h(1),OutTS_step3.d123,'.-');
    irf_plot(h(1),OutTS_step3.k123,'.-');
    irf_plot(h(1),OutTS_step3.del123,'.-');

    irf_zoom(h(1),'y');
    legend(h(1),'d23 (Volts)','k23', 'd123 (Volts)','k123, V1=k123*V23+del123','del123')
    title(h(1),'Step 3')
  end

  % Smooth the data
  d123 = irf.ts_scalar(OutTS_step1.d123.time,movmedian(OutTS_step3.del123.data,5,'omitnan','Endpoints','fill'));

  %% Save results for the given subinterval
  K123 = irf.ts_scalar(OutTS_step1.d12.time,[k123.data, d123.data]);
  filenamestr = [strname,'_subint',num2str(isub),'_of_',num2str(length(sub_int_times)-1)];
  save(filenamestr, 'K123','k23', 'd23');

  if validate_cal_plot
    %% Validation plot 1
    % downsample to 10 sec resolution to make plotting easier
    Tstart = VDC.time.start+3600;
    TstartS = Tstart.toUtc;
    TstartS = [TstartS(1:11) '00:00:00Z'];
    Tstart = EpochTT(TstartS);

    Tstop = VDC.time.stop+3600;
    TstopS = Tstop.toUtc;
    TstopS = [TstopS(1:11) '00:00:00Z'];
    Tstop = EpochTT(TstopS);

    DT = 10; nSteps =  (Tstop-Tstart)/DT;
    outTime = Tstart + ((1:nSteps)-0.5)*DT;
    %% Apply calibration on data in the interval
    [DCE_SRF_temp,PSP_temp,SCPOT_temp] = solo.vdccal(VDC,[],filenamestr);
    DCE_SRF_10s_temp = DCE_SRF_temp.resample(outTime);
    PSP_10s_temp = PSP_temp.resample(outTime);

    %Combine the different parts.
    DCE_SRF=DCE_SRF.combine(DCE_SRF_temp);
    DCE_SRF_10s=DCE_SRF_10s.combine(DCE_SRF_10s_temp);
    PSP=PSP.combine(PSP_temp);
    PSP_10s=PSP_10s.combine(PSP_10s_temp);

    h = irf_figure(948273,4,'reset');
    irf_plot({DCE_SRF_10s,k123,d123,PSP_10s})
  end

end % Exit sub interval loop


%% Combine the calibration files of the subintervals

c_eval('tmpstrings{?} = [strname,''_subint?_of_'',num2str(NrOfSubints)];',1:NrOfSubints);
load(tmpstrings{1});
if NrOfSubints > 1
  for ii=1:length(tmpstrings)
    a=load(tmpstrings{ii});
    K123 = K123.combine(a.K123);
    d23 = d23.combine(a.d23);
    k23 = k23.combine(a.k23);
  end
end
save(strname, 'K123', 'k23','d23');

%% Plot the resulting calibration parameters.
if cal_param_plot
  Tint=[d23.time.start;d23.time.stop];

  BLL = solo.get_data('LL_B_RTN',Tint);
  VLL = solo.get_data('LL_V_RTN',Tint);

  % load bias currents (gives BIAS_current TSeries)
  load('/Volumes/solo/BJOR_files/Bias_currents.mat');

  posSolO = solo.get_position(Tint,'frame','ECLIPJ2000');
  [radius, lon, lat] = cspice_reclat(posSolO.data');
  solopos = irf.ts_vec_xyz(posSolO.time,[radius',lon',lat']);
  AU=1.4960e+11*10^-3; % km

  discrete_events = EpochTT(['2020-12-21T01:39:23.24Z';...
    '2021-01-10T23:33:05.37Z';'2021-01-17T00:19:26.99Z';...
    '2021-03-01T00:45:20.39Z';'2021-08-17T16:57:38.00Z';...
    '2021-08-17T17:16:30.00Z';'2021-08-17T17:52:53.00Z']);

  d123=irf.ts_scalar(K123.time,K123.data(:,2));
  k123=irf.ts_scalar(K123.time,K123.data(:,1));

  h=irf_plot(6,'newfigure');
  fig=h.Parent;
  fig.Position=[10 10 2500 1000];


  irf_plot(h(1),d23,'-o','MarkerFaceColor','k','MarkerSize',2);
  ylabel(h(1),'d23 [V]','interpreter','tex');


  irf_plot(h(2),k23,'-o','MarkerFaceColor','k','MarkerSize',2);
  ylabel(h(2),'k23','interpreter','tex');

  irf_plot(h(3),d123,'-o','MarkerFaceColor','k','MarkerSize',2);
  ylabel(h(3),'d123 [V]','interpreter','tex');

  irf_plot(h(4),k123,'-o','MarkerFaceColor','k','MarkerSize',2);
  ylabel(h(4),'k123','interpreter','tex');

  irf_plot(h(5),BLL)
  ylabel(h(5),'B [nT] LL ','interpreter','tex');

  irf_plot(h(6),VLL)
  ylabel(h(6),'V [km/s] LL ','interpreter','tex');

  irf_zoom(h(1:6),'x',Tint);

  h(1).YLim=[-2,2];
  h(2).YLim=[0.5 1.5];
  h(3).YLim=[-5,0.5];
  h(4).YLim=[0.3,1.8];
  h(5).YLim=[-60,60];
  h(6).YLim=[-1000 500];

  hold(h(1),'on');
  text(h(1),0.70,1.1,'Discrete events','fontsize',26,'units','normalized','color',[1,0.65,0]);
  text(h(1),0.85,1.1,'Bias changes','fontsize',26,'units','normalized','color',[0,0,1]);


  hold(h(2),'on');
  hold(h(3),'on');
  hold(h(4),'on');

  yyaxis(h(3),'right');
  irf_plot(h(3),solopos.x/AU,'r');
  h(3).YColor=[1,0,0];
  ylabel(h(3),'R [AU]','interpreter','tex');
  yyaxis(h(3),'left');
  for ii=1:length(discrete_events)
    markTint = discrete_events(ii)+[-10,10]*60*60;
    irf_pl_mark(h(1:6),markTint,[1,0.65,0]);
  end

  for ii=1:length(BIAS_current)
    markTint = BIAS_current.time(ii)+[-10,10]*30*60;
    irf_pl_mark(h(1:6),markTint,[0,0,1]);
  end
  irf_timeaxis(h(4), 'nodate' );
end


%% Combine calfiles
% If you wish to combine the new calibration file with the old one run this
% code section after filling in the name of the old cal file.

if 0
  oldfile = 'd23K123_20220124'; % Needs to be updated
  newfile = strname;
  combinedfile = 'd23K123_YYMMDD'; %Needs to be specified

  oldcal = load(oldfile);
  newcal = load(strname);

  K123 = oldcal.K123.combine(newcal.K123);
  d23 = oldcal.d23.combine(newcal.d23);
  k23 = oldcal.k23.combine(newcal.k23);

  save(combinedfile, 'K123', 'k23','d23');
end



