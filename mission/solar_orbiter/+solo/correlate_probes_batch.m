% solo.correlate_probes_batch
%
% Script to get basic calibration coefficients for BIAS E and SCPOT
%
%% Times for probe potentials discontinuities
discontTimes=EpochTT(solo.ProbePotDiscontinuities);

%% Specify output calibration file name
strname = 'd23K123_20211021';
%% Load data
% One month, but we take +/- 3 days as a margin to better match with the
% nearby month.
margin = 3*24*60*60; %seconds.
Tint = irf.tint('2020-03-01T00:00:00Z/2021-09-30T23:59:59.99Z')+[-1,1]*margin;

% If there is a discontinuity in the data, e.g. potential jumps due to the
% solar panels, as in late 2020, early 2021, generate subintervals and
% apply the analysis on them separately.
sub_int_times = EpochTT(solo.split_tint(Tint,discontTimes));

%Pre-define output variables:
DCE_SRF=irf.ts_vec_xyz(EpochTT([]),double.empty(0,3));
DCE_SRF_10s=irf.ts_vec_xyz(EpochTT([]),double.empty(0,3));
PSP=irf.ts_scalar(EpochTT([]),[]);
PSP_10s=irf.ts_scalar(EpochTT([]),[]);

% Generate calibration file for each subinterval separately
for isub=1:length(sub_int_times)-1
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
    
    %%
    h=irf_plot(1,'newfigure');
    irf_plot(h,OutTS_step1.d23,'.-')
    hold on
    irf_plot(h,OutTS_step1.k23,'.-')
    irf_plot(h,OutTS_step1.d123,'.-');
    irf_plot(h,OutTS_step1.k123,'.-');
    irf_plot(h,OutTS_step1.del123,'.-');
    
    irf_zoom(h,'y');
    
    legend(h,'d23 (Volts)', 'k23', 'd123 (Volts)','k123, V1=k123*V23+del123','del123')
    title(h,'Step 1')
    %%
    d23 = irf.ts_scalar(OutTS_step1.d23.time,movmedian(OutTS_step1.d23.data,5,'omitnan','Endpoints','fill'));
    k23 = irf.ts_scalar(OutTS_step1.k23.time,movmedian(OutTS_step1.k23.data,5,'omitnan','Endpoints','fill'));
    %% Step 2 - get k123
    OutTS_step2 = solo.correlate_probes(VDC,EDC,d23,k23);
    
    
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
    
    k123 = irf.ts_scalar(OutTS_step1.k123.time,movmedian(OutTS_step2.k123.data,5,'omitnan','Endpoints','fill'));
    
    %% Step 3 - get d123
    OutTS_step3 = solo.correlate_probes(VDC,EDC,d23,k23,k123);
    
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
    
    d123 = irf.ts_scalar(OutTS_step1.d123.time,movmedian(OutTS_step3.del123.data,5,'omitnan','Endpoints','fill'));
    
    %% Save result into a monthly file
    K123 = irf.ts_scalar(OutTS_step1.d12.time,[k123.data, d123.data]);
    filenamestr = [strname,'_subint',num2str(isub),'_of_',num2str(length(sub_int_times)-1)];
    save(filenamestr, 'K123','k23', 'd23');
  if 1
    continue, %skip validation...
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
    
end


%% Combine the calibration files of the subintervals (this can be done in a better way.)
if 1
  
tmpstr1 = [strname,'_subint1_of_5'];
tmpstr2 = [strname,'_subint2_of_5'];
tmpstr3 = [strname,'_subint3_of_5'];
tmpstr4 = [strname,'_subint4_of_5'];
tmpstr5 = [strname,'_subint5_of_5'];

load(tmpstr1);
a=load(tmpstr2);


K123_new = K123.combine(a.K123);
d23_new = d23.combine(a.d23);
k23_new = k23.combine(a.k23);

a=load(tmpstr3);
K123 = K123_new.combine(a.K123);
d23 = d23_new.combine(a.d23);
k23 = k23_new.combine(a.k23);

a=load(tmpstr4);
K123 = K123.combine(a.K123);
d23 = d23.combine(a.d23);
k23 = k23.combine(a.k23);

a=load(tmpstr5);
K123 = K123.combine(a.K123);
d23 = d23.combine(a.d23);
k23 = k23.combine(a.k23);

save(strname, 'K123', 'k23','d23');

end


%% The code below can be used to compare the "piece-wise" calibration using
% the subinterval calibration files done above, with the calibration
% obtained using the combined calibration files. The two should be exactly
% the same.
if 0
    [DCE_SRF_comb,PSP_comb,~] = solo.vdccal(VDC,[],'d23K123_20210511_combined_parts');
    DCE_SRF_10s=DCE_SRF_comb.resample(outTime);
    PSP_10s = PSP_comb.resample(outTime);
end
%% Validation plot 2 - reuires B from SOAR
B_RTN = solo.db_get_ts('solo_L2_mag-rtn-normal-1-minute','B_RTN', Tint);
% B_SRF = solo.db_get_ts('solo_L2_mag-srf-normal','B_SRF', Tint);
B_SRF = B_RTN ;
B_SRF.data(:,1:2) = -B_SRF.data(:,1:2);

E_sw = irf_e_vxb([-250 0 0],B_SRF);

h = irf_figure(948274,4,'reset');

hca = irf_panel('B');
hl = irf_plot(hca,B_SRF)
ylabel(hca,'B [nT]')

%hca = irf_panel('Ex');
%hl = irf_plot(hca,DCE_SRF.x);
%hold(hca,'on')
%hl = irf_plot(hca,E_sw.x,'r');
%hold(hca,'off')
%ylabel(hca,'Ex SRF [mV/m]')
%legend(hca,'E','-Vsw x B')

hca = irf_panel('Ey');
hl = irf_plot(hca,DCE_SRF_10s.y);
hold(hca,'on')
hl = irf_plot(hca,E_sw.y,'r');
hold(hca,'off')
ylabel(hca,'Ey SRF [mV/m]')

hca = irf_panel('Ez');
hl = irf_plot(hca,DCE_SRF_10s.z);
hold(hca,'on')
hl = irf_plot(hca,E_sw.z,'r');
hold(hca,'off')
ylabel(hca,'Ez SRF [mV/m]')

set(h(2:3), 'YLim', 4.99*[-1 1])


hca = irf_panel('PSP');
hl = irf_plot(hca,PSP_10s);
ylabel(hca,'PSP [V]')

irf_zoom(h,'x',Tint)


%% Load all months and combine together
% (This is only needed if a the calibration procedure has changed,
% otherwise see sections below)
mm = {'january','february','march','april','may','june','july','august','september','october'};
K123 = []; d23 = []; CC = []; Gamma0 = []; Gamma1 = [];
for m=3:10
    dd = load(['d23K123_20201110_' mm{m} '.mat']);
    tStart = EpochTT(sprintf('2020-%02d-01T00:00:00Z',m));
    tStop = EpochTT(sprintf('2020-%02d-01T00:00:00Z',m+1));
    if m==3, tStart = dd.d23.time(1);
    elseif m==10, tStop = dd.d23.time(end);
    end
    tCC=dd.CC.tlim(irf.tint(tStart,tStop));
    tGamma0=dd.Gamma0.tlim(irf.tint(tStart,tStop));
    tGamma1=dd.Gamma1.tlim(irf.tint(tStart,tStop));
    td23 = dd.d23.tlim(irf.tint(tStart,tStop));
    tK123 = dd.K123.tlim(irf.tint(tStart,tStop));
    if isempty(d23)
        d23 = td23;
        K123 = tK123;
        CC = tCC;
        Gamma0 = tGamma0;
        Gamma1 = tGamma1;
    else
        d23 = d23.combine(td23);
        K123 = K123.combine(tK123);
        CC = CC.combine(tCC);
        Gamma0 = Gamma0.combine(tGamma0);
        Gamma1 = Gamma1.combine(tGamma1);
    end
end

h = irf_figure(948278,1,'reset');
irf_plot(K123);
hold on
irf_plot(d23);
legend('k123, V1=k123*V23+del123','del123 (Volts)', 'd23 (Volts)')

filename = 'd23K123_YYYYMMDD';
save(filename, 'K123', 'd23', 'Gamma0', 'Gamma1', 'CC');

%% Append one month manually
% This is just an adaptation of the above code section. Instead of
% combining many months at the same time, we take them one at a time. This
% is more convenient when we just append to the old calibration file.

%% ====
% Alternative 1: If we need to replace parts in the old calibration file we run
% the following. (This was needed in Dec. 2020, when we improved the code
% to take into account V discontinuities due to solar panels.
if 0
    % Load the calibration file we want to append variables to.
    old_cal = load('d23K123_20210129');
    oCC = old_cal.CC;
    oGamma0 = old_cal.Gamma0;
    oGamma1 = old_cal.Gamma1;
    oK123 = old_cal.K123;
    od23 = old_cal.d23;
    
    % Load the calibration file we want to append to the old file.
    new_cal = load('d23K123_20210511_december');
    tStart = '2020-12-01T00:00:00.00Z';
    tStop = '2021-01-01T00:00:00.00Z';
    nTint=irf.tint(tStart,tStop);
    nCC = new_cal.CC.tlim(nTint);
    nGamma0 = new_cal.Gamma0.tlim(nTint);
    nGamma1 = new_cal.Gamma1.tlim(nTint);
    nK123 = new_cal.K123.tlim(nTint);
    nd23 = new_cal.d23.tlim(nTint);
    
    % tend = the time at which we want to append the new calibration.
    tend = EpochTT('2020-11-30T23:00:00.000000000Z');
    % Cut the old TimeSeries
    CC = oCC.tlim(irf.tint(oCC.time(1),tend));
    Gamma0 = oGamma0.tlim(irf.tint(oGamma0.time(1),tend));
    Gamma1 = oGamma1.tlim(irf.tint(oGamma1.time(1),tend));
    K123 = oK123.tlim(irf.tint(oK123.time(1),tend));
    d23 = od23.tlim(irf.tint(od23.time(1),tend));
    
    %Append the new ones.
    CC=CC.combine(nCC);
    Gamma0=Gamma0.combine(nGamma0);
    Gamma1=Gamma1.combine(nGamma1);
    K123=K123.combine(nK123);
    d23=d23.combine(nd23);
    
    %Save data
    filename = 'd23K123_20210521';
    save(filename, 'K123', 'd23', 'Gamma0', 'Gamma1', 'CC');
end

%% ====
% Alternative 2: Append the new data to the calibration file without
% cutting the old one. (Should be the standard way)
if 1
    
    % Load the calibration file we want to append variables to.
    old_cal = load('d23K123_20210521_v2');
    oCC = old_cal.CC;
    oGamma0 = old_cal.Gamma0;
    oGamma1 = old_cal.Gamma1;
    oK123 = old_cal.K123;
    od23 = old_cal.d23;
    
    % Load the calibration file we want to append to the old file.
    new_cal = load('d23K123_20210519_february');
    % Specify time interval of the new data:
    tStart = '2021-02-01T00:00:00.00Z';
    tStop = '2021-03-01T00:00:00.00Z';
    nTint=irf.tint(tStart,tStop);
    
    nCC = new_cal.CC.tlim(nTint);
    nGamma0 = new_cal.Gamma0.tlim(nTint);
    nGamma1 = new_cal.Gamma1.tlim(nTint);
    nK123 = new_cal.K123.tlim(nTint);
    nd23 = new_cal.d23.tlim(nTint);
    
    %Append the new ones.
    CC=oCC.combine(nCC);
    Gamma0=oGamma0.combine(nGamma0);
    Gamma1=oGamma1.combine(nGamma1);
    K123=oK123.combine(nK123);
    d23=od23.combine(nd23);
    
    %Save data
    newfilename = 'd23K123_20210521_v3';
    save(newfilename, 'K123', 'd23', 'Gamma0', 'Gamma1', 'CC');
end

%% BIAS change times
biasUpdateTimes = [...
    '2020-05-18T04:05:55Z';...
    '2020-06-14T00:00:00Z';...
    '2020-06-22T00:00:00Z';...
    '2020-07-06T00:00:00Z';...
    '2020-07-14T00:00:00Z';...
    '2020-07-20T07:39:41Z';...
    '2020-08-11T21:27:03Z';...
    '2020-09-05T00:00:00Z';...
    '2020-09-14T00:00:00Z'];

irf_plot(irf.ts_scalar(EpochTT(biasUpdateTimes),zeros(size(biasUpdateTimes,1),1)),'*')
legend('k123, V1=k123*V23+del123','del123 (Volts)', 'd23 (Volts)', 'Bias current update')
