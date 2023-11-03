% Example showing how to you can work with particle distributions.
% PDist is a subclass of TSeries having the additional properties:
%   type - skymap, pitchangles, omnideflux
%   depend - skymap: energy, phi, theta
%            pitchangles: energy, pitchangles
%            omni - energies
%   ancillary - additional help data, can be energy0, energy0, energysteptable
%
% SUGGESTIONS WELCOME ! :)

%% Load PDist using mms.make_pdist
ic = 3; % spacecraft number
time = irf_time('2015-12-02T01:14:56.300Z','utc>epochtt');
filepath_and_filename = mms.get_filepath(irf_ssub('mms?_fpi_brst_l2_des-dist',ic),time);
c_eval('[ePDist?,ePDistError?] = mms.make_pdist(filepath_and_filename);',ic)
c_eval('[iPDist?,iPDistError?] = mms.make_pdist(filepath_and_filename);',ic)

%% Make skymap directly with PDist
c_eval('ePDist? = PDist(desDist?.time,desDist?.data,''skymap'',energy,ephi?.data,etheta?.data);',ic)
c_eval('ePDist?.userData = desDist?.userData; ePDist?.name = desDist?.name; ePDist?.units = desDist?.units;',ic)
c_eval('ePDist?.units = ''s^3/cm^6'';',ic)
c_eval('ePDist?',ic)

%% Make skymap with irf.ts_skymap
% If energy table is NOT passed, energy0, energy1 and energysteptable is necessary
c_eval('ePDist? = irf.ts_skymap(desDist?.time,desDist?.data,[],ephi?.data,etheta?.data,''energy0'',eenergy0?.data,''energy1'',eenergy1?.data,''esteptable'',estepTable?.data);',ic)
c_eval('ePDist?.units = ''s^3/cm^6'';',ic)
c_eval('ePDist?',ic)
% If energy table is passed, energy0, energy1 and energysteptable is not necessary
c_eval('ePDist? = irf.ts_skymap(desDist?.time,desDist?.data,energy,ephi?.data,etheta?.data);',ic)
c_eval('ePDist?.units = ''s^3/cm^6'';',ic)
c_eval('ePDist?',ic)

%% Load supporting data
c_eval('tint = ePDist?.time([1 end]);', ic);
c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
%dmpaB1=mms.db_get_ts('mms1_fgm_brst_l2','mms1_fgm_b_dmpa_brst_l2',tint);

%% Example operations
ePdist_omni = ePDist1.omni; % omnidirectional differential energy flux
ePDist_pa = ePDist1.pitchangles(dmpaB1,24); % construt pitchangle distribution
ePDist_lowE = ePDist1.elim([0 200]); % limit energy range
ePDist_deflux = ePDist1.deflux; % change units to differential energy flux
ePDist_dpflux = ePDist1.dpflux; % change units to particle energy flux
ePDist_skm = ePDist1.convertto('s^3/km^6'); % change units to particle energy flux
ePDist_e64 = ePDist1.e64; % resample energy to 64 energy levels, reduces the time resolution
ePDist_specrec = ePDist1.specrec; % change format to specrec, to be used as input to irf_plot or irf_spectrogram

%% Example plots: time series
nPanels = 7;
h = irf_plot(nPanels);
c_eval('ePDistN = ePDist?; dmpaBN = dmpaB?;',ic)

if 1
  hca = irf_panel('B');
  irf_plot(hca,dmpaBN)
  hold(hca,'on')
  irf_plot(hca,dmpaBN.abs)
  hca.YLabel.String = 'B (nT)';
  title(hca,irf_ssub('MMS?',ic))
end
if 1
  hca = irf_panel('e omni');
  irf_spectrogram(hca,ePDistN.omni.specrec)
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
end
if 0
  hca = irf_panel('e omni 64 energy channels'); %#ok<UNRCH>
  irf_spectrogram(hca,ePDistN.e64.omni.specrec)
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
end
if 0
  hca = irf_panel('e pitchangles'); %#ok<UNRCH>
  irf_spectrogram(hca,ePDistN.pitchangles(dmpaBN,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
end
if 0
  elim = [10 200]; %#ok<UNRCH>
  hca = irf_panel('e pitchangles low');
  irf_spectrogram(hca,ePDistN.elim(elim).pitchangles(dmpaBN,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  elim = [20 200];
  hca = irf_panel('e pitchangles 64 energy channels');
  irf_spectrogram(hca,ePDistN.e64.elim(elim).pitchangles(dmpaBN,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  elim = [200 2000];
  hca = irf_panel('e pitchangles mid');
  irf_spectrogram(hca,ePDistN.e64.elim(elim).pitchangles(dmpaBN,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0
  elim = [400 20000]; %#ok<UNRCH>
  hca = irf_panel('e pitchangles high');
  irf_spectrogram(hca,ePDistN.elim(elim).pitchangles(dmpaBN,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  palim = [0 15];
  hca = irf_panel('e spectrogram parallel');
  irf_spectrogram(hca,ePDistN.pitchangles(dmpaBN,[0 15]).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  palim = [75 105];
  hca = irf_panel('e spectrogram perpendicular');
  irf_spectrogram(hca,ePDistN.pitchangles(dmpaBN,palim).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  palim = [165 180];
  hca = irf_panel('e spectrogram anti-parallel');
  irf_spectrogram(hca,ePDistN.pitchangles(dmpaBN,palim).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h,'x',tint)
irf_plot_axis_align

%% Example plots: particle distributions
ic = 3;
t0 = irf_time('2015-12-02T01:14:54.320Z','utc>epochtt');
c_eval('tInd = find(abs(ePDist?.time-t0)==min(abs(ePDist?.time-t0)));', ic);
c_eval('B0 = dmpaB?.resample(ePDist?.time).data;',ic);
c_eval('E0 = dslE?.resample(ePDist?.time).data;',ic);
c_eval('scpot = scPot?.resample(ePDist?.time);',ic);
c_eval('ePitchN = ePDist?.pitchangles(dmpaB?,[17]);',ic)
c_eval('ePDistN = ePDist?;',ic)
hatB0 = double(irf_norm(B0));
hatE0tmp = double(irf_norm(E0));


% optional input parameters for projection plot
vlim = 12*1e3; % x and ylim
elevlim = 15; % angle over plane to include in slice
strCMap = 'jet'; % colormap
projclim = [0 5]; % colorbar limit


%%
for i=1:80
  idx = tInd + i -1;
  time = ePDistN.time(idx);

  hatExB0 = cross(hatE0tmp(idx,:),hatB0(idx,:));
  hatE0 = cross(hatB0(idx,:),hatExB0);

  x = hatE0;
  y = hatExB0;
  z = hatB0(idx,:);

  % Initialize figure

  % Alt 1: 2x3 plots
  %nRows = 2; nCols = 3;
  %for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
  % Alt 2: 3 + 2 plots
  for ii = 1:3; h(ii) = subplot(2,3,ii); end
  for ii = 4:5; h(ii) = subplot(2,2,ii-1); end
  isub = 1;

  hca = h(isub); isub = isub + 1;
  xyz = [x; y; z * (-1)];
  vlabels = {'v_E','v_{ExB}','v_B'};
  mms.plot_projection(hca,ePDistN.convertto('s^3/km^6'),'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);

  hca = h(isub); isub = isub + 1;
  xyz = [y; z; x * (-1)];
  vlabels = {'v_{ExB}','v_B','v_E'};
  mms.plot_projection(hca,ePDistN.convertto('s^3/km^6'),'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);

  hca = h(isub); isub = isub + 1;
  xyz = [z; x; y * (-1)];
  vlabels = {'v_B','v_E','v_{ExB}'};
  mms.plot_projection(hca,ePDistN.convertto('s^3/km^6'),'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);

  if 0
    hca = h(isub); isub = isub + 1; %#ok<UNRCH>
    mms.plot_skymap(hca,ePDistN,'tint',time,'energy',150,'flat');

    hca = h(isub); isub = isub + 1;
    mms.plot_skymap(hca,ePDistN,'tint',time,'energy',150,'flat','log');
    %hca.CLim = projclim;

    hca = h(isub); isub = isub + 1;
    mms.plot_skymap(hca,ePDistN,'tint',time,'energy',150,'vectors',{hatB0,'B'},'log');
  end

  hca = h(isub); isub = isub + 1;
  plot(hca,ePitchN.depend{1}(idx,:),ePitchN.data(idx,:,1),...
    ePitchN.depend{1}(idx,:),ePitchN.data(idx,:,9),...
    ePitchN.depend{1}(idx,:),ePitchN.data(idx,:,17));
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.YLabel.String = ['f_e (' ePitchN.units ')'];
  hca.XLabel.String = 'E_e (eV)';
  hca.XLim = [1e1 1e3];
  hca.YLim = [1e-31 2e-26];
  hleg = irf_legend(hca,{'0';'90';'180'},[0.98 0.98]);
  title(hca,irf_ssub('MMS?',ic))

  hca = h(isub); isub = isub + 1;
  plot(hca,ePitchN.depend{2},squeeze(ePitchN.data(idx,:,:)));
  hca.YScale = 'log';
  hca.YLabel.String = ['f_e (' ePitchN.units ')'];
  hca.XLabel.String = '\theta (deg)';
  hca.XLim = [0 180];
  hca.YLim = [1e-31 2e-26];
  hca.XTick = [0 45 90 135 180];

  irf_print_fig([irf_ssub('mms?_desdist_',ic) irf_fname(time,4)],'png')
end
