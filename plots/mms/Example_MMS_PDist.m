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
ic = 1; % spacecraft number
time = irf_time('2015-10-16T10:33:00.00Z','utc>epochtt');
%[filepath,filename] = mms.get_filepath('mms1_fpi_brst_l2_des-dist',time);
filepath_and_filename = mms.get_filepath('mms1_fpi_brst_l2_des-dist',time);
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
ePdist1_omni = ePDist1.omni; % omnidirectional differential energy flux
ePDist1_pa = ePDist1.pitchangles(dmpaB1,24); % construt pitchangle distribution
ePDist1_lowE = ePDist1.elim([0 200]); % limit energy range
ePDist1_deflux = ePDist1.deflux; % change units to differential energy flux
ePDist1_dpflux = ePDist1.dpflux; % change units to particle energy flux
ePDist1_skm = ePDist1.convertto('s^3/km^6'); % change units to particle energy flux
ePDist1_e64 = ePDist1.e64; % resample energy to 64 energy levels, reduces the time resolution
ePDist_specrec = ePDist1.specrec; % change format to specrec, to be used as input to irf_plot or irf_spectrogram

%% Example plots: time series
nPanels = 11;
h = irf_plot(nPanels);

if 1
  hca = irf_panel('B');
  irf_plot(hca,dmpaB1)
  hca.YLabel.String = 'B (nT)';
end
if 1
  hca = irf_panel('e omni');
  irf_spectrogram(hca,ePDist1.omni.specrec)
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
end
if 1
  hca = irf_panel('e omni 64 energy channels');
  irf_spectrogram(hca,ePDist1.e64.omni.specrec)
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
end
if 1
  hca = irf_panel('e pitchangles');
  irf_spectrogram(hca,ePDist1.pitchangles(dmpaB1,18).specrec('pa')); 
  hca.YTick = [0 45 90 135];  
end
if 1
  elim = [10 200];
  hca = irf_panel('e pitchangles low');
  irf_spectrogram(hca,ePDist1.elim(elim).pitchangles(dmpaB1,18).specrec('pa')); 
  hca.YTick = [0 45 90 135];  
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  elim = [10 200];
  hca = irf_panel('e pitchangles 64 energy channels');
  irf_spectrogram(hca,ePDist1.e64.elim(elim).pitchangles(dmpaB1,18).specrec('pa')); 
  hca.YTick = [0 45 90 135];  
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  elim = [200 400];
  hca = irf_panel('e pitchangles mid');
  irf_spectrogram(hca,ePDist1.elim(elim).pitchangles(dmpaB1,18).specrec('pa')); 
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  elim = [400 20000];
  hca = irf_panel('e pitchangles high');
  irf_spectrogram(hca,ePDist1.elim(elim).pitchangles(dmpaB1,18).specrec('pa')); 
  hca.YTick = [0 45 90 135];  
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  palim = [0 15];
  hca = irf_panel('e spectrogram parallel');
  irf_spectrogram(hca,ePDist1.pitchangles(dmpaB1,[0 15]).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  palim = [75 105];
  hca = irf_panel('e spectrogram perpendicular');
  irf_spectrogram(hca,ePDist1.pitchangles(dmpaB1,palim).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1
  palim = [165 180];
  hca = irf_panel('e spectrogram anti-parallel');
  irf_spectrogram(hca,ePDist1.pitchangles(dmpaB1,[75 105]).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h,'x',tint)
irf_plot_axis_align

%% Example plots: particle distributions
ic = 1;
time = irf_time('2015-10-16T10:33:30.00Z','utc>epochtt');

ePitch1 = ePDist1.pitchangles(dmpaB1,[17]);
c_eval('B0 = dmpaB?.resample(time).data;',ic); 
hatB0 = double(irf_norm(B0));
c_eval('E0 = dslE?.resample(time).data;',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);
c_eval('scpot = scPot?.resample(time);',ic); 

% optional input parameters for projection plot
vlim = 12*1e3; % x and ylim
elevlim = 15; % angle over plane to include in slice
strCMap = 'jet'; % colormap
projclim = [0 5]; % colorbar limit

x = hatE0;
y = hatExB0;
z = hatB0;

c_eval('tInd = find(abs(ePDist?.time-time)==min(abs(ePDist?.time-time)));', ic);

% Initialize figure
nRows = 3; nCols = 3;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
isub = 1;

hca = h(isub); isub = isub + 1; 
xyz = [hatE0; hatExB0; hatB0 * (-1)];
vlabels = {'v_E','v_{ExB}','v_B'};
mms.plot_projection(hca,ePDist1.convertto('s^3/km^6'),'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);

hca = h(isub); isub = isub + 1; 
xyz = [hatExB0; hatB0; hatE0 * (-1)];
vlabels = {'v_{ExB}','v_B','v_E'};
mms.plot_projection(hca,ePDist1.convertto('s^3/km^6'),'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);

hca = h(isub); isub = isub + 1; 
xyz = [hatB0; hatE0; hatExB0 * (-1)];
vlabels = {'v_B','v_E','v_{ExB}'};
mms.plot_projection(hca,ePDist1.convertto('s^3/km^6'),'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);

hca = h(isub); isub = isub + 1;      
mms.plot_skymap(hca,ePDist1,'tint',time,'energy',100,'flat');

hca = h(isub); isub = isub + 1;      
mms.plot_skymap(hca,ePDist1,'tint',time,'energy',100,'flat','log');
%hca.CLim = projclim;

hca = h(isub); isub = isub + 1;      
mms.plot_skymap(hca,ePDist1,'tint',time,'energy',100,'vectors',{hatB0,'B'},'log');

hca = h(isub); isub = isub + 1; 
plot(hca,ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,1),...
         ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,9),...
         ePitch1.depend{1}(tInd,:),ePitch1.data(tInd,:,17));
hca.XScale = 'log';
hca.YScale = 'log';     
hca.YLabel.String = ['f_e (' ePitch1.units ')'];
hca.XLabel.String = 'E_e (eV)';
hca.XLim = [1e1 3e3];  
hca.YLim = [1e-32 2e-25];
hleg = irf_legend(hca,{'0';'90';'180'},[0.98 0.98]);

hca = h(isub); isub = isub + 1; 
plot(hca,ePitch1.depend{2},squeeze(ePitch1.data(tInd,:,:)));
hca.YScale = 'log';     
hca.YLabel.String = ['f_e (' ePitch1.units ')'];
hca.XLabel.String = '\theta (deg)';
hca.XLim = [0 180];  
