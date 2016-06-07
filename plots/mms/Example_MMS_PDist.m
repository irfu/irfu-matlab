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

%% Example operations
ePdist1_omni = ePDist1.omni; % omnidirectional differential energy flux
ePDist1_pa = g; % construt pitchangle distribution
ePDist1_lowE = ePDist1.elim([0 200]); % limit energy range
ePDist1_deflux = ePDist1.deflux; % change units to differential energy flux
ePDist1_e64 = ePDist1.e64; % resample energy to 64 energy levels, reduces the time resolution

%% Example plot: omnidirectional differential energy flux
h = irf_plot(4);

hca = irf_panel('e omni');
irf_spectrogram(hca,ePDist2.omni('e').specrec,'log')
hca.YScale = 'log';
hca.YLim = [10 30000];

hca = irf_panel('e omni high');
irf_spectrogram(hca,ePDist2.omni('e').elim([3000 40000]).specrec,'log')
hca.YScale = 'log';
hca.YLim = [10 30000];

hca = irf_panel('e omni mid');
irf_spectrogram(hca,ePDist2.omni('e').elim([200 3000]).specrec,'log')
hca.YScale = 'log';
hca.YLim = [10 30000];

hca = irf_panel('e omni low');
irf_spectrogram(hca,ePDist2.omni('e').elim([0 200]).specrec,'log')
hca.YScale = 'log';
hca.YLim = [10 30000];

%% Pitchangles
ePitch = ePDist2.pitchangles(dmpaB2,[]);
ePitchDEF = ePitch.deflux; % change units

%%  Example plot: Pitchangles for different energies
h = irf_plot(5);

hca = irf_panel('e omni high');
irf_spectrogram(hca,ePDist2.omni('e').specrec,'log')
hca.YScale = 'log';
hca.YTick = [1e1 1e2 1e3 1e4];

hca = irf_panel('e pa all');
irf_spectrogram(hca,ePitch.deflux.specrec('pa'))

hca = irf_panel('e pa high');
eint = [0 500];
irf_spectrogram(hca,ePitch.deflux.elim(eint).specrec('pa'))
irf_legend(hca,[num2str(eint(1),'%g') '-' num2str(eint(2),'%g') ' eV'],[0.02 0.95],'color',[1 1 1])

hca = irf_panel('e pa low');
eint = [500 30000];
irf_spectrogram(hca,ePitch.deflux.elim(eint).specrec('pa'))
irf_legend(hca,[num2str(eint(1),'%g') '-' num2str(eint(2),'%g') ' eV'],[0.02 0.95],'color',[1 1 1])

hca = irf_panel('e pa low2');
eint = [0 20];
irf_spectrogram(hca,ePitch.deflux.elim(eint).specrec('pa'))
irf_legend(hca,[num2str(eint(1),'%g') '-' num2str(eint(2),'%g') ' eV'],[0.02 0.95],'color',[1 1 1])

