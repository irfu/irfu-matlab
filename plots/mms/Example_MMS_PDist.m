% Example showing how to you can work with particle distributions.
% PDist is a subclass of TSeries having the additional properties:
%   type - skymap, pitchangles, omni
%   depend - skymap: energy, phi, theta
%            pitchangles: energy, pitchangles
%            omni - energies
%   ancillary - additional help data, can be energy0, energy0, energysteptable   
% 
% SUGGESTIONS WELCOME ! :)

%% Load data WITH DOBJ
db_info = datastore('mms_db');   
ic  = 2;
disp('Loading electron distribution...')
c_eval('etmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l2/des-dist/2015/10/16/mms2_fpi_brst_l2_des-dist_20151016103254_v2.1.0.cdf'']);',ic);
c_eval('desDist? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_dist_brst''));',ic);
c_eval('eenergy0? = get_variable(etmpDataObj?,''mms?_des_energy0_brst'');',ic);
c_eval('eenergy1? = get_variable(etmpDataObj?,''mms?_des_energy1_brst'');',ic);
c_eval('ephi? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_phi_brst''));',ic);
c_eval('etheta? = get_variable(etmpDataObj?,''mms?_des_theta_brst'');',ic);
c_eval('estepTable? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_steptable_parity_brst''));',ic);

% Make energytable from energy0, energy1 and energysteptable
c_eval('energy = repmat(torow(eenergy0?.data),numel(estepTable?.data),1);',ic)
c_eval('energy(estepTable?.data==1,:) = repmat(eenergy1?.data,sum(estepTable?.data),1);',ic)

%% Make skymap directly with PDist
c_eval('ePDist? = PDist(desDist?.time,desDist?.data,''skymap'',energy,ephi?.data,etheta?.data);',ic)
c_eval('ePDist?.userData = desDist?.userData; ePDist?.name = desDist?.name; ePDist?.units = desDist?.units;',ic)
c_eval('ePDist?.units = ''s^3/m^6'';',ic)
c_eval('ePDist?',ic)

%% Make skymap with irf.ts_skymap
% If energy table is NOT passed, energy0, energy1 and energysteptable is necessary
c_eval('ePDist? = irf.ts_skymap(desDist?.time,desDist?.data,[],ephi?.data,etheta?.data,''energy0'',eenergy0?.data,''energy1'',eenergy1?.data,''esteptable'',estepTable?.data);',ic)
c_eval('ePDist?.units = ''s^3/m^6'';',ic)
c_eval('ePDist?',ic)
% If energy table is passed, energy0, energy1 and energysteptable is not necessary
c_eval('ePDist? = irf.ts_skymap(desDist?.time,desDist?.data,energy,ephi?.data,etheta?.data);',ic)
c_eval('ePDist?.units = ''s^3/m^6'';',ic)
c_eval('ePDist?',ic)

%% Make omnidirectional differential energy flux
eOMNI2 = ePDist2.omni('e');
eOMNI2low = eOMNI.elim([0 200]);

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
h = irf_plot(4);

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

