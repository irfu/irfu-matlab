%Example to plot SolO-PAS ion differential energy flux
%Contributed by Sergio Toledo-Redondo on 2021-09-21

%% Set time
Tint = irf.tint('2020-08-02T07:00:00.00Z/2020-08-02T08:59:59.00Z');

%% Load data
B_SRF = solo.db_get_ts('solo_L2_mag-srf-normal','B_SRF', Tint);
ieflux = solo.db_get_ts('solo_L2_swa-pas-eflux','eflux',Tint);
myFile=solo.db_list_files('solo_L2_swa-pas-eflux',Tint);
iEnergy = cdfread([myFile.path '/' myFile.name],'variables','Energy');
iEnergy = iEnergy{1};

%% Plot
h=irf_plot(2,'newfigure');
%
hca=irf_panel('BXYZ');
irf_plot(hca,B_SRF);
ylabel(hca,'B [nT] SRF');
irf_legend(hca,{'B_X';'B_Y';'B_Z'},[1.02 0.98]);
%
hca = irf_panel('ieflux');
specrec   = struct('t', ieflux.time.epochUnix);
specrec.p = ieflux.data;
specrec.p_label='dEF';
specrec.f = repmat(iEnergy,1,numel(specrec.t))';
irf_spectrogram(hca,specrec,'log');
set(hca, 'YScale', 'log');
ylabel(hca,'[eV]');

axis(h,'tight')
irf_plot_axis_align(h)