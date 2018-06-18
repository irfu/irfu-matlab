figure(hcf); % hcf -- figure handle in which to plot
plotId='ttbrowse_plot_1';
numOfPanels=5; 
if diff(tint)==0, return; end % if zero time interval return
hcfTag = get(hcf,'Tag');
if isempty(hcfTag)
	h=irf_plot(numOfPanels,'reset');
	ud.subplot_handles=h;
	set(hcf,'tag',plotId);
	set(hcf,'userdata',ud);
else
	ud=get(hcf,'userdata');
	h=ud.subplot_handles;
end
% tintiso='2003-09-17T08:00:00Z/2003-09-17T10:00:00Z';
% tintiso='2003-09-24T15:30:00Z/2003-09-24T16:30:00Z';
% tint=irf_time(tintiso,'utc>tint');
% h=irf_plot(4,'newfigure');
ic=4;
pwd
% get FGM data in GSM
irf_log('fcal','---- Reading FGM data')
c_eval('B?=local.c_read(''B_vec_xyz_gse__C?_CP_FGM_5VPS'',tint);');
c_eval('B?=irf_abs(B?);');
c_eval('gsmB?=irf_gse2gsm(B?);');
irf_log('fcal','---- Reading RAPID data')
c_eval('[JRAP,dobj_JRAP]=local.c_read(''Electron_Dif_flux__C?_CP_RAP_ESPCT6'',tint,''caa'');',ic);
irf_log('fcal','---- Reading PEACE data')
c_eval('[JPEACE,dobj_JPEACE]=local.c_read(''Data__C?_CP_PEA_PITCH_SPIN_DPFlux'',tint,''caa'');',ic);
irf_log('fcal','---- Reading CIS data')
c_eval('VCIS=local.c_read(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',tint);',ic);
c_eval('gsmVCIS=irf_gse2gsm(VCIS);',ic);

% plot
irf_zoom(h,'x',tint);
zoom reset;
if 1   % PANEL: C?       FGM Bx,By,Bz,B GSM
	hca=irf_panel(h,'C? FGM B GSM');
	c_eval('irf_plot(hca,gsmB?);',ic);
	irf_zoom(hca,'y');
	ylabel(hca,'B [nT] GSM');
	irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.1])
	irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
dd=irf_time(tint(1),'epoch>yyyy-mm-dd hh:mm:ss');
title(hca,dd(1:10));
if 1   % PANEL: C?       CIS Vx,Vy,Vz,V CODIF(HIA) GSM
  hca=irf_panel('C? CIS V GSM');
  if ic ~=2 % on s/c 2 there is no CIS
    irf_plot(hca,gsmVCIS);
    % c_eval('irf_plot(hca,gsmVCIS?);',ic); % HIA 
    ylabel(hca,'V_p [km/s] GSM');
    irf_zoom(hca,'y');
    irf_legend(hca,{'V_X','V_Y','V_Z'},[0.02 0.98])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
  end
end
if 1   % PANEL: C1..C4   FGM BZ GSM
	hca=irf_panel(h,'PANEL: C1..C4, FGM BZ');
	c_pl_tx(hca,'gsmB?',4)
	irf_zoom(hca,'y');
	ylabel(hca,'B_Z [nT] GSM');
	irf_legend(hca,{'C1','C2','C3','C4'},[0.98 0.98],'color','cluster');
end
if 1   % PANEL: C?       RAPID electron spectrogram
	hca=irf_panel(h,'C?_CP_RAP_ESPCT6');
	%varname=irf_ssub('Electron_Dif_flux__C?_CP_RAP_ESPCT6',ic);
	%varunits=irf_get_data(varname,'caa','unit');
	varunits={'log_{10} dF','1/cm^2 s sr keV'};
	plot(hca,dobj_JRAP,JRAP,'colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
	caxis(hca,[0.01 4.49]);
	ylabel(hca,'E [keV]');
	set(hca,'yscale','log');
	set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
	irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
if 1   % PANEL: C?       PEACE PITCH_SPIN_DEFlux spectrogram omni
    hca=irf_panel('C? PEACE energy spectra');
    %varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    %varunits=irf_get_data(varname,'caa','units');
    varunits={'log_{10} dEF','keV/cm^2 s sr keV'};
    plot(hca,dobj_JPEACE,JPEACE,'sum_dim1','colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
    caxis(hca,[5.6 7.9]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E [eV]');
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.2],'color','k')
end


irf_plot_axis_align(h);
irf_timeaxis(h);
r=local.c_read('r',tint);
% add position in GSM
c_eval('add_position(h(end),irf_gse2gsm(r.R?));',ic)
irf_timeaxis(h(end),'nodate');
