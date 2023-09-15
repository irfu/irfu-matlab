% Script to run mms.probe_align_times and optional plots of fk-power spectra
% for all times when the probes satisfy alignment conditions
% Written by D. B. Graham


mms.db_init('local_file_db','/data/mms');

% Select spacecraft number: 1--4
ic = 3;

% Load data from MMS database. Kinda slow; much faster to load from a local
% directory.
if 1
  %load background magnetic field data
  tint = irf.tint('2015-09-02T15:20:00.00Z/2015-09-02T15:30:00.00Z');
  c_eval('Bxyz=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);
  %load spacecraft phase - zphase
  c_eval('zphase = mms.db_get_variable(''mms?_ancillary_defatt'',''zphase'',tint);',ic);
  zphase = irf.ts_scalar(zphase.time, zphase.zphase);
  %load burst mode interval
  tints = irf.tint('2015-09-02T15:26:00.00Z/2015-09-02T15:28:00.00Z');
  c_eval('Exyz=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tints);',ic);
  c_eval('SCpot=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv'',tints);',ic);
  Exyz.data(abs(Exyz.data) > 100) = NaN;
  SCpot.data(abs(SCpot.data) > 100) = NaN;
end

%load brst mode interval from file
%tmpDataObj = dataobj(irf_ssub('data/mms?_edp_brst_ql_dce2d_20150902152504_v0.2.0.cdf',ic));
%Exyz = get_variable(tmpDataObj,irf_ssub('mms?_edp_dce_xyz_dsl',ic));
%Exyz = mms.variable2ts(Exyz);
%Exyz.data(find(abs(Exyz.data) > 100)) = NaN;

%tmpDataObj = dataobj(irf_ssub('data/mms?_edp_brst_l2_scpot_20150902152504_v0.2.0.cdf',ic));
%SCpot = get_variable(tmpDataObj,irf_ssub('mms?_edp_dcv',ic));
%SCpot = mms.variable2ts(SCpot);
%SCpot.data(find(abs(SCpot.data) > 100)) = NaN;

[starttime1,endtime1,starttime3,endtime3] = mms.probe_align_times(Exyz,Bxyz,SCpot,zphase,1);

if 0 % Set to 1 to plot all fk power spectra
  probecomb = 1; %#ok<UNRCH>

  for ii=1:length(starttime1)

    if (endtime1(ii)-starttime1(ii) > 0.2)
      tint = irf.tint(strcat(starttime1(ii).utc,'/',endtime1(ii).utc));
      [fkpower,freq,wavenumber] = mms.fk_powerspectrum(SCpot, Bxyz, zphase, tint, ic, probecomb);

      h=irf_plot(1,'newfigure');
      h(1)=irf_panel('disprel');
      pcolor(h(1),wavenumber,freq*1e-3,log10(fkpower))
      shading(h(1),'flat')
      ylabel(h(1),'f (kHz)','fontsize',20);
      xlabel(h(1),'k_{||} (m^{-1})','fontsize',20);
      set(h(1),'fontsize',20)

      %Overplot phase speed estimate
      vph = mms.estimate_phase_speed(fkpower,freq,wavenumber);
      kfit = 0.0001:0.0001:0.1;
      ffit = abs(vph/(2*pi)*kfit);
      if(vph < 0)
        kfit = -kfit;
      end
      hold(h(1),'on');
      plot(h(1),kfit,ffit*1e-3,'linewidth',3,'color','r')
      hold(h(1),'off');
      irf_legend(h(1),strcat('v = ',num2str(vph(1)/1000),'km s^{-1}'),[0.6 0.9],'fontsize',20,'color','r')

      c=colorbar;
      ylabel(c,'log_{10} P(f,k)/P_{max}','fontsize',20);
      set(c,'fontsize',20)
      set(h(1),'position',[0.09 0.2 0.76 0.78]);

      pictitle = strcat('MMS',num2str(ic),'p',num2str(probecomb),num2str(probecomb+1),'t',starttime1(ii).utc,'to',endtime1(ii).utc,'.eps');
      set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
      print(pictitle,'-depsc','-painters');

    end
  end

  probecomb = 3;

  for ii=1:length(starttime3)

    if (endtime3(ii)-starttime3(ii) > 0.2)
      tint = irf.tint(strcat(starttime3(ii).utc,'/',endtime3(ii).utc));
      [fkpower,freq,wavenumber] = mms.fk_powerspectrum(SCpot, Bxyz, zphase, tint, ic, probecomb);

      h=irf_plot(1,'newfigure');
      h(1)=irf_panel('disprel');
      pcolor(h(1),wavenumber,freq*1e-3,log10(fkpower))
      shading(h(1),'flat')
      ylabel(h(1),'f (kHz)','fontsize',20);
      xlabel(h(1),'k_{||} (m^{-1})','fontsize',20);
      set(h(1),'fontsize',20)

      %Overplot phase speed estimate
      vph = mms.estimate_phase_speed(fkpower,freq,wavenumber);
      kfit = 0.0001:0.0001:0.1;
      ffit = abs(vph/(2*pi)*kfit);
      if(vph < 0)
        kfit = -kfit;
      end
      hold(h(1),'on');
      plot(h(1),kfit,ffit*1e-3,'linewidth',3,'color','r')
      hold(h(1),'off');
      irf_legend(h(1),strcat('v = ',num2str(vph(1)/1000),'km s^{-1}'),[0.6 0.9],'fontsize',20,'color','r')

      c=colorbar;
      ylabel(c,'log_{10} P(f,k)/P_{max}','fontsize',20);
      set(c,'fontsize',20)
      set(h(1),'position',[0.09 0.2 0.76 0.78]);

      pictitle = strcat('MMS',num2str(ic),'p',num2str(probecomb),num2str(probecomb+1),'t',starttime3(ii).utc,'to',endtime3(ii).utc,'.eps');
      set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
      print(pictitle,'-depsc','-painters');

    end
  end
end