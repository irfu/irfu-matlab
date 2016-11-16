%%   HPCA brust mode summary plot
%       1. Bxyz
%       2. numder density [10 second resolution]
%       3-7. dFlux [H+, He+, He++, O+, O++]
%   NOTES: 'tistring?': time string for HPCA/ion data; 'tmstring?': time string for HPCA/moments dat;
%   The upper two variable are slightly different for one crossing and between four spacecraft.

%   History:
%       1. 2016-02-19, v0
%       2. 2016-03-21, v1
%       3. 2016-11-01, fix mean0 bug

%%  1. basic
    ic = 1;
    Tint  = irf.tint('2015-11-12T06:19:07Z/2015-11-12T06:20:00Z');
    c_eval('hpcai_dr? = ''/Volumes/WYLIMMS/work/data/mms?/hpca/brst/l1b/ion/2015/11/12/'';', ic);
    c_eval('hpcam_dr? = ''/Volumes/WYLIMMS/work/data/mms?/hpca/brst/l1b/moments/2015/11/12/'';', ic);
    c_eval('tistring? = ''061901'';', ic);
    c_eval('tmstring? = ''061907'';', ic);
    c_eval('hpcai_fn? = [''mms?_hpca_brst_l1b_ion_20151112'' tistring? ''_v0.2.1.cdf''];', ic);
    c_eval('hpcam_fn? = [''mms?_hpca_brst_l1b_moments_20151112'' tmstring? ''_v0.2.1.cdf''];', ic);
%   FGM
    c_eval('fgm_dr? = ''/Volumes/WYLIMMS/work/data/mms?/fgm/srvy/l2/2015/11/'';', ic);
    c_eval('fgm_fn? = ''mms?_fgm_srvy_l2_20151112_v4.18.0.cdf'';', ic);  
    
%%  2. get data
%       2.1. HPCA moments
    c_eval('hpcam_obj? = dataobj([hpcam_dr? hpcam_fn?]);', ic);
    c_eval('NHp? = get_ts(hpcam_obj?, ''mms?_hpca_hplus_number_density'');', ic);
    c_eval('NHep? = get_ts(hpcam_obj?, ''mms?_hpca_heplus_number_density'');', ic);
    c_eval('NHepp? = get_ts(hpcam_obj?, ''mms?_hpca_heplusplus_number_density'');', ic);
    c_eval('NOp? = get_ts(hpcam_obj?, ''mms?_hpca_oplus_number_density'');', ic);
    c_eval('NOpp? = get_ts(hpcam_obj?, ''mms?_hpca_oplusplus_number_density'');', ic);    
    c_eval('VHp? = get_ts(hpcam_obj?, ''mms?_hpca_hplus_ion_bulk_velocity'');', ic);
    c_eval('VHep? = get_ts(hpcam_obj?, ''mms?_hpca_heplus_ion_bulk_velocity'');', ic);
    c_eval('VHepp? = get_ts(hpcam_obj?, ''mms?_hpca_heplusplus_ion_bulk_velocity'');', ic);
    c_eval('VOp? = get_ts(hpcam_obj?, ''mms?_hpca_oplus_ion_bulk_velocity'');', ic);
    c_eval('VOpp? = get_ts(hpcam_obj?, ''mms?_hpca_oplusplus_ion_bulk_velocity'');', ic);    
%       2.2. HPCA ion
    c_eval('hpcai_obj? = dataobj([hpcai_dr? hpcai_fn?]);', ic);
    c_eval('HpfluxV? = get_variable(hpcai_obj?,''mms?_hpca_hplus_flux'');', ic);
    c_eval('HpfluxT? = get_ts(hpcai_obj?,''mms?_hpca_hplus_flux'');', ic);
    c_eval('HepfluxV? = get_variable(hpcai_obj?,''mms?_hpca_heplus_flux'');', ic);
    c_eval('HepfluxT? = get_ts(hpcai_obj?,''mms?_hpca_heplus_flux'');', ic);
    c_eval('HeppfluxV? = get_variable(hpcai_obj?,''mms?_hpca_heplusplus_flux'');', ic);
    c_eval('HeppfluxT? = get_ts(hpcai_obj?,''mms?_hpca_heplusplus_flux'');', ic);    
    c_eval('OpfluxV? = get_variable(hpcai_obj?,''mms?_hpca_oplus_flux'');', ic);
    c_eval('OpfluxT? = get_ts(hpcai_obj?,''mms?_hpca_oplus_flux'');', ic);
    c_eval('OppfluxV? = get_variable(hpcai_obj?,''mms?_hpca_oplusplus_flux'');', ic);
    c_eval('OppfluxT? = get_ts(hpcai_obj?,''mms?_hpca_oplusplus_flux'');', ic);  
    c_eval('energyHp = HpfluxV?.DEPEND_2.data;', ic);
    c_eval('energyHep = HepfluxV?.DEPEND_2.data;', ic);
    c_eval('energyHepp = HeppfluxV?.DEPEND_2.data;', ic);
    c_eval('energyOp = OpfluxV?.DEPEND_2.data;', ic);
    c_eval('energyOpp = OppfluxV?.DEPEND_2.data;', ic);
%       2.3. B    
    c_eval('fgm_obj? = dataobj([fgm_dr? fgm_fn?]);', ic);
    c_eval('Bxyz? = get_ts(fgm_obj?, ''mms?_fgm_b_gse_srvy_l2'');', ic);

%%  3. compute
    % 3.0. fix mean0 error.
    c_eval('HpfluxT?.data(HpfluxT?.data<=0) = NaN;', ic);    
    c_eval('HepfluxT?.data(HepfluxT?.data<=0) = NaN;', ic);
    c_eval('HeppfluxT?.data(HeppfluxT?.data<=0) = NaN;', ic);
    c_eval('OpfluxT?.data(OpfluxT?.data<=0) = NaN;', ic);
%       3.1. HpfluxT1
    c_eval('specHp=struct(''t'',HpfluxT?.time.epochUnix);', ic);
    c_eval('specHp.p = squeeze(irf.nanmean(HpfluxT?.data, 2));', ic);
    specHp.p_label={'flux','[1/cc s sr eV]'};
    specHp.f_label={'Energy', '[eV]'};
    specHp.f = single(energyHp);
%       3.2. HepfluxT1
    c_eval('specHep=struct(''t'',HepfluxT?.time.epochUnix);', ic);
    c_eval('specHep.p = squeeze(irf.nanmean(HepfluxT?.data, 2));', ic);
    specHep.p_label={'flux','[1/cc s sr eV]'};
    specHep.f_label={'Energy', '[eV]'};
    specHep.f = single(energyHep);
%       3.3. HeppfluxT1
    c_eval('specHepp=struct(''t'',HeppfluxT?.time.epochUnix);', ic);
    c_eval('specHepp.p = squeeze(irf.nanmean(HeppfluxT?.data, 2));', ic);
    specHepp.p_label={'flux','[1/cc s sr eV]'};
    specHepp.f_label={'Energy', '[eV]'};
    specHepp.f = single(energyHepp);
%       3.4. OpfluxT1
    c_eval('specOp=struct(''t'',OpfluxT?.time.epochUnix);', ic);
    c_eval('specOp.p = squeeze(irf.nanmean(OpfluxT?.data, 2));', ic);
    specOp.p_label={'flux','[1/cc s sr eV]'};
    specOp.f_label={'Energy', '[eV]'};
    specOp.f = single(energyOp);    
    specHepp.f = single(energyHepp);
%       3.5. OppfluxT1
    c_eval('specOpp=struct(''t'',OppfluxT?.time.epochUnix);', ic);
    c_eval('specOpp.p = squeeze(irf.nanmean(OppfluxT?.data, 2));', ic);
    specOpp.p_label={'flux','[1/cc s sr eV]'};
    specOpp.f_label={'Energy', '[eV]'};
    specOpp.f = single(energyOpp);   
    
%%  4. plot
    h=irf_plot(7, 'newfigure');
%       4.1. Bxyz
    h(1)=irf_panel('Bxyz');
    c_eval('irf_plot(h(1), Bxyz?);', ic);
    ylabel(h(1), {'B_{GSE}', '[nT]'},'Interpreter','tex');
    irf_legend(h(1),{'B_{X}','B_{Y}','B_{Z}'},[0.5 0.1])
    irf_zoom(h(1),'y',[-30 50]);

%       4.2. Ni
    h(2)=irf_panel('n');
    c_eval('irf_plot(h(2), NHp?);', ic);
    hold(h(2),'on');
    c_eval('irf_plot(h(2), NHep?);', ic);
    c_eval('irf_plot(h(2), NHepp?);', ic);
    c_eval('irf_plot(h(2), NOp?);', ic);
    c_eval('irf_plot(h(2), NOpp?);', ic);
    hold(h(2),'off');
    set(h(2),'yscale','log');
    irf_zoom(h(2),'y',[0.00001 5]);
    ylabel(h(2),{'N', '[cm^{-3}]'},'Interpreter','tex');
    irf_legend(h(2),{'H+','He+', 'He++', 'O+', 'O++'},[0.5 0.9]);

%       4.3. Hp flux    
    h(3)=irf_panel('Hp');
    irf_spectrogram(h(3),specHp,'log');
    set(h(3),'yscale','log');
    set(h(3),'ytick',[1e1 1e2 1e3 1e4]);
    %caxis(h(5),[-15 -5])
    irf_zoom(h(3),'y',[0.9 5e4]);    
    ylabel(h(3),{'E_{H+}', '[eV]'},'Interpreter','tex');

%       4.4. Hep flux    
    h(4)=irf_panel('Hep');
    irf_spectrogram(h(4),specHep,'log');
    set(h(4),'yscale','log');
    set(h(4),'ytick',[1e1 1e2 1e3 1e4]);
    %caxis(h(4),[10^0.5 10^2])
    irf_zoom(h(4),'y',[0.9 5e4]);        
    ylabel(h(4),{'E_{He+}', '[eV]'},'Interpreter','tex');
 
%       4.5. Hepp flux    
    h(5)=irf_panel('Hepp');
    irf_spectrogram(h(5),specHepp,'log');
    set(h(5),'yscale','log');
    set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
    %caxis(h(5),[-15 -5])
    irf_zoom(h(5),'y',[0.9 5e4]);        
    ylabel(h(5),{'E_{He++}', '[eV]'},'Interpreter','tex');

%       4.6. Op flux    
    h(6)=irf_panel('Op');
    irf_spectrogram(h(6),specOp,'log');
    set(h(6),'yscale','log');
    set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
    %caxis(h(5),[-15 -5])
    irf_zoom(h(6),'y',[0.9 5e4]);    
    ylabel(h(6),{'E_{O+}', '[eV]'},'Interpreter','tex');
 
%       4.7. Opp flux    
    h(7)=irf_panel('Opp');
    irf_spectrogram(h(7),specOpp,'log');
    set(h(7),'yscale','log');
    set(h(7),'ytick',[1e1 1e2 1e3 1e4]);
    %caxis(h(5),[-15 -5])
    irf_zoom(h(7),'y',[0.9 5e4]);        
    ylabel(h(7),{'E_{O++}', '[eV]'},'Interpreter','tex');    
%       4.X.
    title(h(1),strcat('MMS',num2str(ic)))
    irf_plot_axis_align(h(1:7));
    irf_zoom(h(1:7),'x', Tint);
    
%%