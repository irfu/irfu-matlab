%   check FPI/fast/ql data

%%  1. basic
    ic = 1;
    time = irf_time('2016-11-07T03:00:00.000Z','utc>epochtt');
    c_eval('ifile? = ''/Volumes/mms/mms?/fpi/fast/ql/dis/2016/11/mms?_fpi_fast_ql_dis_20161109200000_v3.1.0.cdf'';', ic);
    c_eval('efile? = ''/Volumes/mms/mms?/fpi/fast/ql/des/2016/11/mms?_fpi_fast_ql_des_20161109200000_v3.1.0.cdf'';', ic);
    %c_eval('ifile? = mms.get_filepath([''mms?_fpi_fast_ql_dis''], time);', ic);
    %c_eval('efile? = mms.get_filepath([''mms?_fpi_fast_ql_des''], time);', ic);
    
%%  2. load
    % 2.1. load FPI
    c_eval('iobj? = dataobj(ifile?);', ic);
    c_eval('Ni = get_ts(iobj?, ''mms?_dis_numberdensity_fast'');', ic);
    c_eval('dbcsVi = get_ts(iobj?, ''mms?_dis_bulkv_dbcs_fast'');', ic);
    c_eval('paraTi = get_ts(iobj?, ''mms?_dis_temppara_fast'');', ic);
    c_eval('perpTi = get_ts(iobj?, ''mms?_dis_tempperp_fast'');', ic);
    c_eval('enfluxi = get_ts(iobj?, ''mms?_dis_energyspectr_omni_fast'');', ic);
    c_eval('ienergy = get_ts(iobj?, ''mms?_dis_energy_fast'');', ic);    
    c_eval('eobj? = dataobj(efile?);', ic);
    c_eval('Ne = get_ts(eobj?, ''mms?_des_numberdensity_fast'');', ic);
    c_eval('dbcsVe = get_ts(eobj?, ''mms?_des_bulkv_dbcs_fast'');', ic);
    c_eval('paraTe = get_ts(eobj?, ''mms?_des_temppara_fast'');', ic);
    c_eval('perpTe = get_ts(eobj?, ''mms?_des_tempperp_fast'');', ic);
    c_eval('enfluxe = get_ts(eobj?, ''mms?_des_energyspectr_omni_fast'');', ic);  
    c_eval('eenergy = get_ts(eobj?, ''mms?_des_energy_fast'');', ic);
    Tint = irf.tint(Ne.time.start, Ne.time.stop);
    
    % 2.2. load dfg and edp data
    dmpaB = mms.db_get_ts('mms1_dfg_srvy_ql', 'mms1_dfg_srvy_dmpa', Tint);
    Exyz = mms.db_get_ts('mms1_edp_fast_ql_dce', 'mms1_edp_dce_xyz_dsl', Tint);
   
%%  3. make omni flux
    % 3.1. ion omni flux;
    ispec=struct('t', enfluxi.time.epochUnix);
    ispec.p = enfluxi.data;
    ispec.p_label={'flux', enfluxi.units};
    ispec.f_label={'Energy', '[eV]'};
    ispec.f = single(ienergy.data);
    % 3.2. electron omni flux;
    espec=struct('t', enfluxe.time.epochUnix);
    espec.p = enfluxe.data;
    espec.p_label={'flux', enfluxe.units};
    espec.f_label={'Energy', '[eV]'};
    espec.f = single(eenergy.data);

%%  4. plot
    h = irf_plot(7,'newfigure');
    xSize=800; ySize=650;
    set(gcf,'Position',[10 10 xSize ySize]);
    xwidth = 0.82;          ywidth = 0.128;
    c_eval('set(h(?), ''position'', [0.12 0.960-? * ywidth xwidth ywidth]);', 1: 7);
    
    % 4.1. Bxyz
    h(1) = irf_panel('dmpaB');
    irf_plot(h(1), dmpaB, 'LineWidth', 1.5);
    ylabel(h(1),{'B','[nT]'},'Interpreter','tex');
    irf_legend(h(1),{'B_{x} ',' B_{y} ',' B_{z}'},[0.88 0.15])
    irf_legend(h(1),'(a)',[0.99 0.98],'color','k')
    
    % 4.2. Ni & Ne
    h(2) = irf_panel('ni');
    irf_plot(h(2), Ni, 'LineWidth', 1.5);
    hold(h(2), 'on');
    irf_plot(h(2), Ne, 'LineWidth', 1.5);
    hold(h(2), 'off');
    ylabel(h(2),{'N','[cm^{-3}]'},'Interpreter','tex');
    irf_legend(h(2),{'N_{i} ',' N_{e} '},[0.88 0.15])
    irf_legend(h(2),'(b)',[0.99 0.98],'color','k')    
    
    % 4.3. Vixyz
    h(3) = irf_panel('Vigse');
    irf_plot(h(3), dbcsVi, 'LineWidth', 1.5);
    ylabel(h(3),{'V_i','[km/s]'},'Interpreter','tex');
    irf_legend(h(3),{'V_{x} ',' V_{y} ',' V_{z}'},[0.88 0.15])
    irf_legend(h(3),'(c)',[0.99 0.98],'color','k')    
    
    % 4.4. ion energy flux
    h(4) = irf_panel('iEnflux');
    irf_spectrogram(h(4), ispec, 'log');
    hold(h(4), 'on');
    irf_plot(h(4), paraTi, 'LineWidth', 1.5);
    irf_plot(h(4), perpTi, 'LineWidth', 1.5);
    hold(h(4), 'off');
    h(4).YScale = 'log';
    h(4).YTick = 10.^[1 2 3 4];
    irf_legend(h(4),'(d)',[0.99 0.98],'color','k')
    irf_legend(h(4),'T_{i, ||}',[0.75 0.7],'color','k')
    irf_legend(h(4),'T_{i, \perp}',[0.8 0.7],'color','w')
    ylabel(h(4),{'W_i','[eV]'},'Interpreter','tex');
    
    % 4.6. Vexyz
    h(5) = irf_panel('Vegse');
    irf_plot(h(5), dbcsVe, 'LineWidth', 1.5);
    ylabel(h(5),{'V_e','[km/s]'},'Interpreter','tex');
    irf_legend(h(5),{'V_{x} ',' V_{y} ',' V_{z}'},[0.88 0.15])
    irf_legend(h(5),'(e)',[0.99 0.98],'color','k')     
    
    % 4.7. electron energy flux
    h(6) = irf_panel('eEnflux');
    irf_spectrogram(h(6), espec, 'log');
    hold(h(6), 'on');
    irf_plot(h(6), paraTe, 'LineWidth', 1.5);
    irf_plot(h(6), perpTe, 'LineWidth', 1.5);
    hold(h(6), 'off');
    h(6).YScale = 'log';
    h(6).YTick = 10.^[1 2 3 4];
    irf_legend(h(6),'(f)',[0.99 0.98],'color','k')
    irf_legend(h(6),'T_{e, ||}',[0.75 0.7],'color','k')
    irf_legend(h(6),'T_{e, \perp}',[0.8 0.7],'color','w')
    ylabel(h(6),{'W_e','[eV]'},'Interpreter','tex');
    
    % 4.9. Exyz
    h(7) = irf_panel('Egse');
    irf_plot(h(7), Exyz, 'LineWidth', 1.5);
    ylabel(h(7),{'E','[mV/m]'},'Interpreter','tex');
    irf_legend(h(7),{'E_{x} ',' E_{y} ',' E_{z}'},[0.88 0.15])
    irf_legend(h(7),'(g)',[0.99 0.98],'color','k')    
    
    % 4.X. global
    load('caa/cmap.mat');
    colormap(h(4),cmap);
    colormap(h(6),cmap); 
    title(h(1), strcat('MMS-', num2str(ic)));
    irf_plot_axis_align(h(1 : 7))
    irf_zoom(h(1 : 7),'x',Tint);

%%