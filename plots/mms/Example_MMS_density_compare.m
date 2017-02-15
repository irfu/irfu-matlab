
%%  1. basic
    ic = 1;
    time = irf_time('2016-12-09T17:55:40.000Z','utc>epochtt');
  
%%  2. load data
    % 2.1. PDist
    c_eval('ifile = mms.get_filepath([''mms?_fpi_brst_l2_dis-dist''], time);', ic);
    [iPDist, iPDistError] = mms.make_pdist(ifile);
    c_eval('efile = mms.get_filepath([''mms?_fpi_brst_l2_des-dist''], time);', ic);
    [ePDist, ePDistError] = mms.make_pdist(efile);
    Tint = irf.tint(iPDist.time.start, iPDist.time.stop);
    % 2.2. moment
    c_eval('brstNi? = mms.get_data(''Ni_fpi_brst_l2'', Tint, ?);', 1: 4)
    c_eval('fastNi? = mms.get_data(''Ni_fpi_fast_l2'', Tint, ?);', 1: 4)
    c_eval('brstNe? = mms.get_data(''Ne_fpi_brst_l2'', Tint, ?);', 1: 4)
    c_eval('fastNe? = mms.get_data(''Ne_fpi_fast_l2'', Tint, ?);', 1: 4)   
    % 2.3. get B & E field data
    c_eval('gseB=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',Tint);',ic);  
    c_eval('scPot=mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_scpot_fast_l2'',Tint);',ic);
    % 2.4. get epsd data
    c_eval('epsd_file? = mms.get_filepath([''mms?_dsp_fast_l2_epsd''], time);', 1: 4);
    c_eval('epsd_obj? = dataobj(epsd_file?);', 1: 4);
    c_eval('epsd_omni? = get_ts(epsd_obj?, ''mms?_dsp_epsd_omni'');', 1: 4);
    c_eval('epsd_freq? = get_variable(epsd_obj?, ''mms?_e_freq'');', 1: 4);
    
%%  3. make data
    % 3.1. make epsd_spectrum data
    c_eval('epsd_spec? = struct(''t'', epsd_omni?.time.epochUnix);', 1: 4);
    c_eval('epsd_spec?.p = epsd_omni?.data;', 1: 4);
    c_eval('epsd_spec?.p_label={''Epsd'', epsd_omni?.units};', 1: 4);
    c_eval('epsd_spec?.f_label={''freq'', ''[Hz]''};', 1: 4);
    c_eval('epsd_spec?.f = single((epsd_freq?.data/9e3).^2);', 1: 4);
    
%%  4. plot
    h = irf_plot(9, 'newfigure');
    xSize=800; ySize=650;
    set(gcf,'Position',[10 10 xSize ySize]);
    xwidth = 0.82;          ywidth = 0.100;
    c_eval('set(h(?), ''position'', [0.12 0.960-? * ywidth xwidth ywidth]);', 1: 9);
    
    % 4.1. Bxyz
    h(1) = irf_panel('Bgse');
    irf_plot(h(1), gseB, 'LineWidth', 1.5);
    hold(h(1), 'on');
    tmpcl = irf_plot(h(1), gseB.abs, 'LineWidth', 1.5);   
    hold(h(1), 'off');
    ylabel(h(1),{'B','[nT]'},'Interpreter','tex');
    irf_legend(h(1),{'B_{x} ',' B_{y} ',' B_{z} ', ' |B|'},[0.88 0.15])
    irf_legend(h(1),'(a)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    
    % 4.2. Ni & Ne
    h(2) = irf_panel('N');
    c_eval('irf_plot(h(2), brstNi?, ''LineWidth'', 1.5);', ic);
    hold(h(2), 'on');
    c_eval('irf_plot(h(2), brstNe?, ''LineWidth'', 1.5, ''color'', ''b'');', ic);
    c_eval('irf_plot(h(2), fastNi?, ''LineWidth'', 1.5, ''color'', ''r'');', ic);
    c_eval('irf_plot(h(2), fastNe?, ''LineWidth'', 1.5, ''color'', [0 0.5 0]);', ic);
    hold(h(2), 'off');
    %h(2).YScale = 'log';
    ylabel(h(2),{'N','[cm^{-3}]'},'Interpreter','tex');
    irf_legend(h(2),{'N_{i}^{brst} '},[1.02 1.0], 'color', 'k')
    irf_legend(h(2),{'N_{e}^{brst} '},[1.02 0.73], 'color', 'b')
    irf_legend(h(2),{'N_{i}^{fast} '},[1.02 0.27], 'color', 'r')
    irf_legend(h(2),{'N_{e}^{fast} '},[1.02 0.0], 'color', [0 0.5 0])    
    irf_legend(h(2),'(b)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    %irf_zoom(h(2), 'y', [0.1, 150])
    %h(2).YTick = 10.^[0 1 2];
    
    % 4.3. scPot
    h(3) = irf_panel('scpot');
    irf_plot(h(3), scPot, 'LineWidth', 1.5);
    ylabel(h(3),{'scPot','[V]'},'Interpreter','tex');
    irf_legend(h(3),'(c)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');     
    
    % 4.4. ion energy flux
    h(4) = irf_panel('iEnflux');
    irf_spectrogram(h(4), iPDist.deflux.omni.specrec, 'log');
    hold(h(4), 'on');
    irf_plot(h(4), scPot, 'LineWidth', 1.5, 'color', 'k');
    hold(h(4), 'off');
    h(4).YScale = 'log';
    h(4).YTick = 10.^[1 2 3 4];
    irf_legend(h(4),'(d)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    ylabel(h(4),{'W_i','[eV]'},'Interpreter','tex');

    % 4.5. electron energy flux
    h(5) = irf_panel('eEnflux');
    irf_spectrogram(h(5), ePDist.deflux.omni.specrec, 'log');
    hold(h(5), 'on');
    irf_plot(h(5), scPot, 'LineWidth', 1.5, 'color', 'k');
    hold(h(5), 'off');
    h(5).YScale = 'log';
    h(5).YTick = 10.^[1 2 3 4];
    irf_legend(h(5),'(e)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    ylabel(h(5),{'W_e','[eV]'},'Interpreter','tex');    
      
    % 4.6. epsd1
    h(6) = irf_panel('epsd1');
    irf_spectrogram(h(6), epsd_spec1, 'log');
    hold(h(6), 'on');
    irf_plot(h(6), brstNi1, 'LineWidth', 0.5, 'color', 'k');
    irf_plot(h(6), brstNe1, 'LineWidth', 0.5, 'color', 'b');    
    irf_plot(h(6), fastNi1, 'LineWidth', 0.5, 'color', 'r');
    irf_plot(h(6), fastNe1, 'LineWidth', 0.5, 'color', [0 0.5 0]);
    hold(h(6), 'off');
    %h(6).YScale = 'log';
    irf_legend(h(6),'(f)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    %caxis(h(6),[-9 -6]);
    ylabel(h(6),{'N_1','[cm^{-3}]'},'Interpreter','tex');    
    irf_zoom(h(6), 'y', [1 50])
    
    % 4.7. epsd2
    h(7) = irf_panel('epsd2');
    irf_spectrogram(h(7), epsd_spec2, 'log');
    hold(h(7), 'on');
    irf_plot(h(7), brstNi2, 'LineWidth', 0.5, 'color', 'k');
    irf_plot(h(7), brstNe2, 'LineWidth', 0.5, 'color', 'b');    
    irf_plot(h(7), fastNi2, 'LineWidth', 0.5, 'color', 'r');
    irf_plot(h(7), fastNe2, 'LineWidth', 0.5, 'color', [0 0.5 0]);
    hold(h(7), 'off');
    %h(6).YScale = 'log';
    irf_legend(h(7),'(g)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    %caxis(h(6),[-9 -6]);
    ylabel(h(7),{'N_2','[cm^{-3}]'},'Interpreter','tex');    
    irf_zoom(h(7), 'y', [1 50])
    
    % 4.8. epsd3
    h(8) = irf_panel('epsd3');
    irf_spectrogram(h(8), epsd_spec3, 'log');
    hold(h(8), 'on');
    irf_plot(h(8), brstNi3, 'LineWidth', 0.5, 'color', 'k');
    irf_plot(h(8), brstNe3, 'LineWidth', 0.5, 'color', 'b');    
    irf_plot(h(8), fastNi3, 'LineWidth', 0.5, 'color', 'r');
    irf_plot(h(8), fastNe3, 'LineWidth', 0.5, 'color', [0 0.5 0]);
    hold(h(8), 'off');
    %h(6).YScale = 'log';
    irf_legend(h(8),'(h)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    %caxis(h(6),[-9 -6]);
    ylabel(h(8),{'N_3','[cm^{-3}]'},'Interpreter','tex');    
    irf_zoom(h(8), 'y', [1 50])
 
    % 4.9. epsd4
    h(9) = irf_panel('epsd4');
    irf_spectrogram(h(9), epsd_spec4, 'log');
    hold(h(9), 'on');
    irf_plot(h(9), brstNi4, 'LineWidth', 0.5, 'color', 'k');
    irf_plot(h(9), brstNe4, 'LineWidth', 0.5, 'color', 'b');    
    irf_plot(h(9), fastNi4, 'LineWidth', 0.5, 'color', 'r');
    irf_plot(h(9), fastNe4, 'LineWidth', 0.5, 'color', [0 0.5 0]);
    hold(h(9), 'off');
    %h(6).YScale = 'log';
    irf_legend(h(9),'(i)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    %caxis(h(6),[-9 -6]);
    ylabel(h(9),{'N_4','[cm^{-3}]'},'Interpreter','tex');    
    irf_zoom(h(9), 'y', [1 50])    
    
    % 4.X. global
    load('caa/cmap.mat');
    c_eval('colormap(h(?),cmap);', 4: 9)
    title(h(1), strcat('MMS-', num2str(ic)));
    irf_plot_axis_align(h(1 : 9))
    irf_zoom(h(1 : 9),'x',Tint);

%%