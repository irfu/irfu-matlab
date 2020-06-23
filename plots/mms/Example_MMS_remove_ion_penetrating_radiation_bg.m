%%  Example script to show how to remove the ion background mostly due to the penetrating radiation
%   

%%  1. basic
    ic = 1;
    ie = 'i';
    Tint = irf.tint('2017-07-06T22:00:00.00Z/2017-07-06T22:50:00.00Z');
    tint_pl = irf.tint('2017-07-06T22:31:00.00Z/2017-07-06T22:37:00.00Z');
    
%%  2. load data    
    % 2.1. load FPI fast mode
    c_eval('Ni = mms.get_data(''Ni_fpi_fast_l2'',Tint, ?);', ic);
    c_eval('Ne = mms.get_data(''Ne_fpi_fast_l2'',Tint, ?);', ic);
    c_eval('gseVi = mms.get_data(''Vi_gse_fpi_fast_l2'',Tint, ?);', ic);
    c_eval('gseVe = mms.get_data(''Ve_gse_fpi_fast_l2'',Tint, ?);', ic);
    c_eval('gsePi = mms.get_data([''Pi_gse_fpi_fast_l2''], Tint, ?);', ic);
    c_eval('idist = mms.get_data(''PDi_fpi_fast_l2'',Tint, ?);',ic);        
    c_eval('edist = mms.get_data(''PDe_fpi_fast_l2'',Tint, ?);',ic);     
    % 2.2. load Bgse data
    c_eval('gseB=mms.get_data(''B_gse_fgm_srvy_l2'', Tint, ?);',ic);        
    % 2.3. load EDP data
    c_eval('gseE=mms.db_get_ts(''mms?_edp_fast_l2_dce'',''mms?_edp_dce_gse_fast_l2'',Tint);',ic);
    c_eval('scPot=mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_scpot_fast_l2'',Tint);',ic);
    % 2.4. get epsd data
    c_eval('epsd_file = mms.get_filepath([''mms?_dsp_fast_l2_epsd''], Tint);', ic);
    c_eval('epsd_obj = dataobj(epsd_file);', ic);
    c_eval('epsd_omni = get_ts(epsd_obj, ''mms?_dsp_epsd_omni'');', ic);
    c_eval('epsd_freq = get_variable(epsd_obj, ''mms?_e_freq'');', ic);  
    
%%  3. compute
    % 3.1. remove idist background
    [enflux_new, enflux_BG, idist_new, idist_BG, Ni_new, gseVi_new, gsePi_new, ...
        Ni_bg, EnergySpectr_bg, Pres_bg, EnergySpectr_bg_self]= mms.remove_ion_penetrating_radiation_bg(idist, 'tint', Tint);

    % 3.2. make epsd structure;
    %c_eval('epsd_omni?.data(epsd_omni?.data < 1e-14) = NaN;', 1: 4);
    tmp = epsd_omni.data;
    epsd_spec = struct('t', epsd_omni.time.epochUnix);
    epsd_spec.p = epsd_omni.data;
    epsd_spec.p_label = {'Epsd', epsd_omni.units};
    epsd_spec.f_label = {'freq', '[Hz]'};
    epsd_spec.f = single((epsd_freq.data/9e3).^2);
    
%%  4. plot
    npanel1 = 6;
    h=irf_plot(npanel1,'reset');
    xSize = 900; ySize = 700;
    set(gcf,'Position',[10 10 xSize ySize]);
    xwidth = 0.85;          ywidth = 0.90/npanel1;
    c_eval('set(h(?), ''position'', [0.11 0.952-? * ywidth xwidth ywidth]);', 1: npanel1);
   
    % 4.1. MMS1 DFG
    ip = 1;
    h(ip) = irf_panel('gseB');
    irf_plot(h(ip), gseB, 'LineWidth', 1.8);    
    ylabel(h(ip), {' B ','[nT]'}, 'Interpreter', 'tex');
    irf_legend(h(ip), {'  B_{X} ', ' B_{Y} ', ' B_{Z}  '},[0.92 0.82])
    irf_legend(h(ip),'(a)',[0.99 0.95],'color','k','fontsize',12, 'fontWeight', 'bold');
    irf_zoom(h(ip), 'y', [-8 25]);
    h(ip).FontSize = 14;
    
    % 4.2. MMS1 Ni
%     ip = ip + 1;
%     h(ip) = irf_panel('Ni');
%     irf_plot(h(ip), Ni, 'LineWidth', 1.8);
%     hold(h(ip), 'on');
%     irf_plot(h(ip), Ne, 'LineWidth', 1.8); 
%     irf_plot(h(ip), Ni_bg, 'LineWidth', 1.8); 
%     irf_plot(h(ip), Ni - Ni_bg, 'LineWidth', 1.8);
%     hold(h(ip), 'off');
%     ylabel(h(ip), {' N ','[cm^{-3}]'}, 'Interpreter', 'tex');
%     %irf_legend(h(2), {'   N_i   ', '   N_e   '}, [0.08 0.70]);    
%     irf_legend(h(ip), {' N_i ', ' N_e ', ' N_{i,BG} ', ' N_i - N_{i,BG}'},[0.92 0.20]);
%     irf_legend(h(ip),'(b)',[0.99 0.95],'color','k','fontsize',12, 'fontWeight', 'bold');    
%     irf_zoom(h(ip), 'y', [0 0.4]);  
%     h(ip).FontSize = 14;
   
    % 4.3. MMS Epsd
    ip = ip + 1;
    h(ip) = irf_panel('epsd');
    [hepsd, hcepsd] = irf_spectrogram(h(ip), epsd_spec, 'log');
    hold(h(ip), 'on');
    irf_plot(h(ip), Ni, 'LineWidth', 1.8);
    irf_plot(h(ip), Ne, 'LineWidth', 1.8); 
    irf_plot(h(ip), Ni_bg, 'LineWidth', 1.8); 
    irf_plot(h(ip), Ni - Ni_bg, 'LineWidth', 1.8);    
    hold(h(ip), 'off');
    irf_legend(h(ip), '(b)', [0.99 0.99], 'color','k','fontsize',12, 'fontWeight', 'bold');
    caxis(h(ip), [-9 -4]);    
    ylabel(h(ip), {'N','[cm^{-3}]'},'Interpreter','tex');    
    irf_zoom(h(ip), 'y', [0 0.2]);  
    irf_legend(h(ip), {' N_i ', ' N_e ', ' N_{i,BG} ', ' N_i - N_{i,BG}'},[0.92 0.20]);
    h(ip).FontSize = 14;    
    
    % 4.5. MMS1 dEFlux ion
    ip = ip + 1;
    h(ip) = irf_panel('ienflux');
    irf_spectrogram(h(ip), idist.deflux.omni.specrec,'log');
    set(h(ip), 'yscale','log');
    set(h(ip), 'ytick',[1e1 1e2 1e3 1e4]);
    ylabel(h(ip), {'W_{i}', '[eV]'},'Interpreter','tex');
    irf_legend(h(ip),'(c)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    %irf_legend(h(ip),'W_{ExB}',[0.18 0.18],'color','k','fontsize',16);
    h(ip).FontSize = 14;
    caxis(h(ip), [1.5 6.5]);    
    
    % 4.2. MMS1 background enflux
    ip = ip + 1;
    h(ip) = irf_panel('ienflux bg profiles');
    irf_plot(h(ip), EnergySpectr_bg, 'LineWidth', 1.8, 'color', 'r');
    hold(h(ip), 'on');
    irf_plot(h(ip), EnergySpectr_bg_self, 'LineWidth', 1.8, 'color', 'k'); 
    hold(h(ip), 'off');
    ylabel(h(ip), {' Enflux ','[keV/(cm^2 s sr keV)]'}, 'Interpreter', 'tex');
    %irf_legend(h(2), {'   N_i   ', '   N_e   '}, [0.08 0.70]);    
    irf_legend(h(ip), {'Self-Computed BG ', ' ', ' Official BG '}, [0.10 0.20]);
    irf_legend(h(ip),'(d)',[0.99 0.95],'color','k','fontsize',12, 'fontWeight', 'bold');    
    h(ip).FontSize = 14;    
    
    % 4.6. MMS1 dEFlux ion BG
    ip = ip + 1;
    h(ip) = irf_panel('ienflux bg');
    irf_spectrogram(h(ip), enflux_BG,'log');
    set(h(ip), 'yscale','log');
    set(h(ip), 'ytick',[1e1 1e2 1e3 1e4]);
    ylabel(h(ip), {'W_{i}', '[eV]'},'Interpreter','tex');
    irf_legend(h(ip),'(e)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    %irf_legend(h(ip),'W_{ExB}',[0.18 0.18],'color','k','fontsize',16);
    h(ip).FontSize = 14;
    caxis(h(ip), [1.5 6.5]);    
    
    % 4.6. MMS1 dEFlux ion NEW    
    ip = ip + 1;
    h(ip) = irf_panel('ienflux new');
    irf_spectrogram(h(ip), enflux_new, 'log');
    set(h(ip), 'yscale','log');
    set(h(ip), 'ytick',[1e1 1e2 1e3 1e4]);
    ylabel(h(ip), {'W_{i}', '[eV]'},'Interpreter','tex');
    irf_legend(h(ip),'(f)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    %irf_legend(h(ip),'W_{ExB}',[0.18 0.18],'color','k','fontsize',16);
    h(ip).FontSize = 14;  
    caxis(h(ip), [1.5 6.5]);
  
    % 4.X. Mark & time
    load('caa/cmap.mat');
    c_eval('colormap(h(?),cmap);', [2, 3, 5, 6]);
    irf_plot_axis_align(h(1 : npanel1));         % align the width of all panels
    irf_zoom(h(1 : npanel1),'x', tint_pl);       % zoom all panels to the same time interval
    title(h(1), 'MMS1 2017-07-06 UTC');
    h(npanel1).XLabel.String = '';

%%