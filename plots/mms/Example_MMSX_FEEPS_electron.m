
%%   Example to plot FEEPS/e data received from FEEPS team (Drew Turner, drew.lawson.turner@gmail.com). 
%   - Data:
%       1. OMNI particle flux of FEEPS-electron: "*_flux_omni_avg.txt";
%       2. energy table [15; keV]: "*_flux_omni_avg_v.txt";
%       3. pitch angle distribution of one energy level (e1 - e15): "*_flux_pad_e?_smth.txt";
%       4. pitch angle [19; -5 - 185 deg]: "*_flux_pad_e?_smth_v.txt";
%   - History:
%       1. wyli, 2016-03-11, irfu.
%       2. wyli, 2016-10-27, irfu; update code, IDL/tplot2txt.pro 

%%  1. basic
    ic = 1;
    Tint = irf.tint('2015-11-01T03:37:00.00Z/2015-11-01T03:41:10.00Z'); 
    feeps_dr = '/Users/wyli/wyli/work/data/feeps/20151101/feeps_ascii/';
    mmsx_omni_avg_fn = 'mmsx_feeps_brst_flux_omni_avg.txt';
    mmsx_omni_energy_fn = 'mmsx_feeps_brst_flux_omni_avg_v.txt';
    mmsx_omni_fn = 'mmsx_feeps_brst_flux_omni.txt';    
    c_eval('mmsx_pad_e?_fn = ''mmsx_feeps_brst_flux_pad_e?_smth.txt'';', 1: 15);
    c_eval('mmsx_pad_e?_PA_fn = ''mmsx_feeps_brst_flux_pad_e?_smth_v.txt'';', 1: 15);
    
%%  2. load MMSX-FEEPS data
    % 2.1. omni        
    formatomni=['%f-%f-%f/%f:%f:%f', repmat('%f', [1, 15])];
    IDomni = fopen([feeps_dr mmsx_omni_fn], 'r');
    IDomni_avg = fopen([feeps_dr mmsx_omni_avg_fn], 'r');    
    tmpData_omni = textscan(IDomni, formatomni);
    tmpData_omni_avg = textscan(IDomni_avg, formatomni);
    fclose(IDomni);         fclose(IDomni_avg);
    % 2.2. omni energy
    energy = load([feeps_dr mmsx_omni_energy_fn]);
    % 2.3. PAD E1-15
    formatpad=['%f-%f-%f/%f:%f:%f', repmat('%f', [1, 19])];
    c_eval('IDpad_e? = fopen([feeps_dr mmsx_pad_e?_fn], ''r'');', 1: 15);
    c_eval('tmpData_pad_e? = textscan(IDpad_e?, formatpad);', 1: 15);
    c_eval('fclose(IDpad_e?);', 1: 15);
    % 2.4. PAD E1-15 pitch angle
    c_eval('pa_e? = load([feeps_dr mmsx_pad_e?_PA_fn]);', 1: 15);

%%  3. load Bxyz
    c_eval('gseB=mms.get_data(''B_gse_fgm_srvy_l2'', Tint, ?);',ic);
    
%%  4. make
    % 4.1. omni & omni_avg data
    tomni = irf_time([tmpData_omni{1, 1}, tmpData_omni{1, 2}, tmpData_omni{1, 3}, ...
        tmpData_omni{1, 4}, tmpData_omni{1, 5}, tmpData_omni{1, 6}], 'vector6>epochtt');
    Data_omni = [tmpData_omni{1, 7: 21}];
    tomni_avg = irf_time([tmpData_omni_avg{1, 1}, tmpData_omni_avg{1, 2}, tmpData_omni_avg{1, 3}, ...
        tmpData_omni_avg{1, 4}, tmpData_omni_avg{1, 5}, tmpData_omni_avg{1, 6}], 'vector6>epochtt');
    Data_omni_avg = [tmpData_omni_avg{1, 7: 21}];    
    % 4.2. pad data
    c_eval('tpad_e? = irf_time([tmpData_pad_e?{1, 1}, tmpData_pad_e?{1, 2}, tmpData_pad_e?{1, 3}, tmpData_pad_e?{1, 4}, tmpData_pad_e?{1, 5}, tmpData_pad_e?{1, 6}], ''vector6>epochtt'');', 1: 15);
    c_eval('Data_pad_e? = [tmpData_pad_e?{1, 7: 25}];', 1: 15);
    % 4.3. make omni spectrum structure;
    Data_omni((Data_omni < 0.)) = NaN;
    specomni =struct('t', tomni.epochUnix);
    specomni.p = double(Data_omni); 
    specomni.p_label={'dF','#/(cm^2 s sr keV)'};
    specomni.f_label={''};
    specomni.f = single(energy);
    % 4.4. make omni_avg spectrum structure;  
    Data_omni_avg((Data_omni_avg < 0.)) = NaN;    
    specomni_avg =struct('t', tomni_avg.epochUnix);
    specomni_avg.p = double(Data_omni_avg); 
    specomni_avg.p_label={'dF','#/(cm^2 s sr keV)'};
    specomni_avg.f_label={''};
    specomni_avg.f = single(energy);
    % 4.5. make pitch angle spectrum structure;  
    c_eval('Data_pad_e?(Data_pad_e?<0) = NaN;', 1: 15);
    c_eval('specPAD_e? = struct(''t'', tpad_e?.epochUnix);', 1: 15);
    c_eval('specPAD_e?.p = double(Data_pad_e?);', 1: 15);
    c_eval('specPAD_e?.p_label={''dF'',''#/(cm^2 s sr keV)''};', 1: 15);
    c_eval('specPAD_e?.f_label={''''};', 1: 15);
    c_eval('specPAD_e?.f = single(pa_e?);', 1: 15);
    % 4.6. 52 keV flux
    dflux52 = irf.ts_scalar(tomni, Data_omni(:, 2));
   
%%  5. plot
    h = irf_plot(6,'newfigure');
    xSize=900; ySize=650;
    set(gcf,'Position',[10 10 xSize ySize]);
    xwidth = 0.82;          ywidth = 0.15;
    c_eval('set(h(?), ''position'', [0.12 0.965-? * ywidth xwidth ywidth]);', 1: 6);
    
    % 5.1. Bxyz
    h(1) = irf_panel('Bgse');
    irf_plot(h(1), gseB, 'LineWidth', 1.5);
    ylabel(h(1),{'B','[nT]'},'Interpreter','tex');
    irf_legend(h(1),{'B_{x} ',' B_{y} ',' B_{z}'},[0.88 0.15])
    irf_legend(h(1),'(a)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    
    % 5.2. omni;
    h(2) = irf_panel('omni');
    irf_spectrogram(h(2), specomni, 'log');
    caxis(h(2),[1 5]);
    ylabel(h(2), {'W_{e}', '[keV]'},'Interpreter','tex');
    irf_legend(h(2),'(b)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    h(2).YScale = 'log';
    irf_zoom(h(2), 'y', [45 400]);    
    
    % 5.3. omni;
    h(3) = irf_panel('omni_avg');
    irf_spectrogram(h(3), specomni_avg, 'log');
    caxis(h(3), [0 1.8]);
    ylabel(h(3), {'W_{e}^{avg}', '[keV]'},'Interpreter','tex');
    irf_legend(h(3),'(c)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    h(3).YScale = 'log';
    irf_zoom(h(3), 'y', [45 400]);
    
    % 5.4. PAD - e2;
    h(4) = irf_panel('PAD_e2');
    irf_spectrogram(h(4), specPAD_e2, 'log');
    caxis(h(4),[0 2]);
    ylabel(h(4), {'52 keV', '[deg]'},'Interpreter','tex');
    irf_legend(h(4),'(d)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');    
    h(4).YTick = [0 45 90 135 180];
    
    % 5.5. PAD - e3;
    h(5) = irf_panel('PAD_e3');    
    irf_spectrogram(h(5), specPAD_e3, 'log');
    caxis(h(5),[0 1.5]);
    ylabel(h(5), {'70 keV', '[deg]'},'Interpreter','tex');
    irf_legend(h(5),'(e)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    h(5).YTick = [0 45 90 135 180];

    % 5.6. PAD - e5;
    h(6) = irf_panel('PAD_e5');    
    irf_spectrogram(h(6), specPAD_e5, 'log');
    caxis(h(6),[0 1.2]);
    ylabel(h(6), {'107 keV', '[deg]'},'Interpreter','tex');
    irf_legend(h(6),'(f)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
    h(6).YTick = [0 45 90 135 180];    
    
    % 5.X. global set
    load('caa/cmap.mat');
    colormap(h(2),cmap);
    colormap(h(3),cmap); 
    colormap(h(4),cmap); 
    colormap(h(5),cmap); 
    colormap(h(6),cmap);     
    irf_plot_axis_align(h(1: 6));         % align the width of all panels
    irf_zoom(h(1: 6),'x',Tint);       % zoom all panels to the same time interval
    tinthl1 = irf_time('2015-11-01T03:37:30.0000Z', 'utc>epoch');    
    tinthl2 = irf_time('2015-11-01T03:40:50.8000Z', 'utc>epoch');
    irf_pl_mark(h(1:6), tinthl1, 'red', 'LineStyle', '--', 'LineWidth', 2.5);    
    irf_pl_mark(h(1:6), tinthl2, 'blue', 'LineStyle', '--', 'LineWidth', 2.5);
    title(h(1), 'MMSX FEEPS/e');
%%

% IDL/tplot2txt.pro: convert FEEPS-tplot data to *.txt data; need SPEDAS toolbox
%   PRO FEEPS_20151101
%   filenames = 'C:\Users\wyli\Desktop\20151101\mmsx_data_4Wenya_01Nov2015.tplot'
%   cd, 'C:\Users\wyli\Desktop\20151101\feeps_ascii'        ; windows system
%   tplot_restore ,filenames=filenames, all=all, sort=sort, get_tvars=get_tvars
%   tplot_ascii, 'mmsx_feeps_brst_flux_omni_avg'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_omni'   
%   tplot_ascii, 'mmsx_fgm_Bcomps_gsm'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e1'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e2'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e3'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e4'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e5'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e6'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e7'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e8'   
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e9'              ; ... 10 - 15  
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e1_smth'       ; 1-15
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e2_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e3_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e4_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e5_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e6_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e7_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e8_smth'
%   tplot_ascii, 'mmsx_feeps_brst_flux_pad_e9_smth'        ; ... 10 - 15
%   tplot_ascii, 'mmsx_feeps_brst_elec_pas'
%   end
%%