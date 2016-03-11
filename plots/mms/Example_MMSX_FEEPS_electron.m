
%%   Example to plot FEEPS data received from FEEPS team (Drew Turner). 
%   - Data:
%       1. OMNI particle flux of FEEPS-electron: "*_flux_omni_avg.txt";
%       2. energy table [15; keV]: "*_flux_omni_avg_v.txt";
%       3. pitch angle distribution of one energy level (e3: 70 keV; e5: 107 keV, e9: 200 keV): "*_flux_pad_e?_smth.txt";
%       4. pitch angle [19; -5 - 185 deg]: "*_flux_pad_e?_smth_v.txt";
%   - History:
%       1. wyli, 2016-03-11, irfu.

%%  1. basic
    %Tint = irf.tint('2015-09-08T10:29:26.00Z/2015-09-08T10:29:33.00Z');        
    %Tint = irf.tint('2015-09-08T10:44:00.00Z/2015-09-08T10:45:00.00Z');    
    Tint = irf.tint('2015-09-08T10:44:30.00Z/2015-09-08T10:44:40.00Z'); 
    %Tint = irf.tint('2015-12-06T23:38:26.00Z/2015-12-06T23:38:36.00Z');            
    feeps_dr = '.../feeps/...';
    mmsx_omni_fn = 'mmsx_feeps_brst_flux_omni_avg.txt';
    mmsx_energy_fn = 'mmsx_feeps_brst_flux_omni_avg_v.txt';
    mmsx_pad_e3_fn = 'mmsx_feeps_brst_flux_pad_e3_smth.txt';
    mmsx_pad_e3_PA_fn = 'mmsx_feeps_brst_flux_pad_e3_smth_v.txt';
    mmsx_pad_e5_fn = 'mmsx_feeps_brst_flux_pad_e5_smth.txt';
    mmsx_pad_e5_PA_fn = 'mmsx_feeps_brst_flux_pad_e5_smth_v.txt';    
    mmsx_pad_e9_fn = 'mmsx_feeps_brst_flux_pad_e9_smth.txt';
    mmsx_pad_e9_PA_fn = 'mmsx_feeps_brst_flux_pad_e9_smth_v.txt';

%%  2. load
    % 2.1. omni        
    formatomni=['%f-%f-%f/%f:%f:%f', repmat('%f', [1, 15]),  '\n'];
    IDomni = fopen([feeps_dr mmsx_omni_fn], 'r');
    tmpData_omni = textscan(IDomni, formatomni);
    fclose(IDomni);
    % 2.2. omni energy
    energy = load([feeps_dr mmsx_energy_fn]);
    % 2.3. PAD E3/5/9
    formatpad=['%f-%f-%f/%f:%f:%f', repmat('%f', [1, 19]),  '\n'];
    c_eval('IDpad_e? = fopen([feeps_dr mmsx_pad_e?_fn], ''r'');', [3, 5, 9]);
    c_eval('tmpData_pad_e? = textscan(IDpad_e?, formatpad);', [3, 5, 9]);
    c_eval('fclose(IDpad_e?);', [3, 5, 9]);
    % 2.4. PAD E3/5/9 pitch angle
    c_eval('pa_e? = load([feeps_dr mmsx_pad_e?_PA_fn]);', [3, 5, 9]);
    
%%  3. make
    Data_pad_e3(Data_pad_e3 < 0.1) = NaN; 
    % 3.1. omni
    tomni = irf_time([tmpData_omni{1, 1}, tmpData_omni{1, 2}, tmpData_omni{1, 3}, ...
        tmpData_omni{1, 4}, tmpData_omni{1, 5}, tmpData_omni{1, 6}], 'vector6>epochtt');
    Data_omni = [tmpData_omni{1, 7: 21}];
    % 3.2. pad
    c_eval('tpad_e? = irf_time([tmpData_pad_e?{1, 1}, tmpData_pad_e?{1, 2}, tmpData_pad_e?{1, 3}, tmpData_pad_e?{1, 4}, tmpData_pad_e?{1, 5}, tmpData_pad_e?{1, 6}], ''vector6>epochtt'');', [3, 5, 9]);
    c_eval('Data_pad_e? = [tmpData_pad_e?{1, 7: 25}];', [3, 5, 9]);
    % 3.3. make omni spectrum structure;    
    specomni =struct('t', tomni.epochUnix);
    specomni.p = double(Data_omni); 
    specomni.p_label={'dF','#/(cm^2 s sr keV)'};
    specomni.f_label={''};
    specomni.f = single(energy);

    % 3.4. make pitch angle spectrum structure;    
    c_eval('specPAD_e? =struct(''t'', tpad_e?.epochUnix);', [3, 5, 9]);
    c_eval('specPAD_e?.p = double(Data_pad_e?);', [3, 5, 9]);
    c_eval('specPAD_e?.p_label={''dF'',''#/(cm^2 s sr keV)''};', [3, 5, 9]);
    c_eval('specPAD_e?.f_label={''''};', [3, 5, 9]);
    c_eval('specPAD_e?.f = single(pa_e?);', [3, 5, 9]);    
    
%%  4. plot
    h = irf_plot(4, 'newfigure');
    xSize=900; ySize=650;
    set(gcf,'Position',[10 10 xSize ySize]); %position and size of the figure in the screen
% 4.1. omni;
    irf_spectrogram(h(1), specomni, 'log');
    %irf_legend((8),'(d)',[0.99 0.98],'color','w','fontsize',12);
    %caxis(h(1),[1 2000]);
    %irf_legend(hca,{'V'},[0.5 0.6]);
    %irf_legend(hca,'V_{ExB}',[0.45 0.6],'color','r');
    ylabel(h(1), {'E_{i}', '[eV]'},'Interpreter','tex');
    irf_legend(h(1),'(h)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
% 4.2. PAD - e3;
    irf_spectrogram(h(2), specPAD_e3);
    %irf_legend((8),'(d)',[0.99 0.98],'color','w','fontsize',12);
    caxis(h(2),[0 10]);
    %irf_legend(hca,{'V'},[0.5 0.6]);
    %irf_legend(hca,'V_{ExB}',[0.45 0.6],'color','r');
    ylabel(h(2), {'70 keV PAD', '[deg]'},'Interpreter','tex');
    %irf_legend(h(2),'(h)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
% 4.3. PAD - e5;
    irf_spectrogram(h(3), specPAD_e5);
    %irf_legend((8),'(d)',[0.99 0.98],'color','w','fontsize',12);
    %caxis(hca,[4 8]);
    %irf_legend(hca,{'V'},[0.5 0.6]);
    caxis(h(3),[0 4]);
    
    %irf_legend(hca,'V_{ExB}',[0.45 0.6],'color','r');
    ylabel(h(3), {'107 keV PAD', '[deg]'},'Interpreter','tex');
    %irf_legend(h(2),'(h)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');
% 4.4. PAD - e9;
    irf_spectrogram(h(4), specPAD_e9);
    %caxis(hca,[4 8]);
    caxis(h(4),[0 1]);
    
    %irf_legend(hca,{'V'},[0.5 0.6]);
    %irf_legend(hca,'V_{ExB}',[0.45 0.6],'color','r');
    ylabel(h(4), {'200 keV PAD', '[deg]'},'Interpreter','tex');
    %irf_legend(h(2),'(h)',[0.99 0.99],'color','k','fontsize',12, 'fontWeight', 'bold');

    irf_plot_axis_align(h(1: 4));         % align the width of all panels
    irf_zoom(h(1: 4),'x',Tint);       % zoom all panels to the same time interval
    title(h(1), 'MMSX');
    
%%