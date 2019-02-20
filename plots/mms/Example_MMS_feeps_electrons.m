
%   Example code to plot FEEPS-electron burst-l2 data
%       - The current burst-l2 data is not helpful to see the pitch angle distributino. A better way is to ask
%         help from FEEPS team, e.g. Drew Turner (drew.lawson.turner@gmail.com).
%       - Refer to FEEPS and EPD documents for better understanding of the instrument.
%       - Wenya, 2016-03-11, irfu.

% from SPEDAS: mms_feeps_omni.pro - correction factors per SC
eEcorr = [14.0, -1.0, -3.0, -3.0];
iEcorr = [0.0, 0.0, 0.0, 0.0];
eGfact = [1.0, 1.0, 1.0, 1.0];
iGfact = [0.84, 1.0, 1.0, 1.0]; 

ic = 1;
Tint = irf.tint('2017-07-23T16:40:00.000Z/2017-07-23T17:10:00.000Z');
    

%%  2. Old way
%       2.1. load data

    c_eval('feeps_dr = ''/Volumes/WYLIMMS/work/data/mms?/feeps/brst/l2/electron/2015/09/08/'';', ic);
    tstring = '102834';    
    c_eval('feeps_fn = [''mms?_feeps_brst_l2_electron_20150908'' tstring ''_v5.4.1.cdf''];', ic);
    feeps_obj = dataobj([feeps_dr feeps_fn]);            % [ion moments]
    eLow = get_variable(feeps_obj, 'electron_energy_lower_bound');
    eUp = get_variable(feeps_obj, 'electron_energy_upper_bound');
    eval(['energies = (eLow.data + eUp.data)/2. + ' specie(1) 'Ecorr(ic);'])
    c_eval('eTit1 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_1'');', ic);
    c_eval('eTit2 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_2'');', ic);
    c_eval('eTit3 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_3'');', ic);
    c_eval('eTit4 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_4'');', ic);
    c_eval('eTit5 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_5'');', ic);
    c_eval('eTit9 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_9'');', ic);
    c_eval('eTit10 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_10'');', ic);
    c_eval('eTit11 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_11'');', ic);
    c_eval('eTit12 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_top_intensity_sensorid_12'');', ic);
    c_eval('eBit1 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_1'');', ic);
    c_eval('eBit2 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_2'');', ic);
    c_eval('eBit3 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_3'');', ic);
    c_eval('eBit4 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_4'');', ic);
    c_eval('eBit5 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_5'');', ic);
    c_eval('eBit9 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_9'');', ic);
    c_eval('eBit10 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_10'');', ic);
    c_eval('eBit11 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_11'');', ic);
    c_eval('eBit12 = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_bottom_intensity_sensorid_12'');', ic);
    c_eval('ePA = get_ts(feeps_obj, ''mms?_epd_feeps_brst_l2_electron_pitch_angle'');', ic);
%       2.2. tlim data
    sensors = [1:5, 9:12];
    c_eval('eTit? = eTit?.tlim(Tint);', sensors);   
    c_eval('eBit? = eBit?.tlim(Tint);', sensors);
    ePA = ePA.tlim(Tint);
    
    %% Alt
    specie = 'electron'; % alt specie = 'ion';
    mode = 'brst'; % alt mode = 'fast';
    dsetName = sprintf('mms%d_feeps_%s_l2_%s',ic,mode,specie);
    dsetPref = sprintf('mms%d_epd_feeps_%s_l2_%s',ic,mode,specie);
    switch specie
      case 'ion', sensors = 6:8;
      case 'electron'
        switch mode
          case 'brst', sensors = [1:5, 9:12];
          case 'fast', sensors = [3:5, 11:12];
          otherwise, error('invalid mode')
        end
      otherwise, error('invalid specie')
    end
    
    eLow = mms.db_get_variable(dsetName, [specie '_energy_lower_bound'],Tint);
    eUp = mms.db_get_variable(dsetName, [specie '_energy_upper_bound'],Tint);
    eval(['energies = (eLow.data + eUp.data)/2. + ' specie(1) 'Ecorr(ic);'])
    
    nSensors = length(sensors);
    for iSen = 1:nSensors
      sen = sensors(iSen); suf = sprintf('intensity_sensorid_%d',sen);
      sufMask = sprintf('sector_mask_sensorid_%d',sen);
      top = mms.db_get_ts(dsetName,[dsetPref '_top_' suf],Tint);
      mask = mms.db_get_ts(dsetName,[dsetPref '_top_' sufMask],Tint);
      if all(size((mask.data))==size((top.data)))
        top.data(logical(mask.data)) = NaN;
      else
        top.data(logical(repmat(mask.data,1,length(energies)))) = NaN; % obsolete?
      end
      bot = mms.db_get_ts(dsetName,[dsetPref '_bottom_' suf],Tint);
      mask = mms.db_get_ts(dsetName,[dsetPref '_bottom_' sufMask],Tint);
      if all(size((mask.data))==size((bot.data)))
        bot.data(logical(mask.data)) = NaN;
      else
        bot.data(logical(repmat(mask.data,1,length(energies)))) = NaN;
      end
      c_eval([specie(1) 'Tit?=top;'  specie(1) 'Bit?=bot;'],sen)
    end
    
    % cleanup
    if 0 && strcmp(mode,'brst') 
      switch specie
        case 'electron'
          eTit1.data(:,:) = NaN;
          eBit1.data(:,:) = NaN;
          eTit2.data(:,1) = NaN;
          eTit5.data(:,1) = NaN;
        case 'ion'
          iTit6.data(:,1:2) = NaN;
          iBit6.data(:,1:2) = NaN;
      end
    end
    
    %% Omni
nSensors = length(sensors);
eval(['dTmp=' specie(1) 'Tit' num2str(sensors(1)) ';'])
omniD = NaN( [size(dTmp.data) nSensors*2]);
for iSen = 1:nSensors
  c_eval(['omniD(:,:,iSen) = ' specie(1) 'Tit?.data;'...
    'omniD(:,:,nSensors+iSen) = ' specie(1) 'Bit?.data;'],sensors(iSen))
end
eval([specie(1) 'Omni = dTmp; ' specie(1) 'Omni.data =' ...
  'mean(double(omniD),3,''omitnan'')*' specie(1) 'Gfact(ic);'])
    
%%  3. make    
%       3.2. structure for Top
    c_eval(['speTit? = struct(''t'', ' specie(1) 'Tit?.time.epochUnix);'], sensors);
    c_eval(['speTit?.p = double(' specie(1) 'Tit?.data);'], sensors);
    c_eval('speTit?.p_label = {''intensity''};', sensors);
    c_eval('speTit?.f_label = {''Energy''};', sensors);
    c_eval('speTit?.f = double(energies);', sensors);
    
%       3.3. structure for Bottom
    c_eval(['speBit? = struct(''t'', ' specie(1) 'Bit?.time.epochUnix);'], sensors);
    c_eval(['speBit?.p = double(' specie(1) 'Bit?.data);'], sensors);
    c_eval('speBit?.p_label = {''intensity''};', sensors);
    c_eval('speBit?.f_label = {''Energy''};', sensors);
    c_eval('speBit?.f = double(energies);', sensors);
    
    % omni
    eval(['speOmni = struct(''t'', ' specie(1) 'Omni.time.epochUnix);'...
    'speOmni.p = ' specie(1) 'Omni.data;'...
    'speOmni.p_label = {[''log('' ' specie(1) 'Omni.units '')'']};']);
    
    speOmni.f_label = {'Energy'};
    speOmni.f = double(energies);
    
    %% plot
    
    h = irf_plot(nSensors*2+1,'newfigure');
    
    hca=irf_panel('Omni');
    irf_spectrogram(hca, speOmni, 'log', 'donotfitcolorbarlabel');
    caxis(hca,[-0.5, 3])
    ylabel(hca,'Omni [keV]','Interpreter','tex');
    
    for iSen = 1:nSensors
      sen = sensors(iSen); senS = num2str(sen);
      lookS = 'TB';
      for look = lookS
        hca=irf_panel([look 'op' senS]); eval(['dTmp = spe' look 'it' senS ';'])
        irf_spectrogram(hca, dTmp, 'log', 'donotfitcolorbarlabel');
        caxis(hca,[-0.5, 3])
        ylabel(hca,['s' senS look],'Interpreter','tex');
      end
    end

    irf_zoom(h,'y', [20 200])
    irf_zoom(h,'x',Tint)
    irf_plot_axis_align(h)
    irf_colormap('space')    
    tmpts0 = Tint.start.utc;
    tmpts1 = Tint.stop.utc;
    title(h(1),['MMS' int2str(ic) ' ' tmpts0(12:22) '~' tmpts1(12:22) ]);
    
%%  4. plot
%       4.0.
    h = irf_plot(11,'newfigure');
    
    hca=irf_panel('Omni');
    irf_spectrogram(hca, speOmni, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'Omni [keV]','Interpreter','tex');
    
%       4.1. electron Top 5
    hca=irf_panel('Top1');
    irf_spectrogram(hca, speTit1, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H1U [keV]','Interpreter','tex');

%       4.2. electron Top 3
    hca=irf_panel('Top3');
    irf_spectrogram(hca, speTit3, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H2U','Interpreter','tex');
    
%       4.3. electron Top 9
    hca=irf_panel('Top9');
    irf_spectrogram(hca, speTit9, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H5U','Interpreter','tex');

%       4.4. electron Top 11
    hca=irf_panel('Top11');
    irf_spectrogram(hca, speTit11, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H6U','Interpreter','tex'); 
    
%       4.5. electron Top 5
    hca=irf_panel('Top5');
    irf_spectrogram(hca, speTit5, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H3T','Interpreter','tex');
    
%       4.6. electron Bottom 5
    hca=irf_panel('Bottom5');
    irf_spectrogram(hca, speBit5, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H3B','Interpreter','tex');      
    
%       4.7. electron Bottom 2
    hca=irf_panel('Top2');
    irf_spectrogram(hca, speTit2, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H1D','Interpreter','tex');

%       4.8. electron Top 4
    hca=irf_panel('Top4');
    irf_spectrogram(hca, speTit4, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H2D','Interpreter','tex');     

%       4.9. electron Top 10
    hca=irf_panel('Top10');
    irf_spectrogram(hca, speTit10, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H3D','Interpreter','tex');

%       4.10. electron Top 12
    hca=irf_panel('Top12');
    irf_spectrogram(hca, speTit12, 'log', 'donotfitcolorbarlabel');
    %set(hca,'yscale','log');
    %set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    caxis(hca,[-0.5, 3])
    ylabel(hca,'H4D','Interpreter','tex');     
    
    
%       4.X.
    irf_zoom(h,'y', [20 200])
    irf_zoom(h,'x',Tint)
    irf_plot_axis_align(h)
    %set(h ,'position',[0.15 0.2 0.66 0.71]);    
    %irf_colormap('poynting')    
    %add_position(h(end),Rgse)
    %xlabel(h(end),'');
    tmpts0 = Tint.start.utc;
    tmpts1 = Tint.stop.utc;
    title(h(1),['MMS' int2str(ic) ' ' tmpts0(12:22) '~' tmpts1(12:22) ]);