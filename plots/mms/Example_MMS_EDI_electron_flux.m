% Example_MMS_EDI_electron_flux
mms_id = 1;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');

% dpflux - differential particle flux, 1/(cm^2 s sr keV)
% deflux - differential energy flux, keV/(cm^2 s sr keV)
% pflux - particle flux, 1/(cm^2 s sr)
% eflux - energy flux, keV/(cm^2 s sr)

%% Load data
c_eval('ePitch?_pflux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',mms_id)
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',mms_id);
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
%c_eval('[iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0]));',mms_id)
c_eval('ePitch?_deflux_fpi = mms.get_data(''Enflux-pitchangle-e_fpi_brst_l2'',tint,?)',mms_id)

%% Do stuff with data
% Pitchangle distribution
c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,16);',mms_id)
c_eval('ePitch?_30 = ePDist?.pitchangles(dmpaB?,0:30:180);',mms_id)

% Rebin to different energy range
c_eval('ePDist?_rebin = ePDist?.rebin(''sph'',{25:50:1025,[],[]});',mms_id)
c_eval('fpi_edges_test = ePDist?.elim([480]).depend{1}(1,:) + [-ePDist1.elim([480]).ancillary.delta_energy_minus(1,:)  ePDist1.elim([480]).ancillary.delta_energy_plus(1,:)];',mms_id)
c_eval('ePDist?_rebin_test = ePDist?.rebin(''sph'',{fpi_edges_test+[30 -30],[],[]});',mms_id)
c_eval('ePitch?_rebin_test = ePDist?_rebin_test.pitchangles(dmpaB?,0:11.25:180);',mms_id)

% Rebin to EDI energy range
edi_energy_range = ediFlux1.depend{1}(1,:) + [- ediFlux1.ancillary.delta_energy_minus ediFlux1.ancillary.delta_energy_plus];
c_eval('ePDist?_rebin_edi = ePDist?.rebin(''sph'',{edi_energy_range,[],[]});',mms_id);
c_eval('ePitch?_rebin_edi = ePDist?_rebin_edi.pitchangles(dmpaB?,0:11.25:180);',mms_id)

% Fluxes
c_eval('ePDist?_pflux = ePDist?.flux;',mms_id)
c_eval('ePDist?_dpflux = ePDist?.dpflux;',mms_id)
c_eval('ePDist?_deflux = ePDist?.deflux;',mms_id)

c_eval('ePitch?_pflux = ePitch?.flux;',mms_id)
c_eval('ePitch?_dpflux_ = ePitch?.flux(''diff'');',mms_id)
c_eval('ePitch?_dpflux = ePitch?.dpflux;',mms_id)
c_eval('ePitch?_deflux = ePitch?.deflux;',mms_id)

c_eval('ePitch?_pflux_30 = ePitch?_30.flux;',mms_id)
c_eval('ePitch?_dpflux_30 = ePitch?_30.dpflux;',mms_id)
c_eval('ePitch?_deflux_30 = ePitch?_30.deflux;',mms_id) % compare to ePitch?_deflux_fpi

c_eval('ePitch?_pflux_rebin_edi = ePitch?_rebin_edi.flux;',mms_id)
c_eval('ePitch?_dpflux_rebin_edi = ePitch?_rebin_edi.dpflux;',mms_id)
c_eval('ePitch?_deflux_rebin_edi = ePitch?_rebin_edi.deflux;',mms_id) % compare to ePitch?_edi_pflux

% Calculate FPI flux in varying energy range
c_eval('ePDist?_rebin_var = ePDist?.rebin(''sph'',{[450 550],[],[]});',mms_id);
c_eval('ePitch?_rebin_var = ePDist?_rebin_var.pitchangles(dmpaB?,0:11.25:180);',mms_id)
c_eval('eFlux?_rebin_var = ePitch?_rebin_var.flux;',mms_id)

%c_eval('ePDist?_rebin = ePDist?(1:100).rebin(''sph'',{25:50:1000,[],[]});',1);
%c_eval('ePitch?_rebin = ePDist?_rebin.pitchangles(dmpaB?,0:11.25:180);',mms_id)
%c_eval('eFlux?_rebin = ePitch?_rebin.flux;',mms_id)

%% Plot 1, time series, compare, EDI flux to FPI flux, original energy channel 2*65, and rebinned to 2*25
c_eval('fluxEDI = ePitch?_pflux_edi.palim([135 180]);',mms_id)
c_eval('fluxEDI_resample_fpi = fluxEDI.resample(ePDist?);',mms_id)
c_eval('fluxFPI = ePitch?_pflux.palim([135 180]).elim(ePitch?_pflux_edi.ancillary.energy0);',mms_id)
c_eval('fluxFPI_rebin_edi = ePitch?_pflux_rebin_edi.palim([135 180]).elim(ePitch?_pflux_edi.ancillary.energy0);',mms_id)
%fluxFPI_apar = eFlux1_apar.elim(ediFlux2.ancillary.energy0);
%fluxFPI_apar_30 = eFlux1_apar_30.elim(ediFlux2.ancillary.energy0);
%fluxFPI = eFlux1.palim([135 180]).elim([300 600]);
%fluxFPI_rebin = eFlux1_rebin.palim([135 180]).elim(ediFlux2.ancillary.energy0);
%fluxFPI_rebin_var = eFlux1_rebin_var.palim([135 180]);  
h = irf_plot({fluxEDI,fluxEDI_resample_fpi,fluxFPI,fluxFPI_rebin_edi},'comp'); % ,fluxFPI_rebin_var
h(1).Title.String = sprintf('Electron flux, MMS %',mms_id);
irf_legend(h(1),{sprintf('EDI: [%.0f,%.0f] eV',fluxEDI.depend{1}(1)-fluxEDI.ancillary.delta_energy_minus(1),fluxEDI.depend{1}(1)+fluxEDI.ancillary.delta_energy_plus(1));...
                 sprintf('EDI resampled to FPI time: [%.0f,%.0f] eV',fluxEDI_resample_fpi.depend{1}(1)-fluxEDI_resample_fpi.ancillary.delta_energy_minus(1),fluxEDI_resample_fpi.depend{1}(1)+fluxEDI_resample_fpi.ancillary.delta_energy_plus(1));...
                 sprintf('FPI closest energy channel: [%.0f,%.0f] eV',fluxFPI.depend{1}(1)-fluxFPI.ancillary.delta_energy_minus(1),fluxFPI.depend{1}(1)+fluxFPI.ancillary.delta_energy_plus(1));...
                 sprintf('FPI rebinned to EDI energy channel: [%.0f,%.0f] eV',fluxFPI_rebin_edi.depend{1}(1)-fluxFPI_rebin_edi.ancillary.delta_energy_minus(1),fluxFPI_rebin_edi.depend{1}(1)+fluxFPI_rebin_edi.ancillary.delta_energy_plus(1))},...
                 [0.02 0.98])
%h = irf_plot({fluxEDI,fluxFPI,fluxFPI_rebin_500});
linkprop(h,'YLim');
%%
pa_minus = fluxFPI.depend{2}(1,:)-fluxFPI.ancillary.delta_pitchangle_minus(1,:);
pa_plus = fluxFPI.depend{2}(1,:)+fluxFPI.ancillary.delta_pitchangle_plus(1,:);
c_eval('h(?).YLabel.String = {''Flux'',sprintf(''(%s)'',fluxFPI.units),[''\theta = ''  sprintf(''[%g,%g] deg'',pa_minus(?),pa_plus(?))]};',1:4)
for ii = 1:4
  irf_legend(h(ii),{...
    sprintf('EDI: [%g,%g] eV', fluxEDI.depend{1}(1,:)-fluxEDI.ancillary.delta_energy_minus(1),fluxEDI.depend{1}(1,:)+fluxEDI.ancillary.delta_energy_plus(1));...
    sprintf('FPI channel closest to EDI energy: [%g,%g] eV', fluxFPI.depend{1}(1,:)-fluxFPI.ancillary.delta_energy_minus(1),fluxFPI.depend{1}(1,:)+fluxFPI.ancillary.delta_energy_plus(1));...
    sprintf('FPI rebinned to EDI energy channel: [%g,%g] eV', fluxFPI_rebin_edi.depend{1}(1,:)-fluxFPI_rebin_edi.ancillary.delta_energy_minus(1),fluxFPI_rebin_edi.depend{1}(1,:)+fluxFPI_rebin_edi.ancillary.delta_energy_plus(1));...
    %sprintf('FPI rebinned to energy channel: [%g,%g] eV', fluxFPI_rebin_var.depend{1}(1,:)-fluxFPI_rebin_var.ancillary.delta_energy_minus(1),fluxFPI_rebin_var.depend{1}(1,:)+fluxFPI_rebin_var.ancillary.delta_energy_plus(1))
    },[0.02 0.98],'fontsize',10);
end
irf_zoom(h,'x',tint)
%irf_legend(hca,{sprintf('[%.2f,%.2f]',fluxFPI.ancillary.pitchangle_edges(1:2));...
%                sprintf('[%.2f,%.2f]',fluxFPI.ancillary.pitchangle_edges(2:3));...
%                sprintf('[%.2f,%.2f]',fluxFPI.ancillary.pitchangle_edges(3:4));...
%                sprintf('[%.2f,%.2f]',fluxFPI.ancillary.pitchangle_edges(4:5))},[0.02 0.98])

%% Plot 2, scatter plot with fits, including pitch angle diagrams of average flux
tint_plot = tint;% + 6.5*[10 -10];
doFit = 1;
for mms_id = 1
  % antiparallel flux  
  c_eval('fluxFPI = ePitch?_pflux.palim([135 180]).elim(ediFlux?.ancillary.energy0).tlim(tint_plot);',mms_id)
  %c_eval('fluxFPI = ePitch?_pflux_rebin_edi.palim([135 180]).elim(ediFlux?.ancillary.energy0).tlim(tint_plot);',mms_id)
  c_eval('fluxEDI = ePitch?_pflux_edi.palim([135 180]).resample(fluxFPI).tlim(tint_plot);',mms_id)
  
  if doFit % Make linear fit, for each pitchangle, ignoring the zero and NaN values
    indNotZero = cell(1,4);
    indNotNaN = cell(1,4);
    
    ifit = 1;
    ft{ifit} = fittype({'x'}); ifit = ifit + 1;
    ft{ifit} = fittype({'x','1'}); ifit = ifit + 1;
    %ft{ifit} = fittype({'x^2','x'}); ifit = ifit + 1;
    %ft{ifit} = fittype({'x^2','x','1'}); ifit = ifit + 1;
    nfit = ifit - 1;
    curve_fit = cell(4,nfit);
    goodness_fit = cell(4,nfit);
    
    for ipitchangle = 1:4
      indNotZero{ipitchangle} = find(and(not(fluxFPI.data(:,ipitchangle)==0),not(fluxEDI.data(:,ipitchangle)==0)));
      indNotNaN{ipitchangle} = find(and(not(isnan(fluxFPI.data(:,ipitchangle))),not(isnan(fluxEDI.data(:,ipitchangle)))));
      indKeep = intersect(indNotZero{ipitchangle},indNotNaN{ipitchangle});
      fitFluxFPI = fluxFPI.data(indKeep,ipitchangle);
      fitFluxEDI = fluxEDI.data(indKeep,ipitchangle);
      
      for ifit = 1:nfit
        [curve_fit{ipitchangle,ifit}, goodness_fit{ipitchangle,ifit}] = fit( double(fitFluxFPI), double(fitFluxEDI), ft{ifit} );        
      end
    end
  end
  h = setup_subplots(2,4);
  for ipitchangle = 1:4
    hca = h(ipitchangle);
    scatter(hca,fluxFPI.data(:,ipitchangle),fluxEDI.data(:,ipitchangle),'.')    
    hca.XLabel.String = sprintf('flux FPI (%s)',fluxFPI.units);
    hca.YLabel.String = sprintf('flux EDI (%s)',fluxEDI.units);
    pitch_minus = fluxFPI.depend{2}(ipitchangle)-fluxFPI.ancillary.delta_pitchangle_minus(ipitchangle);
    pitch_plus = fluxFPI.depend{2}(ipitchangle)+fluxFPI.ancillary.delta_pitchangle_plus(ipitchangle);
    hca.Title.String = sprintf('MMS %g,\npitch angle = [%g %g] ',mms_id,pitch_minus,pitch_plus);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Box = 'on';
    if doFit
      hold(hca,'on')
      for ifit = 1:nfit
        str_legend = sprintf('y = %s',formula(curve_fit{ipitchangle,ifit}));
        replaceCoeffWithVals = 1;
        if replaceCoeffWithVals
          coeffs = coeffnames(curve_fit{ipitchangle,ifit});
          coeff_vals = coeffvalues(curve_fit{ipitchangle,ifit});
          for icoeff = 1:numcoeffs(curve_fit{ipitchangle,ifit})
            str_legend = strrep(str_legend,coeffs{icoeff},sprintf('%.1e',coeff_vals(icoeff)));
          end
        end
        fit_legends{ifit} = str_legend; %,coeffvalues(curve_fit{ipitchangle,ifit})
        coeffvals = coeffvalues(curve_fit{ipitchangle,ifit});
        ncoeff = numel(coeffvals);
        
        % Remove values that go outside the top of the box
        xfit = linspace(hca.XLim(1),hca.XLim(2),10);
        yfit = feval(curve_fit{ipitchangle,ifit},xfit);
        xfit(yfit>hca.YLim(2)) = NaN;
        yfit(yfit>hca.YLim(2)) = NaN;
        
        h_fit(ifit) = plot(hca,xfit,yfit);   
      end     
      legend(h_fit(1:nfit), fit_legends{1:nfit})
      hold(hca,'off')
    end
  end  
  %drawnow
  hl_scat = linkprop(h(1:4),{'XLim','YLim'});
%   xlim_max = 0;
%   c_eval('xlim_max = max([xlim_max max(hl_scat.Targets(?).Children.XData)]);',1:4)
%   ylim_max = 0;
%   c_eval('ylim_max = max([ylim_max max(hl_scat.Targets(?).Children.YData)]);',1:4)
  
  axis(h(1:4),'equal')
  h(1).XLim(1) = 0;
  h(1).YLim(1) = 0;
  %h(1).XLim(2) = xlim_max;
  %h(1).YLim(2) = ylim_max;
  
  isub = 5;
  hca = h(isub); isub = isub + 1;
  fluxFPI.plot_pad_polar(hca);
  %ePitch1.palim([135 180]).elim(500).plot_pad_polar(hca);
  hca.Title.String = 'FPI flux';
  
  hca = h(isub); isub = isub + 1;
  %dv_fpi/dv_edi
  fluxEDI.mtimes(1).plot_pad_polar(hca);
  hca.Title.String = 'EDI flux';
  
  hl_pad = linkprop(h(5:6),{'XLim','YLim','CLim'});
  clim_min = Inf;
  c_eval('clim_min = min([clim_min min(hl_pad.Targets(?).Children.CData)]);',1:2)
  clim_max = 0;
  c_eval('clim_max = max([clim_max max(hl_pad.Targets(?).Children.CData)]);',1:2)
  axis(h(5:6),'equal')
  h(5).CLim = [clim_min clim_max];
  %h(1).XLim(1) = 0;
  %h(1).YLim(1) = 0;
  %pause
end

%% Plot 2, scatter plot with fits, including pitch angle diagrams of average flux
tint_plot = tint;% + 6.5*[10 -10];
doFit = 1;
for mms_id = 1
  % antiparallel flux
  c_eval('fluxFPI = eFlux?.palim([135 180]).elim(ediFlux?.ancillary.energy0).tlim(tint_plot);',mms_id)
  c_eval('fluxEDI = ediFlux?.palim([135 180]).resample(fluxFPI).tlim(tint_plot);',mms_id)
  
  if doFit % Make linear fit, for each pitchangle, ignoring the zero and NaN values
    indNotZero = cell(1,4);
    indNotNaN = cell(1,4);
    
    ifit = 1;
    ft{ifit} = fittype({'x'}); ifit = ifit + 1;
    ft{ifit} = fittype({'x','1'}); ifit = ifit + 1;
    %ft{ifit} = fittype({'x^2','x'}); ifit = ifit + 1;
    %ft{ifit} = fittype({'x^2','x','1'}); ifit = ifit + 1;
    nfit = ifit - 1;
    curve_fit = cell(4,nfit);
    goodness_fit = cell(4,nfit);
    
    for ipitchangle = 1:4
      indNotZero{ipitchangle} = find(and(not(fluxFPI_rebin_edi.data(:,ipitchangle)==0),not(fluxEDI.data(:,ipitchangle)==0)));
      indNotNaN{ipitchangle} = find(and(not(isnan(fluxFPI_rebin_edi.data(:,ipitchangle))),not(isnan(fluxEDI.data(:,ipitchangle)))));
      indKeep = intersect(indNotZero{ipitchangle},indNotNaN{ipitchangle});
      fitFluxFPI = fluxFPI_rebin_edi.data(indKeep,ipitchangle);
      fitFluxEDI = fluxEDI.data(indKeep,ipitchangle);
      
      for ifit = 1:nfit
        [curve_fit{ipitchangle,ifit}, goodness_fit{ipitchangle,ifit}] = fit( double(fitFluxFPI), double(fitFluxEDI), ft{ifit} );        
      end
    end
  end
  h = setup_subplots(3,4);
  for ipitchangle = 1:4
    hca = h(ipitchangle);
    scatter(hca,fluxFPI.data(:,ipitchangle),fluxEDI.data(:,ipitchangle),'.')    
    hca.XLabel.String = sprintf('flux FPI (%s)',fluxFPI.units);
    hca.YLabel.String = sprintf('flux EDI (%s)',fluxEDI.units);
    pitch_minus = fluxFPI.depend{2}(ipitchangle)-fluxFPI.ancillary.delta_pitchangle_minus(ipitchangle);
    pitch_plus = fluxFPI.depend{2}(ipitchangle)+fluxFPI.ancillary.delta_pitchangle_plus(ipitchangle);
    hca.Title.String = sprintf('MMS %g,\npitch angle = [%g %g] ',mms_id,pitch_minus,pitch_plus);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Box = 'on';
    if doFit
      hold(hca,'on')
      for ifit = 1:nfit
        str_legend = sprintf('y = %s',formula(curve_fit{ipitchangle,ifit}));
        replaceCoeffWithVals = 1;
        if replaceCoeffWithVals
          coeffs = coeffnames(curve_fit{ipitchangle,ifit});
          coeff_vals = coeffvalues(curve_fit{ipitchangle,ifit});
          for icoeff = 1:numcoeffs(curve_fit{ipitchangle,ifit})
            str_legend = strrep(str_legend,coeffs{icoeff},sprintf('%.1e',coeff_vals(icoeff)));
          end
        end
        fit_legends{ifit} = str_legend; %,coeffvalues(curve_fit{ipitchangle,ifit})
        coeffvals = coeffvalues(curve_fit{ipitchangle,ifit});
        ncoeff = numel(coeffvals);
        
        % Remove values that go outside the top of the box
        xfit = linspace(hca.XLim(1),hca.XLim(2),10);
        yfit = feval(curve_fit{ipitchangle,ifit},xfit);
        xfit(yfit>hca.YLim(2)) = NaN;
        yfit(yfit>hca.YLim(2)) = NaN;
        
        h_fit(ifit) = plot(hca,xfit,yfit);   
      end     
      legend(h_fit(1:nfit), fit_legends{1:nfit})
      hold(hca,'off')
    end
  end
  for ipitchangle = 1:4
    hca = h(ipitchangle);
    scatter(hca,fluxFPI_rebin_edi.data(:,ipitchangle),fluxFPI_rebin_edi.data(:,ipitchangle),'.')    
    hca.XLabel.String = sprintf('flux FPI (%s)',fluxFPI_rebin_edi.units);
    hca.YLabel.String = sprintf('flux EDI (%s)',fluxFPI_rebin_edi.units);
    pitch_minus = fluxFPI_rebin_edi.depend{2}(ipitchangle)-fluxFPI_rebin_edi.ancillary.delta_pitchangle_minus(ipitchangle);
    pitch_plus  = fluxFPI_rebin_edi.depend{2}(ipitchangle)+fluxFPI_rebin_edi.ancillary.delta_pitchangle_plus(ipitchangle);
    hca.Title.String = sprintf('MMS %g,\npitch angle = [%g %g] ',mms_id,pitch_minus,pitch_plus);
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.Box = 'on';
    if doFit
      hold(hca,'on')
      for ifit = 1:nfit
        str_legend = sprintf('y = %s',formula(curve_fit{ipitchangle,ifit}));
        replaceCoeffWithVals = 1;
        if replaceCoeffWithVals
          coeffs = coeffnames(curve_fit{ipitchangle,ifit});
          coeff_vals = coeffvalues(curve_fit{ipitchangle,ifit});
          for icoeff = 1:numcoeffs(curve_fit{ipitchangle,ifit})
            str_legend = strrep(str_legend,coeffs{icoeff},sprintf('%.1e',coeff_vals(icoeff)));
          end
        end
        fit_legends{ifit} = str_legend; %,coeffvalues(curve_fit{ipitchangle,ifit})
        coeffvals = coeffvalues(curve_fit{ipitchangle,ifit});
        ncoeff = numel(coeffvals);
        
        % Remove values that go outside the top of the box
        xfit = linspace(hca.XLim(1),hca.XLim(2),10);
        yfit = feval(curve_fit{ipitchangle,ifit},xfit);
        xfit(yfit>hca.YLim(2)) = NaN;
        yfit(yfit>hca.YLim(2)) = NaN;
        
        h_fit(ifit) = plot(hca,xfit,yfit);   
      end     
      legend(h_fit(1:nfit), fit_legends{1:nfit})
      hold(hca,'off')
    end
  end
  linkprop(h(1:4),{'XLim','YLim'});
  axis(h(1:4),'equal')
  h(1).XLim(1) = 0;
  h(1).YLim(1) = 0;
  
  isub = 5;
  hca = h(isub); isub = isub + 1;
  fluxFPI.plot_pad_polar(hca);
  hca.Title.String = 'FPI flux';
  
  hca = h(isub); isub = isub + 1;
  %dv_fpi/dv_edi
  fluxEDI.mtimes(1).plot_pad_polar(hca);
  hca.Title.String = 'EDI flux';
  
  linkprop(h(5:6),{'XLim','YLim','CLim'});
  axis(h(5:6),'equal')
  %h(1).XLim(1) = 0;
  %h(1).YLim(1) = 0;
  %pause
end

%% Plot 2, scatter plot
fluxFPI = eFlux1.palim([135 180]).elim(ediFlux1.ancillary.energy0);
fluxEDI = ediFlux1.palim([135 180]).resample(fluxFPI);
h = subplot(1,1,1);
for ipitchangle = 1:4
  hca = h(1);
  hold(hca,'on')
  scatter(hca,fluxFPI.data(:,ipitchangle),fluxEDI.data(:,ipitchangle),'.')
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Box = 'on';
end
linkprop(h,{'XLim','YLim'});
axis(h,'equal')
h(1).XLim(1) = 0;
h(1).YLim(1) = 0;

%% Aid functions
function h = setup_subplots(nrows,ncols,orientation)

  orientation_default = 'horizontal';
  if not(exist('orientation','var'))
    orientation = orientation_default;
  elseif not(any(strcmp(orientation,{'horizontal','vertical'})))
    warning(sprintf('Orientation: %s, not recognized. Using default orientation %s',orientation,orientation_default))
    orientation = orientation_default;
  end

  npanels = nrows*ncols;
  ipanel = 0;
  for icol = 1:ncols
    for irow = 1:nrows
      ipanel = ipanel + 1;
      h(irow,icol) = subplot(nrows,ncols,ipanel);
    end
  end

  if strcmp(orientation,'vertical')
    h = h';
  end
  h = h(:);
end
