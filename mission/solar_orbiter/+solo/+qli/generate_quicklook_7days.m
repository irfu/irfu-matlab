function generate_quicklook_7days(Data, OutputPaths, Tint, logoPath)
%
% Generates ONE quicklook (file) for covering ONE UTC week of data.
%
%
% ARGUMENTS
% =========
% Data
%     Struct with various time series of data extracted from SPICE and datasets.
%     See the call from solo.qli.generate_quicklooks_all_types().
% OutputPaths
%     Struct with paths to separate output directories for the different types
%     of quicklooks (see solo.qli.generate_quicklooks_all_types).
% Tint
%     Should be a 7-day time interval consistent with the time series in "data"
%     e.g.
%     irf.tint('2020-06-03T00:00:00.00Z', '2020-06-10T00:00:00.00Z');
% logoPath
%     Either path to IRF logo, or empty.
%
%
% NOTES
% =====
% * The function obtains some data by reading CDF files directly (cdfread;
%   solo_L2_swa-pas-eflux).
% * The function also locates CDF files (reads files via SolO DB?) which are
%   used to obtain a time interval which is passed to solo.read_TNR(). ==>
%   Relies on more files.



tBeginSec = tic();



% Setup figure:
LINE_WIDTH       = 1.0;   % irf_plot() line width.
FONT_SIZE        = 18;    % Font size
LEGEND_FONT_SIZE = 22;    % irf_legend() font size.
COLORS           = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];

Units = irf_units;
Me    = Units.me;      % Electron mass [kg]
epso  = Units.eps0;    % Permitivitty of free space [Fm^-1]
qe    = Units.e;       % Elementary charge [C]

h            = irf_plot(9, 'newfigure');
fig          = gcf;
fig.Position = [1, 1, 1095, 800];



%===================================
% Fill panel 1: B vector components
%===================================
if ~isempty(Data.B)
  irf_plot(h(1), Data.B.tlim(Tint), 'linewidth', LINE_WIDTH);
  hold(    h(1), 'on');
  irf_plot(h(1), Data.B.abs.tlim(Tint), 'linewidth', LINE_WIDTH);
end
irf_legend(h(1), {'B_{R}', 'B_{T}', 'B_{N}', '|B|'}, [0.98 0.18], 'Fontsize', LEGEND_FONT_SIZE);
ylabel(    h(1), {'B_{RTN}';'(nT)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
irf_zoom(  h(1), 'y');

tBeginSec = solo.qli.utils.log_time('End panel 1', tBeginSec);



%======================
% Fill panel 2: abs(B)
%======================
if ~isempty(Data.B)
  %fci = qe*data.B.abs*10^-9/mp/(2*pi);
  irf_plot(h(2), Data.B.abs.tlim(Tint), 'linewidth', LINE_WIDTH);
end
ylabel(h(2), {'|B|';'(nT)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
h(2).YScale = 'log';
h(2).YTick  = [10, 100];
%h(2).YLim   = [0.1, 200];

tBeginSec = solo.qli.utils.log_time('End panel 2', tBeginSec);



%=========================
% Fill panel 3: Densities
%=========================
hold(h(3), 'on');
if ~isempty(Data.Ne)
  irf_plot(h(3), Data.Ne.tlim(Tint), 'color', COLORS(1,:), 'linewidth', LINE_WIDTH);
end
if ~isempty(Data.Npas)
  irf_plot(h(3), Data.Npas.tlim(Tint), 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
end
ylabel(    h(3), {'N';'(cm^{-3})'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
irf_legend(h(3), {'N_{e,RPW} ', ' N_{i,PAS}'}, [0.98 0.16], 'Fontsize', LEGEND_FONT_SIZE);
h(3).YScale = 'log';
h(3).YTick  = [10, 100];
%h(3).YLim   = [0.8, 200];

tBeginSec = solo.qli.utils.log_time('End panel 3', tBeginSec);



%===============================
% Fill panel 4: Ion temperature
%===============================
if ~isempty(Data.Tpas)
  irf_plot(h(4), Data.Tpas.tlim(Tint), 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
end
ylabel(h(4), {'T_i';'(eV)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
h(4).YScale = 'log';
h(4).YTick  = [1, 10, 100];
h(4).YLim   = [0.5, 300];

tBeginSec = solo.qli.utils.log_time('End panel 4', tBeginSec);



%==============
% Fill panel 5
%==============
% y,z PAS velocities
if ~isempty(Data.Vpas)
  irf_plot(h(5), Data.Vpas.y.tlim(Tint), 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
  hold(    h(5), 'on');
  irf_plot(h(5), Data.Vpas.z.tlim(Tint), 'color', COLORS(3,:), 'linewidth', LINE_WIDTH);
end
irf_legend(h(5), {'', 'v_{T}', 'v_{N}'}, [0.98 0.18], 'Fontsize', LEGEND_FONT_SIZE);
irf_zoom(  h(5), 'y');
ylabel(    h(5), {'v_{T,N}';'(km/s)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 5', tBeginSec);



%==============
% Fill panel 6
%==============
hold(h(6), 'on');
if ~isempty(Data.Vrpw)
  irf_plot(h(6), -Data.Vrpw, 'o', 'color', COLORS(1,:));
end
if ~isempty(Data.Vpas)
  irf_plot(h(6), Data.Vpas.x.tlim(Tint), 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
end
irf_legend(h(6), {'V_{RPW}', 'V_{PAS}'}, [0.98 0.15], 'Fontsize', LEGEND_FONT_SIZE);
%h(6).YLim=[150, 950];
ylabel(h(6), {'v_{R}'; '(km/s)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 6', tBeginSec);



%==============
% Fill panel 7
%==============
if ~isempty(Data.E)
  irf_plot(h(7), Data.E.y, 'color', COLORS(2,:), 'linewidth', LINE_WIDTH)
  hold(    h(7), 'on');
  %irf_plot(h(7), data.E.z, 'color', COLORS(3,:), 'linewidth', LWIDTH)
end
irf_legend(h(7), {'', 'E_y'}, [0.98 0.20], 'Fontsize', LEGEND_FONT_SIZE);
irf_zoom(  h(7), 'y');
ylabel(    h(7), {'E_{SRF}'; '(mV/m)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 7', tBeginSec);



%===================================
% Fill panel 8: Ion energy spectrum
% ---------------------------------
% NOTE: READS CDF FILES!
%===================================
if ~isempty(Data.ieflux)
  SwaFileArray = solo.db_list_files('solo_L2_swa-pas-eflux', Tint);
  iDEF   = struct('t', Data.ieflux.tlim(Tint).time.epochUnix);
  %for ii = 1:round((myFile(end).stop-myFile(1).start)/3600/24)
  for iFile = 1:length(SwaFileArray)
    iEnergy = cdfread(...
      fullfile(SwaFileArray(iFile).path, SwaFileArray(iFile).name), ...
      'variables', 'Energy');
    iEnergy = iEnergy{1};
    iDEF.p  = Data.ieflux.data;
  end
  iDEF.f       = repmat(iEnergy, 1, numel(iDEF.t))';
  iDEF.p_label = {'dEF', 'keV/', '(cm^2 s sr keV)'};
  irf_spectrogram(h(8), iDEF, 'log', 'donotfitcolorbarlabel');
  % set(h(1), 'ytick', [1e1 1e2 1e3]);
  hold(h(8), 'on');
  h8_clims = h(8).CLim;
  % Fix color axis
  h8_medp = mean(iDEF.p);
  h8_medp = min(h8_medp(h8_medp>0));
  if h8_medp > 0 && h8_medp > h8_clims(1) && log10(h8_medp)+2<(max(max(log10(iDEF.p))))
    caxis(h(8), [log10(h8_medp)+2 (max(max(log10(iDEF.p))))])
  end
  set(     h(8), 'YScale', 'log');
  colormap(h(8), jet)
  ylabel(  h(8), '[eV]')
end

tBeginSec = solo.qli.utils.log_time('End panel 8', tBeginSec);



%======================================
% Fill panel 9: E-field spectrum (TNR)
% ------------------------------------
% NOTE: READS CDF FILES!
%======================================
% NOTE: Panel takes much more time than other panels.
if ~isempty(Data.Etnr)
  % Electron plasma frequency
  TnrFileArray = solo.db_list_files('solo_L2_rpw-tnr-surv-cdag', Tint);
  tp = [];
  pp = [];
  warning('off', 'fuzzy:general:warnDeprecation_Combine');
  TNR = [];
  %for iii = 1:round((myFile2(end).stop-myFile2(1).start)/3600/24)
  for iFile = 1:length(TnrFileArray)
    tt     = [TnrFileArray(iFile).start, TnrFileArray(iFile).stop];
    [TNRp] = solo.read_TNR(tt);
    if isa(TNRp, 'struct')
      % NOTE: MATLAB documentation (R2019b):
      % "combine will be removed in a future release"
      TNR.t = combine(tp, TNRp.t);
      tp    = TNR.t;
      TNR.p = combine(pp, TNRp.p);
      pp    = TNR.p;

      % IMPLEMENTATION NOTE: Only read from TNRp from within this if
      % clause, since it might not be a struct if read from elsewhere,
      % even if it in principle means overwriting the value multiple
      % times as for TNRp.f and TNRp.p_label.
      TNR.f       = TNRp.f;
      TNR.p_label = TNRp.p_label;
    end
  end
  if isstruct(TNR)
    % TNR.f       = TNRp.f;
    % TNR.p_label = TNRp.p_label;
    sz_tnr = size(TNR.p);
    if sz_tnr(1) == length(TNR.t) && sz_tnr(2) == length(TNR.f)
      irf_spectrogram(h(9), TNR, 'log', 'donotfitcolorbarlabel')
      hold(           h(9), 'on');
    end
    if ~isempty(Data.Ne)
      wpe_sc       = (sqrt(((Data.Ne.tlim(Tint)*1000000)*qe^2)/(Me*epso)));
      fpe_sc       = (wpe_sc/2/pi)/1000;
      fpe_sc.units = 'kHz';
      fpe_sc.name  = 'f [kHz]';
      irf_plot(h(9), fpe_sc, 'r', 'linewidth', LINE_WIDTH);
    end
    text(    h(9), 0.01, 0.3, 'f_{pe,RPW}', 'units', 'normalized', 'fontsize', FONT_SIZE, 'Color', 'r');
    %set(h(9), 'YScale', 'log');
    colormap(h(9), jet)
    %ylabel(h(9), 'f [kHz]')
    set(     h(9), 'ColorScale', 'log')
    %caxis([.01 10]*10^-12)
    yticks(  h(9), [10^1 10^2]);
  end
end

tBeginSec = solo.qli.utils.log_time('End panel 9', tBeginSec);



%======================
% Other, miscellaneous
%======================
irf_plot_axis_align(h(1:9));
irf_zoom(h(1:9), 'x', Tint);
irf_zoom(h(1),   'y');
irf_zoom(h(5:9), 'y');

h(2).YLabel.Position = [1.05, 0.5, 0];
%yyaxis(h(2), 'left');
h(2).YLabel.Units    = 'normalized';
h(2).YLabel.Position = h(3).YLabel.Position;

h(9).XLabel.Visible = 'off';

% Add spacecraft position as text.
[soloStr, earthStr] = solo.qli.utils.get_context_info_strings(Data.soloPos, Data.earthPos, Tint);
text(h(9), -0.11, -0.575, soloStr, 'units', 'normalized', 'fontsize', FONT_SIZE);
% Add Earth longitude as text.
text(h(9), -0.11, -0.925, earthStr, 'units', 'normalized', 'fontsize', FONT_SIZE);



xtickangle(h(9), 0)

%======================================================
% Add IRF logo and data source information info string
%======================================================
logoPos = h(1).Position;    %  [left, bottom, width, height]
logoPos(1) = logoPos(1) + logoPos(3) + 0.01;
logoPos(2) = logoPos(2) + 0.06;
logoPos(3) = 0.05;
logoPos(4) = logoPos(3) * 1095/800;
hLogoAxes = axes('position', logoPos);
if ~isempty(logoPath)
  [x, ~] = imread(logoPath);
  image(x)
end
% colormap (map)
set(hLogoAxes, 'handlevisibility', 'off', 'visible', 'off')

str = solo.qli.utils.get_data_source_info_string();
text(h(1), 0, 1.2, str, 'Units', 'normalized')



%===============
% Adjust panels
%===============
% Remove overlapping ticks.
solo.qli.utils.ensure_axes_data_tick_margins(h)

%yyaxis(h(2), 'left');
%oldlims2 = h(2).YLim;
%oldticks2 = h(2).YTick;
h(2).YScale = 'log';
h(2).YTick  = [1, 10, 100];
h(2).YLim   = [0.8, 200];

% yyaxis(h(2), 'right');
% oldlims2_r=h(2).YLim;
% oldticks2_r = h(2).YTick;
% h(2).YScale='log';
% h(2).YTick=[1, 10, 100];
%h(2).YLim=[0.1, 200];

%oldlims5 = h(5).YLim;
%oldticks5 = h(5).YTick;
h(5).YScale = 'log';
h(5).YTick  = [1, 10, 100];
%h(5).YLim   = [0.5, 300];

% c_eval('h(?).FontSize=18;', 1:9);
for i = 1:9
  h(i).FontSize = FONT_SIZE;
end


irf_plot_axis_align(h(1:9));
irf_zoom(           h(1:9), 'x', Tint);
% irf_zoom(h(1:7), 'y');

% Plot complete, print figure.
fig = gcf;
fig.PaperPositionMode = 'auto';

%=====================
% Save figure to file
%=====================
solo.qli.utils.save_figure_to_file(OutputPaths.path_1w, Tint)
% TODO-NI: Why are there any commands (except close()) after this?
%          Did the code use to iterate over 24h, 6h, 2h plots too?

% h(2).YScale='lin';
% h(2).YTick=oldticks2_r;
% h(2).YLim=oldlims2_r;
%yyaxis(h(2), 'left');
% h(2).YScale='lin';
% h(2).YLim=oldlims2;
% h(2).YTick=oldticks2;
%
% h(5).YScale='lin';
% h(5).YLim=oldlims5;
% h(5).YTick=oldticks5;

close(fig);

[~] = solo.qli.utils.log_time('End of generate_quicklook_7days.m', tBeginSec);

end
