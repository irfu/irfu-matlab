function generate_quicklooks_24h_6h_2h(Data, OutputPaths, Tint24h, logoPath)
%
% Generates quicklooks (files) for covering ONE UTC day of data:
% 1x24h quicklooks, 4x6h quicklooks, 12x2h quicklooks.
%
%
% ARGUMENTS
% =========
% Data
%     Struct with various time series of data extracted from SPICE and datasets.
%     See the call from solo.qli.generate_quicklooks_24h_6h_2h_using_DB_SPICE().
% OutputPaths
%     Struct with paths to separate output directories for the different types
%     of quicklooks.
%     Has fields: .dir2h, .dir6h, .dir24h
% Tint24h
%     Should be a 24-hour time interval consistent with the time series in
%     "data", e.g.
%     irf.tint('2020-06-01T00:00:00.00Z', '2020-06-02T00:00:00.00Z');
% logoPath
%     Either path to IRF logo, or empty.
%
%
% NOTES
% =====
% * Computes the spectrum for B when magnetic field data is available. When it
%   does, the code takes a lot of time.



% ==========
% KNOWN BUGS
% ==========
% NOTE: Several bug descriptions here may obsolete.
%
% BUG?: Panel 2/density/abs(B): Sometimes has no left-hand ticks (for density?).
%   /EJ 2023-05-10
%   Ex: 20220329T04_20220329T06.png
%   EJ 2023-05-11: Should be fixed.
% BUG?: Panel 6: V_T, V_N: Y limits seem bad (recently gotten worse from earlier implementation).
%   /EJ 2023-05-10
%   Ex: 2022-03-23T10-T12. -- Seems wrong timestamp.
%   Mistakenly looked at wrong plots?
%
% BUG: Cirka dfe637c8 (2023-05-12 11:02:05 +0200)
%   Panels can get extra wide, and IRF logo is partly outside image. Presumably
%   when there is no color bar (due to missing data).
%   NOTE: Might not be able to observe bug later if relevant data gaps are later
%         filled in.
%   Ex: 6h plots, as of 2023-05-23:
%      2023-02-05 18-00: Normal. Has one colorbar for "f (kHz)"
%      2023-02-06 00-06: Wider panels. Has no colorbar for "f (kHz)"
%   Ex: 24h plots, as of 2023-05-23:
%      2023-02-05: Normal. Has one colorbar for "f (kHz)"
%      2023-02-06: Wider panels. Has no colorbar for "f (kHz)"
%
% BUG: 24h, panel 2, |B| (right axis) is scaled badly on y axis. Too much extra
% space. Old quicklooks were better.
% a44b3127 Erik P G Johansson (2024-03-21 12:56:51 +0100) SolO QLI: Change terms: official {processing-->generation}
% Ex: 2023-01-05
% irf_zoom(h(2), 'y', [minAbsB-1, maxAbsB+1]); at plotting gives good zoom, but
% set_YLim_YTick(h([]), h([2]), h([])) later re-zooms in worse way.
% /Erik P G Johansson 2024-03-21
% Should be solved.
% /Erik P G Johansson 2024-03-25
%
% BUG?: 24h, panel 10, TNR spectrum (irf_spectrogram()):
% Spectral data (irfu-matlab version 2024-03-22) freuently looks different
% compared to earlier spectras (cron job, irfu-matlab version circa 2023-04).
% More sections with constant spectrum, some new rectangular "holes" (white,
% unfilled) in spectra. Some sections are unaltered. Unclear if this si due to
% changes in underlying data or in code (e.g. recent change in solo.read_TNR()).
% Ex: 24h quicklook for 2023-01-27:
%   Large part of spectrum has changed from time changing to constant.
%   New "hole" in spectrum.
%   Last section is the same as before.
% /Erik P G Johansson 2024-03-25
%
% BUG: There is no date label for the data under the lowest panel, if there is
% no data. Do not think this was a problem previously(?).
% NOTE: irf_spectrogram() sets the date label (sic!). There may also be other
% functions which also set it.
% Ex: 2023-03-24, solo.qli.generate_quicklooks_24h_6h_2h___UTEST.test_no_data().
% /Erik P G Johansson 2024-03-27



% TODO-NI: Panels 2 & 5 are logarithmic for 24h plots and linear for 6h & 2h?
%          Should they be?
% TODO-NI: Old panels 5 had a constant y axis range (YLim) for 24h, but dynamic
%          (changing) for 6h & 2h. Was that intentional?
% TODO-NI Panel 10 (log) is hardcoded to YLim~[10, 100] (because that is what
%         it used to be). This does not cover the entire interval of data
%         (there is more data at lower y). Should it be that way?
%
% PROPOSAL: Eliminate isempty(Data.tnrBand).
%   PRO: Seems unnecessary.



tBeginSec = tic();



LINE_WIDTH       = 1.0;   % irf_plot() line width.
FONT_SIZE        = 18;    % Font size
LEGEND_FONT_SIZE = 22;    % irf_legend() font size.
% Different colors to reuse. One color per row.
COLORS           = [
  0 0 0;    % Black
  0 0 1;    % Blue
  1 0 0;    % Red
  0 0.5 0;
  0 1 1;
  1 0 1;
  1 1 0];

% NOTE: Unclear units. Magnitude implies pixels, but actual quicklooks are
% 2281x1667 pixels.
FIG_WIDTH  = 1095;
FIG_HEIGHT =  800;

UNITS = irf_units;
Me    = UNITS.me;      % Electron mass [kg]
eps0  = UNITS.eps0;    % Permittivity of free space [Fm^-1]
mp    = UNITS.mp;      % Proton mass [km]
qe    = UNITS.e;       % Elementary charge [C]

% Setup figure
h            = irf_plot(10, 'newfigure');
fig          = gcf;
fig.Position = [1, 1, FIG_WIDTH, FIG_HEIGHT];



%fig.set('units', 'normalized', 'outerposition', [0.5 0 0.5 1])  % DEBUG



%===================================
% Fill panel 1: B vector components
%===================================
if ~isempty(Data.B)
  irf_plot(h(1), Data.B.tlim(Tint24h), 'linewidth', LINE_WIDTH);
  hold(    h(1), 'on');
  irf_plot(h(1), Data.B.abs.tlim(Tint24h), 'linewidth', LINE_WIDTH);
end
irf_legend(h(1), {'B_{R}', 'B_{T}', 'B_{N}', '|B|'}, [0.98 0.18], 'Fontsize', LEGEND_FONT_SIZE);
ylabel(    h(1), {'B_{RTN}'; '(nT)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 1', tBeginSec);



%=======================
% Fill panel 2: N & |B|
%=======================
hold(h(2), 'on');
if ~isempty(Data.Ne)
  irf_plot(h(2), Data.Ne.tlim(Tint24h), '-', 'color', COLORS(1,:), 'linewidth', LINE_WIDTH);
end
if ~isempty(Data.Npas)
  irf_plot(h(2), Data.Npas.tlim(Tint24h), '-', 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
end
ylabel(h(2), {'N'; '(cm^{-3})'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
h(2).ColorOrder = COLORS;
irf_legend(h(2), {'N_{e,RPW}', 'N_{i,PAS}', '|B|'}, [0.98 0.16], 'Fontsize', LEGEND_FONT_SIZE);

yyaxis(h(2), 'right');
ABS_B_COLOR = COLORS(3,:);
if ~isempty(Data.B)
  irf_plot(h(2), Data.B.abs.tlim(Tint24h), 'color', ABS_B_COLOR, 'linewidth', LINE_WIDTH);
  %Bnan = rmmissing(data.B.abs.data);
  %if ~isempty(Bnan)
  %    h(2).YLim = [floor(min(abs(Bnan))),ceil(max(abs(Bnan)))];
  %end
  minAbsB = min(Data.B.tlim(Tint24h).abs.data);
  maxAbsB = max(Data.B.tlim(Tint24h).abs.data);
  if ~isnan(minAbsB) && ~isnan(maxAbsB)
    % Only zoom if min & max are not NaN (==> Avoid crash).
    % NOTE: NaN should only happen if there is NaN in the data. Absence of data
    % should yield 0x1 array.
    irf_zoom(h(2), 'y', [minAbsB-1, maxAbsB+1]);
  end
else
  % Values in the absence of data.
  minAbsB = NaN;
  maxAbsB = NaN;
end
ylabel(h(2), {'|B|'; '(nT)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
h(2).YColor = ABS_B_COLOR;

tBeginSec = solo.qli.utils.log_time('End panel 2', tBeginSec);



%==========================================
% Fill panel 3 & 4: Spectra derived from B
%==========================================
IRF_EBSP_FREQ_MIN_HZ            = 0.05;
B_SAMPLING_PERIOD_THRESHOLD_SEC = 0.1250*0.95;  % 0.1250*0.95 = 0.1187
if ~isempty(Data.B) && solo.qli.const.B_SPECTRA_ENABLED
  if  ~isempty(rmmissing(Data.B.data))
    B = Data.B;

    fci = qe*(B.abs*10^-9)/mp/(2*pi);    % Proton gyration frequency [cycles/s]

    medianSamplingPeriodSec = median(diff((B.time.epochUnix)));
    if medianSamplingPeriodSec < B_SAMPLING_PERIOD_THRESHOLD_SEC
      fMag = 128; freqMaxHz = 7;
    else
      fMag =   8; freqMaxHz = 3;
    end

    % Create filtered version of B.
    %  TSeries.filt() --> irf_filt() --> filtfilt().
    B_0 = B.filt(0, 0.01, fMag, 5);

    %---------------------------------------------------------------
    % IMPORTANT NOTE: The call to irf_ebsp() is very time-consuming
    %---------------------------------------------------------------
    tBeginSec = solo.qli.utils.log_time('irf_ebsp(): Begin call', tBeginSec);
    Ebsp      = irf_ebsp([], B, [], B_0, [], [IRF_EBSP_FREQ_MIN_HZ, freqMaxHz], 'fullB=dB', 'polarization', 'fac');
    tBeginSec = solo.qli.utils.log_time('irf_ebsp(): End call', tBeginSec);

    frequency   = Ebsp.f;
    time        = Ebsp.t;
    Bsum        = Ebsp.bb_xxyyzzss(:, :, 4);
    ellipticity = Ebsp.ellipticity;
    dop         = Ebsp.dop;    % DOP = Degree Of Polarization

    % Remove points with very low degree of polarization
    DEGREE_OF_POLARIZATION_THRESHOLD = 0.7;
    % iRemove = find(dop < DEGREE_OF_POLARIZATION_THRESHOLD);
    bRemove = dop < DEGREE_OF_POLARIZATION_THRESHOLD;
    ellipticity(bRemove) = NaN;

    % Remove "lonely" pixels
    msk              = ellipticity;
    msk(~isnan(msk)) = 1;
    msk( isnan(msk)) = 0;
    % "BW2 = bwareaopen(BW,P) removes from a binary image all connected
    % components (objects) that have fewer than P pixels, producing another
    % binary image BW2."
    msk_denoise                 = bwareaopen(msk, 8);
    ellipticity(msk_denoise==0) = NaN;

    %-----------------
    % Panel 3: "Bsum"
    %-----------------
    Specrec         = struct('t', time);
    Specrec.f       = frequency;
    Specrec.p       = Bsum;
    Specrec.f_label = '';
    Specrec.p_label = {'log_{10}B^{2}', 'nT^2 Hz^{-1}'};
    irf_spectrogram(h(3), Specrec, 'log', 'donotfitcolorbarlabel');
    set(     h(3), 'yscale', 'log');
    % set(h(1), 'ytick', [1e1 1e2 1e3]);
    % caxis(h(3), [-8 -1])
    hold(    h(3), 'on');
    irf_plot(h(3), fci, 'k', 'linewidth', LINE_WIDTH);
    text(    h(3), 0.01, 0.3, 'f_{ci}', 'units', 'normalized', 'fontsize', FONT_SIZE);
    colormap(h(3), 'jet');

    %----------------------
    % Panel 4: Ellipticity
    %----------------------
    Specrec         = struct('t', time);
    Specrec.f       = frequency;
    Specrec.p       = ellipticity;
    Specrec.f_label = '';
    Specrec.p_label = {'Ellipticity', 'DOP>0.7'};
    irf_spectrogram(h(4), Specrec, 'log', 'donotfitcolorbarlabel');
    set(     h(4), 'yscale', 'log');
    % set(h(1), 'ytick', [1e1 1e2 1e3]);
    caxis(   h(4), [-1 1])
    hold(    h(4), 'on');
    irf_plot(h(4), fci, 'k', 'linewidth', LINE_WIDTH);
    text(    h(4), 0.01, 0.3, 'f_{ci}', 'units', 'normalized', 'fontsize', FONT_SIZE);

    crr     = interp1([1 64 128 192 256], [0.0  0.5 0.75 1.0 0.75], 1:256);
    cgg     = interp1([1 64 128 192 256], [0.0  0.5 0.75 0.5 0.00], 1:256);
    cbb     = interp1([1 64 128 192 256], [0.75 1.0 0.75 0.5 0.00], 1:256);
    bgrcmap = [crr' cgg' cbb'];
    colormap(h(4), bgrcmap);
  end
end
ylabel(h(3), {'f'; '(Hz)'}, 'fontsize', FONT_SIZE);
ylabel(h(4), {'f'; '(Hz)'}, 'fontsize', FONT_SIZE);
tBeginSec = solo.qli.utils.log_time('End panel 3 & 4', tBeginSec);



%===============================
% Fill panel 5: Ion temperature
%===============================
if ~isempty(Data.Tpas)
  irf_plot(h(5), Data.Tpas.tlim(Tint24h), 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
end
irf_zoom(h(5), 'y');
ylabel(  h(5), {'T_i'; '(eV)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 5', tBeginSec);



%==================================
% Fill panel 6: y,z PAS velocities
%==================================
if ~isempty(Data.Vpas)
  irf_plot(h(6), Data.Vpas.y.tlim(Tint24h), 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
  hold(    h(6), 'on');
  irf_plot(h(6), Data.Vpas.z.tlim(Tint24h), 'color', COLORS(3,:), 'linewidth', LINE_WIDTH);
end
irf_legend(h(6), {'', 'v_{T}', 'v_{N}'}, [0.98 0.18], 'Fontsize', LEGEND_FONT_SIZE);
irf_zoom(  h(6), 'y');
ylabel(    h(6), {'V_{T,N}'; '(km/s)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 6', tBeginSec);



%=====================================
% Fill panel 7: Vrpw, Vpas velocities
%=====================================
hold(h(7), 'on');
if ~isempty(Data.Vrpw)
  irf_plot(h(7),-Data.Vrpw, 'o-', 'color', COLORS(1,:));
end
if ~isempty(Data.Vpas)
  irf_plot(h(7), Data.Vpas.x.tlim(Tint24h), 'color', COLORS(2,:), 'linewidth', LINE_WIDTH);
end
irf_legend(h(7), {'V_{RPW}', 'V_{PAS}'}, [0.98 0.18], 'Fontsize', LEGEND_FONT_SIZE);
irf_zoom(  h(7), 'y');
ylabel(    h(7), {'V_R'; '(km/s)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 7', tBeginSec);



%==============================
% Fill panel 8: Electric field
%==============================
if ~isempty(Data.E)
  irf_plot(h(8), Data.E.y, 'color', COLORS(2,:), 'linewidth', LINE_WIDTH)
  hold(    h(8), 'on');
  %irf_plot(h(8), data.E.z, 'color', COLORS(3,:), 'linewidth', LINE_WIDTH)

  minEy = min(rmmissing(Data.E.y.data));
  maxEy = max(rmmissing(Data.E.y.data));
  if ~isempty(minEy) && ~isempty(maxEy)
    irf_zoom(h(8), 'y', [minEy-5 maxEy+5]);
  end
end
irf_legend(h(8), {'', 'E_y'}, [0.98 0.20], 'Fontsize', LEGEND_FONT_SIZE);
ylabel(    h(8), {'E_{SRF}'; '(mV/m)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 8', tBeginSec);



%============================================================================
% Fill panel 9: Ion energy spectrum
% ---------------------------------
% NOTE: READS CDF FILES! -- OBSOLETE INFO. REFACTORED AWAY.
% NOTE: Essentially the same as solo.qli.generate_quicklook_7days(): Panel 8
%============================================================================
if ~isempty(Data.ieflux)
  %   SwaFileArray = solo.db_list_files('solo_L2_swa-pas-eflux', Tint24h);
  iDEF         = struct('t', Data.ieflux.tlim(Tint24h).time.epochUnix);
  %for ii = 1:round((myFile(end).stop-myFile(1).start)/3600/24)
  %   for iFile = 1:length(SwaFileArray)
  %     % NOTE: Reads CDFs using cdfread() which is a MATLAB function (i.e. not
  %     %       dataobj(), not spdfcdfread()).
  %     % NOTE: zVariable "Energy" seems to be metadata (not science data).
  %     %       zVariable attributes CATDESC="Center of energy bins",
  %     %       VAR_TYPE="support_data". No DEPEND_0, so not time-dependent.
  %     % NOTE: Can not load this variable using
  %     %       solo.qli.utils.db_get_ts('solo_L2_swa-pas-eflux', 'eflux', Tint);
  %     %       Gets error message: "Data does not contain DEPEND_0 or DATA"
  %     iEnergy = cdfread(...
  %       fullfile(SwaFileArray(iFile).path, SwaFileArray(iFile).name), ...
  %       'variables', 'Energy');
  %     iEnergy = iEnergy{1};
  %   end
  iEnergy      = Data.swaEnergyMetadata;
  iDEF.p       = Data.ieflux.data;
  iDEF.p_label = {'dEF', 'keV/', '(cm^2 s sr keV)'};
  iDEF.f       = repmat(iEnergy, 1, numel(iDEF.t))';
  irf_spectrogram(h(9), iDEF, 'log', 'donotfitcolorbarlabel');   % NOTE: Somewhat time-consuming.
  % set(h(1), 'ytick', [1e1 1e2 1e3]);
  %caxis(h(9), [-1 1])

  hold(h(9), 'on');
  h9_clims = h(9).CLim;
  % Fix color axis
  h9_medp       = mean(iDEF.p);              % MxN --> 1xN
  h9_medp       = min(h9_medp(h9_medp>0));
  h9_caxisRange = [log10(h9_medp)+2, max(max(log10(iDEF.p)))];
  %if (h9_medp > 0) && (h9_medp > h9_clims(1)) && (log10(h9_medp)+2 < max(max(log10(iDEF.p))))
  if (h9_medp > 0) && (h9_medp > h9_clims(1)) && (h9_caxisRange(1) < h9_caxisRange(2))
    %caxis(h(9), [log10(h9_medp)+2 max(max(log10(iDEF.p)))])
    caxis(h(9), h9_caxisRange)
  end
end
set(     h(9), 'YScale', 'log');
colormap(h(9), jet)
ylabel(  h(9), {'W_{i}'; '(eV)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);

tBeginSec = solo.qli.utils.log_time('End panel 9', tBeginSec);



%=======================================================
% Fill panel 10: E-field spectrum (TNR)
% -------------------------------------
% NOTE: READS CDF FILES indirectly via solo.read_TNR()! -- OBSOLETE INFO. REFACTORED AWAY.
%=======================================================
% BUG(?): PROBABLY WHAT HAPPENS: Does not create color bar (no call to
% colormap()) when there is no data.
% ==> The panel becomes wider.
% ==> Other panels become wider.
% ==> Moves the IRF logo to the right, and partially outside image.
if ~isempty(Data.tnrBand)
  %   try
  %     [TNR] = solo.read_TNR(Tint24h);
  %   catch Exc
  %     if strcmp(Exc.identifier, 'read_TNR:FileNotFound')
  %       TNR = [];
  %     end
  %   end

  if isstruct(Data.Tnr)
    sz_tnr = size(Data.Tnr.p);
    % NOTE: Condition is a way of determining whether Data.Tnr contains
    %       any real data, despite not being empty.
    if sz_tnr(1) == length(Data.Tnr.t) && sz_tnr(2) == length(Data.Tnr.f)
      irf_spectrogram(h(10), Data.Tnr, 'log', 'donotfitcolorbarlabel')
      hold(           h(10), 'on');
      if ~isempty(Data.Ne)
        % Electron plasma frequency
        wpe_sc = (sqrt(((Data.Ne.tlim(Tint24h)*1000000)*qe^2)/(Me*eps0)));   % TSeries
        fpe_sc = (wpe_sc/2/pi)/1000;                                         % TSeries --> TSeries (sic!)
        irf_plot(h(10), fpe_sc, 'r', 'linewidth', LINE_WIDTH);
        fpe_sc.units = 'kHz';
        fpe_sc.name  = 'f [kHz]';
      end
      hold(h(10), 'off');
      text(h(10), 0.01, 0.3, 'f_{pe,RPW}', 'units', 'normalized', 'fontsize', FONT_SIZE, 'Color', 'r');
      set( h(10), 'YScale', 'log');
      %set(h(10), 'ColorScale', 'log')
      %caxis(h(10), [.01 1]*10^-12)
      ylabel(h(10), {'f'; '(kHz)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
      colormap(h(10), jet)
      %yticks(h(10), [10^1 10^2]);
      %irf_zoom(h(10), 'y', [10^1 10^2])
    end
  end
end
ylabel(h(10), {'f'; '(kHz)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
yticks(h(10),      [10^1 10^2]);
% Not set YLim.
%irf_zoom(h(10), 'y', [10^1 10^2])
irf_zoom(h(10), 'y', [9, 110])    % Not overwritten later.

if isempty(Data.Vrpw) ...
    && isempty(Data.E)    && isempty(Data.Ne)   && isempty(Data.B) ...
    && isempty(Data.Tpas) && isempty(Data.Npas) && isempty(Data.ieflux) ...
    && isempty(Data.tnrBand)
  nanPlot = irf.ts_scalar(Tint24h,ones(1, 2)*NaN);
  irf_plot(h(10), nanPlot);    % No LINE_WIDTH?
  grid(    h(10), 'off');
  ylabel(  h(10), {'f'; '(kHz)'}, 'interpreter', 'tex', 'fontsize', FONT_SIZE);
end

tBeginSec = solo.qli.utils.log_time('End panel 10', tBeginSec);



%======================
% Other, miscellaneous
%======================
irf_plot_axis_align(h(1:10));      % Make the panels have the same width.
irf_zoom(h(1:10), 'x', Tint24h);   % Make the panels cover the same x (time) range.
irf_zoom(h(1),    'y');

% Correct the position of the right y label.
% (It is affected by irf_plot_axis_align(h(1:10))).
yyaxis(h(2), 'right');
h(2).YLabel.Position = [1.05, 0.5, 0];

% Correct/adjust the position of the left y label, so that it lines up with
% another panel's y label
% -------------------------------------------------------------------------
% NOTE: *NOT* using h(3) or h(4) since they do not have ylabels if the relevant
% data is missing. h(1) ylabel always has a position. Using h(3) (old
% implementation) lead to left panel 2 ylabel having the wrong position (too far
% left) when h(3) did not have any label).
yyaxis(h(2), 'left');
h(2).YLabel.Units    = h(1).YLabel.Units;
h(2).YLabel.Position = h(1).YLabel.Position;

% Add Context Info Strings (CIS): Spacecraft position, Earth longitude as text.
[soloStr, earthStr] = solo.qli.utils.get_context_info_strings(Data.soloPos, Data.earthPos, Tint24h);
hCisText1 = text(h(10), -0.11, -0.575, soloStr,  'units', 'normalized', 'fontsize', FONT_SIZE);
hCisText2 = text(h(10), -0.11, -0.925, earthStr, 'units', 'normalized', 'fontsize', FONT_SIZE);



xtickangle(h(10), 0)

%======================================================
% Add IRF logo and data source information info string
%======================================================
panelPos = h(1).Position;    %  [left, bottom, width, height]
logoPos = [
  panelPos(1) + panelPos(3) + 0.06;
  panelPos(2) + 0.06;
  0.05;
  0.05 * FIG_WIDTH/FIG_HEIGHT;
  ];
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
% IMPLEMENTATION NOTE: Not sure why these commands are located here rather
% than where their respective panels are created. Related to
% irf_plot_axis_align()?
yyaxis(h(2), 'left');
h(2).YScale = 'log';       % NOTE: Later changed to LIN for non-24h quicklooks.
h(2).YTick  = [1, 10, 100];

% Remove overlapping ticks.
% Automatically set YLim+YTick, or automatically set YLim, or adjust YLim,
% depending on panel.
%set_YLim_YTick(h([1, 6:9]), h([2]), h([3:5, 10]))
solo.qli.utils.set_YLim_YTick_automatically( h([1, 6:9 ]))
%solo.qli.utils.set_YLim_YTick_automatically( h([2])      )
solo.qli.utils.ensure_axes_data_tick_margins(h([3:5, 10]))



% NOTE: Had difficulties making panel 2, right axis have a sensible range and
% ticks. Therefore using explicit min & max from data. In part because MATLAB
% sets YLim min=0 which is bad for a log scale.
yyaxis(h(2), 'right');
h(2).YScale = 'log';       % NOTE: Later changed to LIN for non-24h quicklooks.
h(2).YTick  = [1, 10, 100];
if ~isnan(minAbsB) && ~isnan(maxAbsB)
  % NOTE: Log scale makes effective visual distance smaller. May need greater
  % deviation from 1 than intuitively expected.
  C_ABS_B_MARGIN = 1.2;
  h(2).YLim = [minAbsB/C_ABS_B_MARGIN, maxAbsB*C_ABS_B_MARGIN];
end

solo.qli.utils.ensure_axes_data_tick_margins(h(2))


% NOTE: h(5).YLim are hardcoded and seem too broad/wide.
% PROPOSAL: Not overwrite automatic YLim?
oldYLimH5   = h(5).YLim;
oldYTickH5  = h(5).YTick;
h(5).YScale = 'log';       % NOTE: Later changed to LIN for non-24h quicklooks.
h(5).YTick  = [1, 10, 100];
h(5).YLim   = [0.5, 300];



fig.PaperPositionMode = 'auto';

%===========================
% Save figure to file (24h)
%===========================
% PROPOSAL: Use modify_save_subinterval_plot() here?
%   NOTE: Indirectly calls solo.qli.utils.ensure_axes_data_tick_margins() which
%   thus should be disabled above.
%   ==> Changing order of commands.
%   ==> Calls to YScale, YTick above become superseded.
%   ==> Unwanted change of behaviour.
solo.qli.utils.save_figure_to_file(OutputPaths.dir24h, Tint24h)



%=============================================
% Modify panels 2 & 5, AFTER saving 24 h plot
%=============================================
% Change panel 2+5 y scales to "lin" (previously "log").
% h(5): Keep old ylimits and ticks!
yyaxis(h(2), 'right');
h(2).YScale    = 'lin';       % NOTE: Previously LOG.
h(2).YTickMode = 'auto';
yyaxis(h(2), 'left');
h(2).YScale    = 'lin';       % NOTE: Previously LOG.
h(2).YTickMode = 'auto';

h(5).YScale = 'lin';          % NOTE: Previously LOG.
h(5).YLim   = oldYLimH5;
h(5).YTick  = oldYTickH5;



I_6H = 0:3;
I_2H = 0:11;
if ~solo.qli.const.NONWEEKLY_ALL_PLOTS_ENABLED
  % For debugging/testing.
  I_6H = [0];
  I_2H = [0];
  %I_6H = [1];
  %I_2H = [5];
end

%===========================
% Iterate over 6h intervals
%===========================
tBeginSec = solo.qli.utils.log_time('Begin iterating over 6 h intervals', tBeginSec);
for i6h = I_6H
  Tint6h = Tint24h(1) + 6*60*60*(i6h+[0, 1]);
  modify_save_subinterval_plot(h, hCisText1, hCisText2, Data, Tint6h, OutputPaths.dir6h)
end

%===========================
% Iterate over 2h intervals
%===========================
tBeginSec = solo.qli.utils.log_time('Begin iterating over 2 h intervals', tBeginSec);
for i2h = I_2H
  Tint2h = Tint24h(1) + 2*60*60*(i2h+[0, 1]);
  modify_save_subinterval_plot(h, hCisText1, hCisText2, Data, Tint2h, OutputPaths.dir2h)
end



close(fig);

[~] = solo.qli.utils.log_time('End of generate_quicklooks_24h_6h_2h.m', tBeginSec);

end



% Function to remove duplicated code.
%
% Presumes pre-existing figure with specific axes. Uses customized code to zoom
% in on the sub-time interval and adjusts the y limits for that interval.
function modify_save_subinterval_plot(hAxesArray, hCisText1, hCisText2, Data, Tint, parentDirPath)
assert(isa(hCisText1, 'matlab.graphics.primitive.Text'))
assert(isa(hCisText2, 'matlab.graphics.primitive.Text'))
assert(isstruct(Data))
assert(isa(Tint,      'EpochTT') && (length(Tint) == 2))

irf_zoom(hAxesArray, 'x', Tint);

%irf_zoom(hAxesArray(1), 'y');
%adjust_panel_ylimits_N_B(  hAxesArray(2), data,      Tint)
%adjust_panel_ylimits_Ti(   hAxesArray(5), data.Tpas, Tint)
%adjust_panel_ylimits_VT_VN(hAxesArray(6), data.Vpas, Tint)
%adjust_panel_ylimits_ESRF( hAxesArray(8), data.E,    Tint)
%irf_zoom(hAxesArray(7), 'y');

% NOTE: Different from for 24h plots.
yyaxis(hAxesArray(2), 'left');
% set_YLim_YTick(hAxesArray([1:2, 5:9]), hAxesArray([]), hAxesArray([3:4, 10]))
solo.qli.utils.set_YLim_YTick_automatically( hAxesArray([1:2, 5:9]))
solo.qli.utils.ensure_axes_data_tick_margins(hAxesArray([3:4, 10 ]))



yyaxis(hAxesArray(2), 'right');
% set_YLim_YTick(hAxesArray([2]), hAxesArray([]), hAxesArray([]))
solo.qli.utils.set_YLim_YTick_automatically(hAxesArray([2]))

% Update text
[hCisText1.String, hCisText2.String] = solo.qli.utils.get_context_info_strings(Data.soloPos, Data.earthPos, Tint);

solo.qli.utils.save_figure_to_file(parentDirPath, Tint)
end



% Function to remove duplicated code.
% function adjust_panel_ylimits_N_B(hAxes, data, Tint)
%     assert(isa(hAxes, 'matlab.graphics.axis.Axes') && isscalar(hAxes))
%     assert(isstruct(data))
%     assert(isa(Tint,  'EpochTT'))
%
%     %=============
%     % Left Y axis
%     %=============
%     Neflag   = ~isempty(data.Ne)   && ~isempty(data.Ne.tlim(  Tint)) && ~all(isnan(data.Ne.tlim(Tint).data));
%     Npasflag = ~isempty(data.Npas) && ~isempty(data.Npas.tlim(Tint));
%     if Neflag && Npasflag
%         yyaxis(hAxes, 'left');
%         hAxes.YLim=[...
%             min(floor([min(data.Npas.tlim(Tint).data),min(data.Ne.tlim(Tint).data)])),...
%             max(ceil( [max(data.Npas.tlim(Tint).data),max(data.Ne.tlim(Tint).data)]))...
%         ];
%     elseif Neflag
%         yyaxis(hAxes,'left');
%         hAxes.YLim=[floor(min(data.Ne.tlim(Tint).data)),ceil(max(data.Ne.tlim(Tint).data))];
%     elseif Npasflag
%         yyaxis(hAxes,'left');
%         hAxes.YLim=[floor(min(data.Npas.tlim(Tint).data)),ceil(max(data.Npas.tlim(Tint).data))];
%     end
%
%     %==============
%     % Right Y axis
%     %==============
%     if ~isempty(data.B) && ~isempty(data.B.tlim(Tint)) && ~all(isnan(data.B.abs.tlim(Tint).data))
%         yyaxis(hAxes,'right');
%         hAxes.YLim=[floor(min(data.B.abs.tlim(Tint).data)),ceil(max(data.B.abs.tlim(Tint).data))];
%     end
% end



% Function to remove duplicated code.
% function adjust_panel_ylimits_Ti(hAxes, TpasTSeries, Tint)
%     Y_MARGIN = 2;
%
%     assert(isa(hAxes,       'matlab.graphics.axis.Axes') && isscalar(hAxes))
%     assert(isa(TpasTSeries, 'TSeries') || isempty(TpasTSeries))
%     assert(isa(Tint,        'EpochTT'))
%
%     if ~isempty(TpasTSeries)
%         minTi = min(rmmissing(TpasTSeries.tlim(Tint).abs.data));
%         maxTi = max(rmmissing(TpasTSeries.tlim(Tint).abs.data));
%
%         if ~isempty(minTi) && ~isempty(maxTi)
%             % Only zoom if min & max are not NaN (==> Avoid crash).
%             irf_zoom(hAxes,'y',[minTi-Y_MARGIN, maxTi+Y_MARGIN]);
%         end
%     end
% end



% Function to remove duplicated code.
% function adjust_panel_ylimits_VT_VN(hAxes, VpasTSeries, Tint)
%     Y_MARGIN = 10;
%
%     assert(isa(hAxes,       'matlab.graphics.axis.Axes') && isscalar(hAxes))
%     assert(isa(VpasTSeries, 'TSeries') || isempty(VpasTSeries))
%     assert(isa(Tint,        'EpochTT'))
%
%     if ~isempty(VpasTSeries)
%         minVy = min(rmmissing(VpasTSeries.y.tlim(Tint).data));
%         minVz = min(rmmissing(VpasTSeries.z.tlim(Tint).data));
%         maxVy = max(rmmissing(VpasTSeries.y.tlim(Tint).data));
%         maxVz = max(rmmissing(VpasTSeries.z.tlim(Tint).data));
%         maxV = max(maxVy,maxVz);
%         minV = min(minVy,minVz);
%         if ~isempty(minV) && ~isempty(maxV)
%             % Only zoom if min & max are not NaN (==> Avoid crash).
%             irf_zoom(hAxes,'y',[minV-Y_MARGIN, maxV+Y_MARGIN]);
%         end
%     end
% end



% Function to remove duplicated code.
% function adjust_panel_ylimits_ESRF(hAxes, ETSeries, Tint)
%     Y_MARGIN = 5;
%
%     assert(isa(hAxes,    'matlab.graphics.axis.Axes') && isscalar(hAxes))
%     assert(isa(ETSeries, 'TSeries') || isempty(ETSeries))
%     assert(isa(Tint,     'EpochTT'))
%
%     if ~isempty(ETSeries)
%         minEy = min(rmmissing(ETSeries.y.tlim(Tint).data));
%         maxEy = max(rmmissing(ETSeries.y.tlim(Tint).data));
%
%         if ~isempty(minEy) && ~isempty(maxEy)
%             irf_zoom(hAxes,'y',[minEy-Y_MARGIN, maxEy+Y_MARGIN]);
%         end
%     end
% end
