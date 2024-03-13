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
%     See the call from solo.qli.generate_quicklooks_all_types().
% OutputPaths
%     Struct with paths to separate output directories for the different types
%     of quicklooks (see solo.qli.generate_quicklooks_all_types).
% Tint24h
%     Should be a 24-hour time interval consistent with the time series in
%     "data", e.g.
%     irf.tint('2020-06-01T00:00:00.00Z','2020-06-02T00:00:00.00Z');
% logoPath
%     Either path to IRF logo, or empty.
%
%
% NOTES
% =====
% * Computes the spectrum for B when magnetic field data is available. When it
%   does, the code takes a lot of time.
% * The function uses solo.read_TNR() which in turns relies on a
%   hardcoded path to "/data/solo/remote/data/L2/thr/" and selected
%   subdirectories.
% * The function obtains some data by reading CDF files directly (cdfread;
%   solo_L2_swa-pas-eflux).



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
% TODO-NI: Panels 2 & 5 are logarithmic for 24h plots and linear for 6h & 2h?
%          Should they be?
% TODO-NI: Old panels 5 had a constant y axis range (YLim) for 24h, but dynamic
%          (changing) for 6h & 2h. Was that intentional?
% TODO-NI Panel 10 (log) is hardcoded to YLim~[10, 100] (because that is what
%         it used to be). This does not cover the entire interval of data
%         (there is more data at lower y). Should it be that way?
%
% PROPOSAL: Make function not directly call solo.read_TNR()
%   PRO: Makes function testable.
%     CON: Function still reads other CDF files.
%   CON: Must understand the solo.read_TNR() return value.
%     CON: Seems feasible.
%       case 0:
%         out = 0;
%       case 1:
%         out = struct('t', time_.epochUnix, 'f', freq_tnr, 'p',vp.^10);
%         out.p_label={'dB'};
%
%   PROPOSAL: Only call solo.read_TNR() via dependency injection.
%     CON: Overkill.
%   PROPOSAL: Submit the return value of solo.read_TNR() as argument instead of
%             calling it.



tBeginSec = tic();



% Setup figure:
LWIDTH   = 1.0;   % irf_plot() line width.
FSIZE    = 18;    % Font size
LEG_SIZE = 22;    % irf_legend() font size.
COLORS   = [0 0 0; 0 0 1; 1 0 0; 0 0.5 0; 0 1 1; 1 0 1; 1 1 0];

Units = irf_units;
Me    = Units.me;      % Electron mass [kg]
epso  = Units.eps0;    % Permitivitty of free space [Fm^-1]
mp    = Units.mp;      % Proton mass [km]
qe    = Units.e;       % Elementary charge [C]

h            = irf_plot(10,'newfigure');
fig          = gcf;
fig.Position = [1,1,1095,800];



%===================================
% Fill panel 1: B vector components
%===================================
if ~isempty(Data.B)
  irf_plot(h(1),Data.B.tlim(Tint24h),'linewidth',LWIDTH);
  hold(    h(1),'on');
  irf_plot(h(1),Data.B.abs.tlim(Tint24h),'linewidth',LWIDTH);
end
irf_legend(h(1),{'B_{R}','B_{T}','B_{N}','|B|'},[0.98 0.18],'Fontsize',LEG_SIZE);
ylabel(    h(1),{'B_{RTN}';'(nT)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 1', tBeginSec);



%=======================
% Fill panel 2: N & |B|
%=======================
hold(h(2),'on');
if ~isempty(Data.Ne)
  irf_plot(h(2),Data.Ne.tlim(Tint24h),'-','color',COLORS(1,:),'linewidth',LWIDTH);
end
if ~isempty(Data.Npas)
  irf_plot(h(2),Data.Npas.tlim(Tint24h),'-','color',COLORS(2,:),'linewidth',LWIDTH);
end
ylabel(h(2),{'N';'(cm^{-3})'},'interpreter','tex','fontsize',FSIZE);
h(2).ColorOrder=COLORS;
irf_legend(h(2),{'N_{e,RPW}','N_{i,PAS}','|B|'},[0.98 0.16],'Fontsize',LEG_SIZE);

yyaxis(h(2),'right');
if ~isempty(Data.B)
  fci = qe*Data.B.abs*10^-9/mp/(2*pi);
  irf_plot(h(2),Data.B.abs.tlim(Tint24h),'color',COLORS(3,:),'linewidth',LWIDTH);
  %Bnan = rmmissing(data.B.abs.data);
  %if ~isempty(Bnan)
  %    h(2).YLim=[floor(min(abs(Bnan))),ceil(max(abs(Bnan)))];
  %end
  minAbsB = min(Data.B.tlim(Tint24h).abs.data);
  maxAbsB = max(Data.B.tlim(Tint24h).abs.data);
  if ~isnan(minAbsB) && ~isnan(maxAbsB)
    % Only zoom if min & max are not NaN (==> Avoid crash).
    irf_zoom(h(2),'y',[minAbsB-1, maxAbsB+1]);
  end
end
ylabel(h(2),{'|B|';'(nT)'},'interpreter','tex','fontsize',FSIZE);
h(2).YColor=[1,0,0];

tBeginSec = solo.qli.utils.log_time('End panel 2', tBeginSec);



%===========================
% Fill panel 3 & 4: Spectra
%===========================
if ~isempty(Data.B) && solo.qli.const.NONWEEKLY_SPECTRA_ENABLED
  if  ~isempty(rmmissing(Data.B.data))
    bb = Data.B;
    if median(diff((bb.time.epochUnix))) < 0.1250*0.95
      fMag = 128; fMax = 7;
    else
      fMag =   8; fMax = 3;
    end
    b0 = bb.filt(0, 0.01,fMag, 5);

    % IMPORTANT NOTE: The call to irf_ebsp() is very time-consuming.
    tBeginSec = solo.qli.utils.log_time('irf_ebsp(): Begin call', tBeginSec);
    ebsp = irf_ebsp([],bb,[],b0,[],[0.05 fMax],'fullB=dB', 'polarization', 'fac');
    tBeginSec = solo.qli.utils.log_time('irf_ebsp(): End call', tBeginSec);

    frequency   = ebsp.f;
    time        = ebsp.t;
    Bsum        = ebsp.bb_xxyyzzss(:,:,4);
    ellipticity = ebsp.ellipticity;
    dop         = ebsp.dop;

    % Remove points with very low degree of polarization
    dopthresh = 0.7;
    removepts = find(dop < dopthresh);
    ellipticity(removepts) = NaN;

    % Remove "lonely" pixels
    msk              = ellipticity;
    msk(~isnan(msk)) = 1;
    msk( isnan(msk)) = 0;
    % "BW2 = bwareaopen(BW,P) removes from a binary image all connected
    % components (objects) that have fewer than P pixels, producing another
    % binary image BW2."
    msk_denoise                 = bwareaopen(msk,8);
    ellipticity(msk_denoise==0) = NaN;

    %---------
    % Panel 3
    %---------
    specrec         = struct('t',time);
    specrec.f       = frequency;
    specrec.p       = Bsum;
    specrec.f_label = '';
    specrec.p_label = {'log_{10}B^{2}', 'nT^2 Hz^{-1}'};
    irf_spectrogram(h(3), specrec, 'log', 'donotfitcolorbarlabel');
    set(     h(3), 'yscale', 'log');
    % set(h(1),'ytick',[1e1 1e2 1e3]);
    % caxis(h(3),[-8 -1])
    hold(    h(3), 'on');
    irf_plot(h(3), fci, 'k', 'linewidth', LWIDTH);
    text(    h(3), 0.01, 0.3, 'f_{ci}', 'units', 'normalized', 'fontsize', FSIZE);
    colormap(h(3), 'jet');

    %---------
    % Panel 4
    %---------
    specrec         = struct('t',time);
    specrec.f       = frequency;
    specrec.p       = ellipticity;
    specrec.f_label = '';
    specrec.p_label = {'Ellipticity', 'DOP>0.7'};
    irf_spectrogram(h(4), specrec, 'log', 'donotfitcolorbarlabel');
    set(     h(4), 'yscale', 'log');
    % set(h(1),'ytick',[1e1 1e2 1e3]);
    caxis(   h(4), [-1 1])
    hold(    h(4), 'on');
    irf_plot(h(4), fci, 'k', 'linewidth', LWIDTH);
    text(    h(4), 0.01, 0.3, 'f_{ci}', 'units', 'normalized', 'fontsize', FSIZE);

    crr     = interp1([1 64 128 192 256], [0.0  0.5 0.75 1.0 0.75], 1:256);
    cgg     = interp1([1 64 128 192 256], [0.0  0.5 0.75 0.5 0.00], 1:256);
    cbb     = interp1([1 64 128 192 256], [0.75 1.0 0.75 0.5 0.00], 1:256);
    bgrcmap = [crr' cgg' cbb'];
    colormap(h(4), bgrcmap);
  end
end
ylabel(h(3), {'f';'(Hz)'}, 'fontsize', FSIZE);
ylabel(h(4), {'f';'(Hz)'}, 'fontsize', FSIZE);
tBeginSec = solo.qli.utils.log_time('End panel 3 & 4', tBeginSec);



%===============================
% Fill panel 5: Ion temperature
%===============================
if ~isempty(Data.Tpas)
  irf_plot(h(5),Data.Tpas.tlim(Tint24h),'color',COLORS(2,:),'linewidth',LWIDTH);
end
irf_zoom(h(5),'y');
ylabel(  h(5),{'T_i';'(eV)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 5', tBeginSec);



%==================================
% Fill panel 6: y,z PAS velocities
%==================================
if ~isempty(Data.Vpas)
  irf_plot(h(6),Data.Vpas.y.tlim(Tint24h),'color',COLORS(2,:),'linewidth',LWIDTH);
  hold(    h(6),'on');
  irf_plot(h(6),Data.Vpas.z.tlim(Tint24h),'color',COLORS(3,:),'linewidth',LWIDTH);
end
irf_legend(h(6),{'','v_{T}','v_{N}'},[0.98 0.18],'Fontsize',LEG_SIZE);
irf_zoom(  h(6),'y');
ylabel(    h(6),{'V_{T,N}';'(km/s)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 6', tBeginSec);



%=====================================
% Fill panel 7: Vrpw, Vpas velocities
%=====================================
hold(h(7),'on');
if ~isempty(Data.Vrpw)
  irf_plot(h(7),-Data.Vrpw,'o-','color',COLORS(1,:));
end
if ~isempty(Data.Vpas)
  irf_plot(h(7),Data.Vpas.x.tlim(Tint24h),'color',COLORS(2,:),'linewidth',LWIDTH);
end
irf_legend(h(7), {'V_{RPW}','V_{PAS}'},[0.98 0.18],'Fontsize',LEG_SIZE);
irf_zoom(  h(7), 'y');
ylabel(    h(7), {'V_R';'(km/s)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 7', tBeginSec);



%==============================
% Fill panel 8: Electric field
%==============================
if ~isempty(Data.E)
  irf_plot(h(8),Data.E.y,'color',COLORS(2,:),'linewidth',LWIDTH)
  hold(h(8),'on');
  %irf_plot(h(8),data.E.z,'color',COLORS(3,:),'linewidth',LWIDTH)

  minEy = min(rmmissing(Data.E.y.data));
  maxEy = max(rmmissing(Data.E.y.data));
  if ~isempty(minEy) && ~isempty(maxEy)
    irf_zoom(h(8),'y',[minEy-5 maxEy+5]);
  end
end
irf_legend(h(8),{'','E_y'},[0.98 0.20],'Fontsize',LEG_SIZE);
ylabel(h(8),{'E_{SRF}';'(mV/m)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 8', tBeginSec);



%===================================
% Fill panel 9: Ion energy spectrum
% ---------------------------------
% NOTE: READS CDF FILES!
%===================================
if ~isempty(Data.ieflux)
  SwaFileArray = solo.db_list_files('solo_L2_swa-pas-eflux',Tint24h);
  iDEF         = struct('t',  Data.ieflux.tlim(Tint24h).time.epochUnix);
  % for ii = 1:round((myFile(end).stop-myFile(1).start)/3600/24)
  for iFile = 1:length(SwaFileArray)
    iEnergy = cdfread(...
      fullfile(SwaFileArray(iFile).path, SwaFileArray(iFile).name), ...
      'variables', 'Energy');
    iEnergy = iEnergy{1};
    iDEF.p  = Data.ieflux.data;
  end
  iDEF.p_label = {'dEF','keV/','(cm^2 s sr keV)'};
  iDEF.f       = repmat(iEnergy,1,numel(iDEF.t))';
  irf_spectrogram(h(9),iDEF,'log','donotfitcolorbarlabel');
  % set(h(1),'ytick',[1e1 1e2 1e3]);
  %caxis(h(9),[-1 1])
  hold(h(9),'on');
  h9_clims = h(9).CLim;
  % Fix color axis
  h9_medp = mean(iDEF.p);
  h9_medp = min(h9_medp(h9_medp>0));
  if h9_medp > 0 && h9_medp > h9_clims(1) && log10(h9_medp)+2<(max(max(log10(iDEF.p))))
    caxis(h(9),[log10(h9_medp)+2 (max(max(log10(iDEF.p))))])
  end
end
set(     h(9), 'YScale', 'log');
colormap(h(9), jet)
ylabel(  h(9), {'W_{i}';'(eV)'},'interpreter','tex','fontsize',FSIZE);
tBeginSec = solo.qli.utils.log_time('End panel 9', tBeginSec);



%=======================================
% Fill panel 10: E-field spectrum (TNR)
%=======================================
% BUG(?): PROBABLY WHAT HAPPENS: Does not create color bar (no call to
% colormap()) when no data.
% ==> The panel becomes wider.
% ==> Other panels become wider.
% ==> Moves the IRF logo to the right, and partially outside image.
if ~isempty(Data.Etnr)    % && false
  try
    [TNR] = solo.read_TNR(Tint24h);
  catch Exc
    if strcmp(Exc.identifier, 'read_TNR:FileNotFound')
      TNR = [];
    end
  end
  if isa(TNR,'struct')
    sz_tnr = size(TNR.p);
    if sz_tnr(1) == length(TNR.t) && sz_tnr(2) == length(TNR.f)
      irf_spectrogram(h(10), TNR, 'log', 'donotfitcolorbarlabel')
      hold(           h(10), 'on');
      if ~isempty(Data.Ne)
        % Electron plasma frequency
        wpe_sc = (sqrt(((Data.Ne.tlim(Tint24h)*1000000)*qe^2)/(Me*epso)));
        fpe_sc = (wpe_sc/2/pi)/1000;
        irf_plot(h(10),fpe_sc,'r','linewidth',LWIDTH);
        fpe_sc.units = 'kHz';
        fpe_sc.name  = 'f [kHz]';
      end
      hold(h(10), 'off');
      text(h(10), 0.01, 0.3, 'f_{pe,RPW}', 'units', 'normalized', 'fontsize', FSIZE, 'Color', 'r');
      set( h(10), 'YScale', 'log');
      %set(h(10),'ColorScale','log')
      %caxis(h(10),[.01 1]*10^-12)
      ylabel(h(10),{'f';'(kHz)'},'interpreter','tex','fontsize',FSIZE);
      colormap(h(10),jet)
      %yticks(h(10),[10^1 10^2]);
      %irf_zoom(h(10),'y',[10^1 10^2])
    end
  end
end
ylabel(h(10), {'f';'(kHz)'},'interpreter','tex','fontsize',FSIZE);
yticks(h(10),      [10^1 10^2]);
% Not set YLim.
%irf_zoom(h(10),'y',[10^1 10^2])
irf_zoom(h(10),'y',[9, 110])    % Not overwritten later.

if isempty(Data.Vrpw) ...
    && isempty(Data.E)    && isempty(Data.Ne)   && isempty(Data.B) ...
    && isempty(Data.Tpas) && isempty(Data.Npas) && isempty(Data.ieflux) ...
    && isempty(Data.Etnr)
  nanPlot = irf.ts_scalar(Tint24h,ones(1,2)*NaN);
  irf_plot(h(10),nanPlot);    % No LWIDTH?
  grid(h(10),'off');
  ylabel(h(10),{'f';'(kHz)'},'interpreter','tex','fontsize',FSIZE);
end

tBeginSec = solo.qli.utils.log_time('End panel 10', tBeginSec);



%======================
% Other, miscellaneous
%======================
irf_plot_axis_align(h(1:10));  % Make panels ("data area") have the same length.
irf_zoom(h(1:10),'x',Tint24h);
irf_zoom(h(1),'y');

h(2).YLabel.Position=[1.05,0.5,0];
yyaxis(h(2),'left');
h(2).YLabel.Units='normalized';
% NOTE: *NOT* using h(3) or h(4) since they do not have ylabels if the relevant
% data is missing. h(1) ylabel always has a position. Using h(3) (old
% implementation) lead to left panel 2 ylabel having the wrong position (too far
% left) when h(3) did not have any label).
h(2).YLabel.Position=h(1).YLabel.Position;

% Add Context Info Strings (CIS): Spacecraft position, Earth longitude as text.
[soloStr, earthStr] = solo.qli.utils.get_context_info_strings(Data.soloPos, Data.earthPos, Tint24h);
hCisText1 = text(h(10), -0.11, -0.575, soloStr,  'units', 'normalized', 'fontsize', FSIZE);
hCisText2 = text(h(10), -0.11, -0.925, earthStr, 'units', 'normalized', 'fontsize', FSIZE);



xtickangle(h(10), 0)

%======================================================
% Add IRF logo and data source information info string
%======================================================
logoPos = h(1).Position;    %  [left, bottom, width, height]
logoPos(1) = logoPos(1) + logoPos(3) + 0.06;
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
yyaxis(h(2), 'left');
h(2).YScale = 'log';       % NOTE: Later changed to LIN for non-24h.
h(2).YTick  = [1, 10, 100];

% NOTE: Panel 2 YTick not auto-adjusted partly because
% solo.qli.utils.ensure_axes_data_tick_margins() can not handle both left &
% right yaxis.
yyaxis(h(2), 'right');
h(2).YScale = 'log';       % NOTE: Later changed to LIN for non-24h.
h(2).YTick  = [1, 10, 100];

% NOTE: h(5).YLim are hardcoded and seem too broad/wide.
% PROPOSAL: Not overwrite automatic YLim?
oldlims5  = h(5).YLim;
oldticks5 = h(5).YTick;
h(5).YScale = 'log';       % NOTE: Later changed to LIN.
h(5).YTick  = [1, 10, 100];
h(5).YLim   = [0.5, 300];



% Remove overlapping ticks.
%solo.qli.utils.ensure_axes_data_tick_margins(h)
% Automatically set YLim+YTick, or automatically set YLim, or adjust YLim,
% depending on panel.
yyaxis(h(2), 'right');
set_YLim_YTick(h([1, 3:4, 6:9]), h([2]), h([5, 10]))
yyaxis(h(2), 'left');
set_YLim_YTick(h([]), h([2]), h([]))



fig = gcf;
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
solo.qli.utils.save_figure_to_file(OutputPaths.path_24h, Tint24h)



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
h(5).YLim   = oldlims5;
h(5).YTick  = oldticks5;



I_6H = 0:3;
I_2H = 0:11;
if ~solo.qli.const.NONWEEKLY_ALL_PLOTS_ENABLED
  % For debugging/testing.
  %I_6H = [0];
  %I_2H = [0];
  I_6H = [1];
  I_2H = [5];
end

%===========================
% Iterate over 6h intervals
%===========================
tBeginSec = solo.qli.utils.log_time('Begin iterating over 6 h intervals', tBeginSec);
for i6h = I_6H
  Tint6h = Tint24h(1) + 6*60*60*(i6h+[0, 1]);
  modify_save_subinterval_plot(h, hCisText1, hCisText2, Data, Tint6h, OutputPaths.path_6h)
end

%===========================
% Iterate over 2h intervals
%===========================
tBeginSec = solo.qli.utils.log_time('Begin iterating over 2 h intervals', tBeginSec);
for i2h = I_2H
  Tint2h = Tint24h(1) + 2*60*60*(i2h+[0, 1]);
  modify_save_subinterval_plot(h, hCisText1, hCisText2, Data, Tint2h, OutputPaths.path_2h)
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
yyaxis(hAxesArray(2), 'right');
set_YLim_YTick(hAxesArray(1:9), hAxesArray([]), hAxesArray([10]))
yyaxis(hAxesArray(2), 'left');
set_YLim_YTick(hAxesArray([2]), hAxesArray([]), hAxesArray([]))

% Update text
[hCisText1.String, hCisText2.String] = solo.qli.utils.get_context_info_strings(Data.soloPos, Data.earthPos, Tint);

solo.qli.utils.save_figure_to_file(parentDirPath, Tint)
end



% Automatically set or adjust y limits and y tick positions.
% "Always" ensures that ticks are not at the min/max (YLim) to avoid overlapping labels.
% Can specify different behaviour for different axes.
%
% NOTE: Function can not simultaneously handle both yyaxis left & right.
% NOTE: MATLAB's automatic setting of y ticks for log scale (and which is used)
%       can be bad. May therefore want to set y ticks for panels with log scale.
%
% ARGUMENTS
% =========
% hAxesAutoYLimYTickArray
%       Array of axes for which to
%       (1) set YLim automatically (from data; with margins)
%       (2) set YTick automatically.
% hAxesAutoYLimArray
%       Array of axes for which to
%       (1) set YLim automatically (from data; with margins),
%       (2) keep YTick as is.
% hAxesMarginYLimArray
%       Array of axes for which to
%       (1) add margins to pre-existing YLim
%       (2) keep YTick as is.
%
function set_YLim_YTick(hAxesAutoYLimYTickArray, hAxesAutoYLimArray, hAxesMarginYLimArray)
% PROPOSAL: Automatically (not MATLAB) set YTick for logarithmic axis to
%           ensure one tick per power of ten, 10^n.
% PROPOSAL: Set YLimMode=manual for YLimYTick axes.
% PROPOSAL: Split into 2/3 separate functions.
%     NOTE: solo.qli.utils.ensure_axes_data_tick_margins)() is called for all
%           axes.

assert(isa(hAxesAutoYLimYTickArray, 'matlab.graphics.axis.Axes'))
assert(isa(hAxesAutoYLimArray,      'matlab.graphics.axis.Axes'))
assert(isa(hAxesMarginYLimArray,    'matlab.graphics.axis.Axes'))

hAxesArray = unique([...
  hAxesAutoYLimYTickArray(:); ...
  hAxesAutoYLimArray(:); ...
  hAxesMarginYLimArray(:) ...
  ]);

% Assert that all axes are unique (no overlap/intersection).
assert(numel(hAxesArray) == (...
  numel(hAxesAutoYLimYTickArray) + ...
  numel(hAxesAutoYLimArray) + ...
  numel(hAxesMarginYLimArray)...
  ))

%=======================================================================
% Automatically set preliminary YLim (y limits) and final YTick (y tick
% positions) for selected axes.
%=======================================================================
% Set axes y range (YLim) to only cover the data (plus rounding outwards
% to ticks).
set(hAxesAutoYLimYTickArray, 'YLimMode', 'auto')
% Auto-generate ticks (YTick; y values at which there should be ticks).
set(hAxesAutoYLimYTickArray, 'YTickMode', 'auto')
%---------------------------------------------------------------------------
% IMPORTANT: Read YLim without using the return result ("do nothing")
% --------------------------------------------------------------------
% IMPLEMENTATION NOTE: THIS COMMAND SHOULD THEORETICALLY NOT BE NEEDED,
% BUT IS NEEDED FOR THE YLim VALUES TO BE SET PROPERLY. MATLAB BUG?!
% This behaviour has been observed on Erik P G Johansson's laptop
% "irony" (MATLAB R2019b, Ubuntu Linux) as of 2023-05-25.
% Ex: (Re-)scaling of panel 5, 2022-02-23T10-12 (2h plot).
get(hAxesAutoYLimYTickArray, 'YLim');
%---------------------------------------------------------------------------
% Prevent later setting of YLim (next command) from generating new ticks.
set(hAxesAutoYLimYTickArray, 'YTickMode', 'manual')

%=========================================================================
% Automatically set YLim (y limits) but keep old YTick (y tick positions)
% for selected axes.
%=========================================================================
set(hAxesAutoYLimArray, 'YTickMode', 'manual')
%get(hAxesManualArray, 'YLim');   % READ ONLY. UNNECESSARY?
set(hAxesAutoYLimArray, 'YLimMode',  'auto')
%get(hAxesManualArray, 'YLim');   % READ ONLY. UNNECESSARY?

%===========================================================================
% If needed, adjust YLim (but not YTick) to ensure there are margins between
% YLim and YTick = No ticks on the panel edges/corners.
%===========================================================================
%i = 10;
%h = hAxesArray(i);
%fprintf('hAxesArray(%i).YLim  = %s\n', i, num2str(h.YLim))
%fprintf('hAxesArray(%i).YTick = %s\n', i, num2str(h.YTick))
solo.qli.utils.ensure_axes_data_tick_margins(hAxesArray)
%fprintf('hAxesArray(%i).YLim  = %s\n', i, num2str(h.YLim))
%fprintf('hAxesArray(%i).YTick = %s\n', i, num2str(h.YTick))
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
%         yyaxis(hAxes,'left');
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