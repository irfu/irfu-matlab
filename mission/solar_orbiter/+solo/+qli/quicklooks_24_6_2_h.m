function quicklooks_24_6_2_h(data,paths,Tint_24h,logoPath)
% Given data in the struct 'data' (see solo.qli.quicklooks_main), generates
% plots and saves them in the paths specified in the struct 'paths' (see
% solo.qli.quicklooks_main). Computes spectrum of B, so takes a while to run.
% Tint_24h should be a 24hour time interval, e.g.
% irf.tint('2020-06-01T00:00:00.00Z','2020-06-02T00:00:00.00Z');

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
% TODO-DEC: Panels 2 & 5 are logarithmic for 24h plots and linear for 6h & 2h?
%           Should they be?
% TODO-DEC: Panel 10 (log) used to be hardcoded to YLim~[10, 100] which does not
%           cover the entire interval of data (there is ore at lower y). Should
%           it be?


tBeginSec = tic();



% Whether to enable/disable panels with time-consuming spetra. Disabling these
% is useful for debugging and testing. Should be enabled by default.
SPECTRA_ENABLED = 1;
% Whether to generate all plots or only some (e.g. one 6h plot, one 2h plot).
% Disabling this is useful for debugging and testing. Should be enabled by
% default.
ALL_PLOTS_ENABLED = 1;



% Setup figure:
LWIDTH   = 1.0;   % irf_plot() Line width
FSIZE    = 18;    % Font size
LEG_SIZE = 22;    % irf_legend() font size
COLORS   = [0 0 0;0 0 1;1 0 0;0 0.5 0;0 1 1 ;1 0 1; 1 1 0];

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
if ~isempty(data.B)
    irf_plot(h(1),data.B.tlim(Tint_24h),'linewidth',LWIDTH);
    hold(h(1),'on');
    irf_plot(h(1),data.B.abs.tlim(Tint_24h),'linewidth',LWIDTH);
end
irf_legend(h(1),{'B_{R}','B_{T}','B_{N}','|B|'},[0.98 0.18],'Fontsize',LEG_SIZE);
ylabel(h(1),{'B_{RTN}';'(nT)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 1', tBeginSec);



%=======================
% Fill panel 2: N & |B|
%=======================
hold(h(2),'on');
if ~isempty(data.Ne)
    irf_plot(h(2),data.Ne.tlim(Tint_24h),'-','color',COLORS(1,:),'linewidth',LWIDTH);
end
if ~isempty(data.Npas)
    irf_plot(h(2),data.Npas.tlim(Tint_24h),'-','color',COLORS(2,:),'linewidth',LWIDTH);
end
ylabel(h(2),{'N';'(cm^{-3})'},'interpreter','tex','fontsize',FSIZE);
h(2).ColorOrder=COLORS;
irf_legend(h(2),{'N_{e,RPW}','N_{i,PAS}','|B|'},[0.98 0.16],'Fontsize',LEG_SIZE);

yyaxis(h(2),'right');
if ~isempty(data.B)
    fci = qe*data.B.abs*10^-9/mp/(2*pi);
    irf_plot(h(2),data.B.abs.tlim(Tint_24h),'color',COLORS(3,:),'linewidth',LWIDTH);
    %Bnan = rmmissing(data.B.abs.data);
    %if ~isempty(Bnan)
    %    h(2).YLim=[floor(min(abs(Bnan))),ceil(max(abs(Bnan)))];
    %end
    minAbsB = min(data.B.tlim(Tint_24h).abs.data);
    maxAbsB = max(data.B.tlim(Tint_24h).abs.data);
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
if ~isempty(data.B) && SPECTRA_ENABLED
    if  ~isempty(rmmissing(data.B.data))
        bb = data.B;
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

        frequency = ebsp.f;
        time = ebsp.t;
        Bsum = ebsp.bb_xxyyzzss(:,:,4);
        ellipticity = ebsp.ellipticity;
        dop = ebsp.dop;

        % Remove points with very low degree of polarization
        dopthresh = 0.7;
        removepts = find(dop < dopthresh);
        ellipticity(removepts) = NaN;

        % Remove "lonely" pixels
        msk = ellipticity;
        msk(~isnan(msk)) = 1;
        msk(isnan(msk)) = 0;
        msk_denoise = bwareaopen(msk,8);
        ellipticity(msk_denoise==0) = NaN;

        % Plot
        specrec=struct('t',time);
        specrec.f=frequency;
        specrec.p=Bsum;
        specrec.f_label='';
        specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
        irf_spectrogram(h(3),specrec,'log','donotfitcolorbarlabel');
        set(h(3),'yscale','log');
        % set(h(1),'ytick',[1e1 1e2 1e3]);
        % caxis(h(3),[-8 -1])
        hold(h(3),'on');
        irf_plot(h(3),fci,'k','linewidth',LWIDTH);
        text(h(3),0.01,0.3,'f_{ci}','units','normalized','fontsize',18);
        colormap(h(3),'jet');

        specrec=struct('t',time);
        specrec.f=frequency;
        specrec.p=ellipticity;
        specrec.f_label='';
        specrec.p_label={'Ellipticity','DOP>0.7'};
        irf_spectrogram(h(4),specrec,'log','donotfitcolorbarlabel');
        set(h(4),'yscale','log');
        % set(h(1),'ytick',[1e1 1e2 1e3]);
        caxis(h(4),[-1 1])
        hold(h(4),'on');
        irf_plot(h(4),fci,'k','linewidth',LWIDTH);
        text(h(4),0.01,0.3,'f_{ci}','units','normalized','fontsize',18);

        crr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
        cgg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
        cbb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
        bgrcmap = [crr' cgg' cbb'];
        colormap(h(4),bgrcmap);
    end
end
ylabel(h(3),{'f';'(Hz)'},'fontsize',FSIZE);
ylabel(h(4),{'f';'(Hz)'},'fontsize',FSIZE);
tBeginSec = solo.qli.utils.log_time('End panel 3 & 4', tBeginSec);



%===============================
% Fill panel 5: Ion temperature
%===============================
if ~isempty(data.Tpas)
    irf_plot(h(5),data.Tpas.tlim(Tint_24h),'color',COLORS(2,:),'linewidth',LWIDTH);
end
irf_zoom(h(5),'y');
ylabel(h(5),{'T_i';'(eV)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 5', tBeginSec);



%==================================
% Fill panel 6: y,z PAS velocities
%==================================
if ~isempty(data.Vpas)
    irf_plot(h(6),data.Vpas.y.tlim(Tint_24h),'color',COLORS(2,:),'linewidth',LWIDTH);
    hold(h(6),'on');
    irf_plot(h(6),data.Vpas.z.tlim(Tint_24h),'color',COLORS(3,:),'linewidth',LWIDTH);
end
irf_legend(h(6),{'','v_{T}','v_{N}'},[0.98 0.18],'Fontsize',LEG_SIZE);
irf_zoom(h(6),'y');
ylabel(h(6),{'V_{T,N}';'(km/s)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 6', tBeginSec);



%=====================================
% Fill panel 7: Vrpw, Vpas velocities
%=====================================
hold(h(7),'on');
if ~isempty(data.Vrpw)
    irf_plot(h(7),-data.Vrpw,'o-','color',COLORS(1,:));
end
if ~isempty(data.Vpas)
    irf_plot(h(7),data.Vpas.x.tlim(Tint_24h),'color',COLORS(2,:),'linewidth',LWIDTH);
end
irf_legend(h(7),{'V_{RPW}','V_{PAS}'},[0.98 0.18],'Fontsize',LEG_SIZE);
irf_zoom(h(7),'y');
ylabel(h(7),{'V_R';'(km/s)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 7', tBeginSec);



%==============================
% Fill panel 8: Electric field
%==============================
if ~isempty(data.E)
    irf_plot(h(8),data.E.y,'color',COLORS(2,:),'linewidth',LWIDTH)
    hold(h(8),'on');
    %irf_plot(h(8),data.E.z,'color',COLORS(3,:),'linewidth',LWIDTH)

    minEy = min(rmmissing(data.E.y.data));
    maxEy = max(rmmissing(data.E.y.data));
    if ~isempty(minEy) && ~isempty(maxEy)
        irf_zoom(h(8),'y',[minEy-5 maxEy+5]);
    end
end
irf_legend(h(8),{'','E_y'},[0.98 0.20],'Fontsize',LEG_SIZE);
ylabel(h(8),{'E_{SRF}';'(mV/m)'},'interpreter','tex','fontsize',FSIZE);

tBeginSec = solo.qli.utils.log_time('End panel 8', tBeginSec);



%===================================
% Fill panel 9: Ion energy spectrum
%===================================
if ~isempty(data.ieflux)
    myFile=solo.db_list_files('solo_L2_swa-pas-eflux',Tint_24h);
    iDEF   = struct('t',  data.ieflux.tlim(Tint_24h).time.epochUnix);
    % for ii = 1:round((myFile(end).stop-myFile(1).start)/3600/24)
    for ii = 1:length(myFile)
        iEnergy = cdfread([myFile(ii).path '/' myFile(ii).name],'variables','Energy');
        iEnergy = iEnergy{1};
        iDEF.p = data.ieflux.data;
    end
    iDEF.p_label={'dEF','keV/','(cm^2 s sr keV)'};
    iDEF.f = repmat(iEnergy,1,numel(iDEF.t))';
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
set(h(9), 'YScale', 'log');
colormap(h(9),jet)
ylabel(h(9),{'W_{i}';'(eV)'},'interpreter','tex','fontsize',FSIZE);
tBeginSec = solo.qli.utils.log_time('End panel 9', tBeginSec);



%=======================================
% Fill panel 10: E-field spectrum (TNR)
%=======================================
if ~isempty(data.Etnr)    % && false
    try
        [TNR] = solo.read_TNR(Tint_24h);
    catch Exc
        if strcmp(Exc.identifier, 'read_TNR:FileNotFound')
            TNR = [];
        end
    end
    if isa(TNR,'struct')
        sz_tnr = size(TNR.p);
        if sz_tnr(1) == length(TNR.t) && sz_tnr(2) == length(TNR.f)
            irf_spectrogram(h(10),TNR,'log','donotfitcolorbarlabel')
            hold(h(10),'on');
            if ~isempty(data.Ne)
                % Electron plasma frequency
                wpe_sc = (sqrt(((data.Ne.tlim(Tint_24h)*1000000)*qe^2)/(Me*epso)));
                fpe_sc = (wpe_sc/2/pi)/1000;
                irf_plot(h(10),fpe_sc,'r','linewidth',LWIDTH);
                fpe_sc.units = 'kHz';
                fpe_sc.name = 'f [kHz]';
            end
            hold(h(10),'off');
            text(h(10),0.01,0.3,'f_{pe,RPW}','units','normalized','fontsize',18,'Color','r');
            set(h(10), 'YScale', 'log');
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
irf_zoom(h(10),'y',[10^1 10^2])

if isempty(data.Vrpw) ...
        && isempty(data.E)    && isempty(data.Ne)   && isempty(data.B) ...
        && isempty(data.Tpas) && isempty(data.Npas) && isempty(data.ieflux) ...
        && isempty(data.Etnr)
    nanPlot = irf.ts_scalar(Tint_24h,ones(1,2)*NaN);
    irf_plot(h(10),nanPlot);    % No LWIDTH?
    grid(h(10),'off');
    ylabel(h(10),{'f';'(kHz)'},'interpreter','tex','fontsize',FSIZE);
end

tBeginSec = solo.qli.utils.log_time('End panel 10', tBeginSec);



%======================
% Other, miscellaneous
%======================
irf_plot_axis_align(h(1:10));  % Make panels ("data area") have the same length.
irf_zoom(h(1:10),'x',Tint_24h);
irf_zoom(h(1),'y');

h(2).YLabel.Position=[1.05,0.5,0];
yyaxis(h(2),'left');
h(2).YLabel.Units='normalized';
% NOTE: *NOT* using h(3) or h(4) since they do not have ylabels if the relevant
% data is missing. h(1) ylabel always has a position. Using h(3) (old
% implementation) lead to left panel 2 ylabel having the wrong position (too far
% left) when h(3) did not have any label).
h(2).YLabel.Position=h(1).YLabel.Position;

% Add context info strings (CIS): Spacecraft position, Earth longitude as text.
[soloStr, earthStr] = solo.qli.context_info_strings(data.solopos, data.earthpos, Tint_24h);
hCisText1 = text(h(10), -0.11, -0.575, soloStr,  'units', 'normalized', 'fontsize', 18);
hCisText2 = text(h(10), -0.11, -0.925, earthStr, 'units', 'normalized', 'fontsize', 18);



xtickangle(h(10),0)
% Add plot information and IRF logo
logopos = h(1).Position;
logopos(1)=logopos(1)+logopos(3)+0.06;
logopos(2)=logopos(2)+0.06;
logopos(3)=0.05;
logopos(4)=logopos(3)*1095/800;
ha2=axes('position',logopos);

if ~isempty(logoPath)
    [x, map]=imread(logoPath);
    image(x)
end
% colormap (map)
set(ha2,'handlevisibility','off','visible','off')
str = solo.qli.utils.generate_data_source_info();
text(h(1), 0, 1.2, str, 'Units', 'normalized')

yyaxis(h(2), 'left');
h(2).YScale = 'log';       % NOTE: Later changed to LIN.
h(2).YTick  = [1, 10, 100];

% NOTE: Panel 2 YTick not auto-adjusted partly because
% solo.qli.utils.ensure_axes_data_tick_margins() can not handle both left &
% right yaxis.
yyaxis(h(2), 'right');
h(2).YScale = 'log';       % NOTE: Later changed to LIN.
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
% Set all YLim and most YTick automatically.
set_YLim_YTick(h([1, 3:4, 6:9]), h([2, 5, 10]))



fig=gcf;
fig.PaperPositionMode='auto';

%===========================
% Save figure to file (24h)
%===========================
% PROPOSAL: Use modify_save_subinterval_plot() here?
%   NOTE: Indirectly calls solo.qli.utils.ensure_axes_data_tick_margins() which
%   thus should be disabled above.
%   ==> Changing order of commands.
%   ==> Calls to YScale, YTick above become superseded.
%   ==> Unwanted change of behaviour.
solo.qli.utils.save_figure_to_file(paths.path_24h, Tint_24h)



%=============================================
% Modify panels 2 & 5, AFTER saving 24 h plot
%=============================================
% Change panel 2+5 y scales to "lin" (previously "log").
% h(5): Keep old ylimits and ticks!
yyaxis(h(2),'right');
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
if ~ALL_PLOTS_ENABLED
    % For debugging/testing.
    %I_6H = [0];
    %I_2H = [0];
    I_6H = [];
    I_2H = [5];
end

%===========================
% Iterate over 6h intervals
%===========================
tBeginSec = solo.qli.utils.log_time('Begin iterating over 6 h intervals', tBeginSec);
for i6h = I_6H
    Tint_6h = Tint_24h(1) + 6*60*60*(i6h+[0, 1]);
    modify_save_subinterval_plot(h, hCisText1, hCisText2, data, Tint_6h, paths.path_6h)
end

%===========================
% Iterate over 2h intervals
%===========================
tBeginSec = solo.qli.utils.log_time('Begin iterating over 2 h intervals', tBeginSec);
for i2h = I_2H
    Tint_2h = Tint_24h(1) + 2*60*60*(i2h+[0, 1]);
    modify_save_subinterval_plot(h, hCisText1, hCisText2, data, Tint_2h, paths.path_2h)
end



close(fig);

[~] = solo.qli.utils.log_time('End of quicklooks_24_6_2_h.m', tBeginSec);

end



% Function to remove duplicated code.
%
% Presumes pre-existing figure with specific axes. Uses customized code to zoom
% in on the sub-time interval and adjusts the y limits for that interval.
function modify_save_subinterval_plot(hAxesArray, hCisText1, hCisText2, data, Tint, parentDirPath)
    assert(isa(hCisText1,  'matlab.graphics.primitive.Text'))
    assert(isa(hCisText2,  'matlab.graphics.primitive.Text'))
    assert(isstruct(data))
    assert(isa(Tint,       'EpochTT') && (length(Tint) == 2))

    irf_zoom(hAxesArray, 'x', Tint);

    %irf_zoom(hAxesArray(1), 'y');
    %adjust_panel_ylimits_N_B(  hAxesArray(2), data,      Tint)
    %adjust_panel_ylimits_Ti(   hAxesArray(5), data.Tpas, Tint)
    %adjust_panel_ylimits_VT_VN(hAxesArray(6), data.Vpas, Tint)
    %adjust_panel_ylimits_ESRF( hAxesArray(8), data.E,    Tint)
    %irf_zoom(hAxesArray(7), 'y');

    % NOTE: Different from for 24h plots.
    set_YLim_YTick(hAxesArray(1:9), hAxesArray(10))

    % Update text
    [hCisText1.String, hCisText2.String] = solo.qli.context_info_strings(data.solopos, data.earthpos, Tint);

    solo.qli.utils.save_figure_to_file(parentDirPath, Tint)
end



% Set y limits and y tick positions.
% Ensure that ticks are not at the min/max (YLim) to avoid overlapping labels.
%
% NOTE: Function can not simultaneously handle both yyaxis left & right.
% NOTE: MATLAB's automatic setting of y ticks for log scale (and which is used)
%       can be bad.
%
% ARGUMENTS
% =========
% hAxesAutoArray
%   Axes for which to set YLim and YTick automatically.
% hAxesManualArray
%   Axes for which to set YLim, but not YTick, automatically.
%
function set_YLim_YTick(hAxesAutoArray, hAxesManualArray)
    % PROPOSAL: Automatically (not MATLAB) set YTick for logarithmic axis to
    %           ensure one tick per power of ten, 10^n.

    assert(isa(hAxesAutoArray,   'matlab.graphics.axis.Axes'))
    assert(isa(hAxesManualArray, 'matlab.graphics.axis.Axes'))
    assert(isempty(intersect(hAxesAutoArray, hAxesManualArray)))

    %=======================================================================
    % Automatically set preliminary YLim (y limits) and final YTick (y tick
    % positions) for selected axes.
    %=======================================================================
    % Set axes y range (YLim) to only cover the data (plus rounding outwards
    % to ticks).
    set(hAxesAutoArray, 'YLimMode', 'auto')
    % Auto-generate ticks (YTick; y values at which there should be ticks).
    set(hAxesAutoArray, 'YTickMode', 'auto')
    %---------------------------------------------------------------------------
    % IMPORTANT: Read YLim without using the return result ("do nothing")
    % --------------------------------------------------------------------
    % IMPLEMENTATION NOTE: THIS COMMAND SHOULD THEORETICALLY NOT BE NEEDED,
    % BUT IS NEEDED FOR THE YLim VALUES TO BE SET PROPERLY. MATLAB BUG?!
    % This behaviour has been observed on Erik P G Johansson's laptop
    % "irony" (MATLAB R2019b, Ubuntu Linux) as of 2023-05-25.
    % Ex: (Re-)scaling of panel 5, 2022-02-23T10-12 (2h plot).
    get(hAxesAutoArray, 'YLim');
    %---------------------------------------------------------------------------
    % Prevent the setting of YLim (next command) from generating new ticks.
    set(hAxesAutoArray,   'YTickMode', 'manual')

    %=========================================================================
    % Automatically set YLim (y limits) but keep old YTick (y tick positions)
    % for selected axes.
    %=========================================================================
    set(hAxesManualArray, 'YTickMode', 'manual')
    %get(hAxesManualArray, 'YLim');   % READ ONLY. UNNECESSARY?
    set(hAxesManualArray, 'YLimMode',  'auto')
    %get(hAxesManualArray, 'YLim');   % READ ONLY. UNNECESSARY?

    %===========================================================================
    % If needed, adjust YLim (but not YTick) to ensure there are margins between
    % YLim and YTick = No ticks on the panel edges/corners.
    %===========================================================================
    hAxesArray = union(hAxesAutoArray, hAxesManualArray);
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