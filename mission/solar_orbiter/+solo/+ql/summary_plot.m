%
% Class with "standard routines" for producing a summary plot.
% An instance of this class models a figure (with multiple panels).
%
% NOTE: Somewhat experimental at this point, but a similar architecture has been
% used successfully for plots elsewhere.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-10-13
%
classdef summary_plot < handle
    % PROPOSAL: Automatically set panel tags (encapsulate; no arguments).
    % PROPOSAL: Assert unique panel tags.
    %
    % PROPOSAL: Consistent naming of ~Epoch/zvTt2000.
    %
    % PROPOSAL: Merge add_panel_spectrogram_SWF_LSF &
    %                 add_panel_spectrogram_CWF.
    %   CON: Too many differences in functionality.
    %       Ex: Enlarge separate snapshots.
    %
    % PROPOSAL: Merge add_panel_spectrogram_SWF_LSF &
    %                 panel_spectrogram_snapshots.
    %   CON: panel_spectrogram_snapshots designed to also be used for TDS
    %        snapshots (in the future).
    %
    % PROPOSAL: Build add_panel_* functions more using interpret_settings_args() for last arguments).
    %   PRO: Can rely on default values more.
    %
    % PROBLEM: Slight problem that one can not truly add functionality when
    % wrapping an add_panel_* method as one would like to. ~Ideally, one would
    % like to wrap the function handle, not the add_panel_* method.
    %   Ex: Can not fade line color in wrapper.
    %
    % TODO-DECISION: Content of figure title
    %   PROPOSAL: Time range
    %   PROPOSAL: DOY.
    %       ~PROBLEM/NOTE: Could span multiple days, or part of day.
    %
    % PROPOSAL: Eliminate use of dataobj as argument (HK). Let caller create
    %           wrapper functions for doing that instead.
    %
    % PROPOSAL: Let plot_* scripts do more of the customizing via e.g. wrappers.
    %
    % PROPOSAL: Rename. Not (MATLAB) plot, but (MATLAB) figure.



    properties(Constant)

        % NOTE: More samples per spectrum is faster (sic!).
        %N_SAMPLES_PER_SPECTRUM = 512;    % Speed test
        N_SAMPLES_PER_SPECTRUM = 128;    % YK request 2020-02-26.
        
        % Fraction of the (minimum) time distance between snapshots (centers)
        % that will be used for displaying the spectra. E.g. 1.0 means that
        % spectras are adjacent between snapshots (for minimum snapshot
        % distance). Less than 1.0 means that there will be some empty space
        % between snapshot spectras.
        SNAPSHOT_WIDTH_FRACTION  = 0.90;
        
        % Overlap in successive time intervals for separate FFTs.
        % Percent, not fraction.
        %SPECTRUM_OVERLAP_PERCENT = 0;
        SPECTRUM_OVERLAP_PERCENT = 50;
            
        % Colormap used for spectras.
        COLORMAP = load('cmap').cmap;
        
        % 0.0 = Fade to white. 1.0 = No fading.
        C_FADE = 0.7;
        
        LEGEND_TOP_LEFT_POSITION  = [0.02 0.98];
        LEGEND_TOP_RIGHT_POSITION = [0.98 0.98];
        
    end
    
    
    
    properties(SetAccess = private)
        
        figureComplete = false;
        
        % PCFC = Panel Creation Function Call
        pcfcCa = {};
        
        minTt2000 = int64( Inf);
        maxTt2000 = int64(-Inf);
        
        % Values used for auto-sizing the height of each panel.
        panelHeightFixedSizeArray = zeros(0,1);
        panelHeightWeightArray    = zeros(0,1);
    end
    
    
    
    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access = public)
        
        
        
        % Constructor
        function obj = summary_plot()
        end



        % Add panel for one bit-valued array. Intended for HK.
        %
        % ARGUMENTS
        % =========
        % channelName :
        %
        function add_panel_plot_bit_series(obj, panelTag, ...
                zvEpoch, zvB, channelName, panelHeight)
            
            assert(~obj.figureComplete)
            
            obj.add_panel_internal_vars(...
                @() (panel_plot_bit_series()), zvEpoch, panelHeight, 0);
            
            
            
            function hAxes = panel_plot_bit_series()
                
                CHANNEL_NAME_POS = solo.ql.summary_plot.LEGEND_TOP_RIGHT_POSITION;
                
                % ASSERTIONS
                assert(all(ismember(zvB, [0,1])))
                
                if isempty(zvB)
                    zvEpoch(:) = [];
                    channelName = ['(', channelName, ': empty zVar)'];
                end
                
                
                
                % Divide b into b0 and b1, which when plotted as lines, form one
                % line identical to b. This way one can plot separate parts with
                % separate colors. Use NaN for parts that are not plotted by
                % either b0 or b1.
                b = double(zvB);   % Must be double to be able to assign NaN. Renaming.
                clear zvB
                % Identify indices just after b has flipped.
                iJumpUp = 1 + find((b(1:end-1)==0) & (b(2:end)==1));
                iJumpDn = 1 + find((b(1:end-1)==1) & (b(2:end)==0));
                % Create b0, b1.
                b0 = b;
                b1 = b;
                b0(b0 == 1) = NaN;
                b1(b1 == 0) = NaN;
                b0(iJumpUp) = 1;
                b1(iJumpDn) = 0;
                
                % Create time series object
                % -------------------------
                Ts = irf.ts_scalar(zvEpoch, [b0, b1]);
                
                % Plot
                hAxes = irf_panel(panelTag);
                hLines = irf_plot(hAxes, Ts, 'LineWidth', 3);
                assert(numel(hLines) == 2)
                hLines(1).Color = [0,0,1];
                hLines(2).Color = [1,0,0];
                
                solo.ql.summary_plot.fade_color(hLines)
                
                
                
                % NOTE: Command puts the text relative to the specified
                % coordinates in different ways depending on coordinates.
                irf_legend(hAxes, EJ_library.graph.escape_str(channelName), CHANNEL_NAME_POS);
                
                set(hAxes, 'Clipping', 'on')   % Should be default, so not necessary.
                set(hAxes, 'ClippingStyle', 'rectangle')
                set(hAxes, 'yTickLabel', {})
                set(hAxes, 'YGrid', 'off')          % Remove horizontal grid lines.
                %set(hAxes, 'TickLength', [0, 0])    % Remove ticks by setting their length.
                
                set(hAxes, 'YLim', [-0.2, 1.2])
                set(hAxes, 'YTick', [])
                %set(hAxes, 'yTickLabel', {'0', '1'})
            end
        end
        
        
        
        % Add panel for one or multiple scalar time series. This method is meant
        % to be very general and to be called from wrapper methods. It can hence
        % be allowed to have many settings & arguments.
        %
        % ARGUMENTS
        % =========
        % panelTag      :
        % zvData        : Array. (iRecord). ~CWF data.
        % yLabelNonUnit : y label with unit.
        % --
        % tlLegend : Top-left  (TL) legend. Empty if not used.
        % trLegend : Top-right (TR) legend. Empty if not used.
        %            Cell array of strings, one per
        %            scalar time series.
        %
        function add_panel_time_series_CWF_general(obj, panelTag, ...
            zvEpoch, zvData, linesPropCa, axesPropCa, yLabel, tlLegend, trLegend)
            
            assert(~obj.figureComplete)
            assert(nargin == 1+8)
            
            % NOTE: Implicitly an assertion on argument sizes.
            Ts = irf.ts_scalar(zvEpoch, zvData);
            
            obj.add_panel_internal_vars(...
                @() (panel_time_series()), zvEpoch, 0, 1);
            
            
            
            function hAxes = panel_time_series()
                
                hAxes = irf_panel(panelTag);
                hLines = irf_plot(hAxes, Ts);
                
                % TEMPORARY HACK for argument "fade"?
                if ~isempty(linesPropCa) && strcmp(linesPropCa{1}, 'fade')
                    solo.ql.summary_plot.fade_color(hLines)
                    linesPropCa = linesPropCa(2:end);
                end
                
                if ~isempty(linesPropCa)   set(hLines, linesPropCa{:});   end
                if ~isempty(axesPropCa)    set(hAxes,  axesPropCa {:});   end
                
                ylabel(hAxes, yLabel)
                if ~isempty(tlLegend)
                    irf_legend(hAxes, tlLegend, ...
                        solo.ql.summary_plot.LEGEND_TOP_LEFT_POSITION, 'color', 'k')
                end
                if ~isempty(trLegend)
                    irf_legend(hAxes, trLegend, ...
                        solo.ql.summary_plot.LEGEND_TOP_RIGHT_POSITION)
                end
            end
        end
        
        
        
        function add_panel_time_series_CWF(obj, panelTag, ...
            zvEpoch, zvData, yLabel, removeMean)
        
            assert(isscalar(removeMean))
            EJ_library.assert.sizes(...
                zvEpoch, [-1], ...
                zvData,  [-1])
        
            if removeMean
                zvData = zvData - mean(zvData, 'omitnan');
            end
        
            obj.add_panel_time_series_CWF_general(panelTag, ...
                zvEpoch, zvData, {}, {}, yLabel, {}, {})
        end
        
        

        % Add one panel for 3 numeric time series (one for each antenna).
        %
        % ARGUMENTS
        % =========
        % D           : dataobj
        %
        function add_panel_time_series1_HK(obj, D, zvName, linesPropCa, axesPropCa)
            % NOTE: Automatically derives panel tag.
            panelTag = zvName;
            
            zvEpoch = D.data.Epoch.data;
            zvData  = D.data.(zvName).data;            
            
            % CALL INSTANCE METHOD            
            obj.add_panel_time_series_CWF_general(...
                panelTag, zvEpoch, zvData, [{'fade'}, linesPropCa], axesPropCa, ...
                '', {}, {EJ_library.graph.escape_str(zvName)})
        end
        
        
        
        % Add one panel for 3 numeric time series (one for each antenna).
        %
        % ARGUMENTS
        % =========
        % D           : dataobj
        % zVarNamesCa : Length 3 cell array of zVar names.
        %
        function add_panel_time_series3_HK(obj, D, zvNamesCa)
            assert(numel(zvNamesCa) == 3)
            
            % NOTE: Automatically derives panel tag.
            panelTag = zvNamesCa{1};
            
            zvEpoch = D.data.Epoch.data;
            zvData  = [...
                D.data.(zvNamesCa{1}).data, ...
                D.data.(zvNamesCa{2}).data, ...
                D.data.(zvNamesCa{3}).data];
            
            
            
            % CALL INSTANCE METHOD
            obj.add_panel_time_series_CWF_general(panelTag, zvEpoch, zvData, ...
                {'fade'}, {}, '', {}, {EJ_library.graph.escape_str(zvNamesCa{1}), '2', '3'});
        end
        


        % Add panel for time series for SWF data (one snapshot per row) for one
        % specified LFR sampling frequency (LSF). Uses only the samples of the
        % specified sampling frequency.
        %
        % NOTE: Removes mean from each snapshot separately (both DC & AC)
        % NOTE: Can also handle TDS snapshots some day?
        %
        function add_panel_time_series_SWF_LSF(obj, panelTagSignalsStr, ...
                zvEpoch, zvDataCa, zvSamplFreqHz, iLsf, trLegend, removeMean)
            
            nChannels = numel(zvDataCa);
            assert(numel(removeMean) == nChannels)
            
            bRecords = (zvSamplFreqHz == EJ_library.so.constants.LSF_HZ(iLsf));
            samplFreqHz = EJ_library.so.constants.LSF_HZ(iLsf);
            lsfName     = EJ_library.so.constants.LSF_NAME_ARRAY{iLsf};
            
            
            
            zvEpoch     = zvEpoch(bRecords);
            
            % IMPLEMENTATION NOTE: Can not obviously use
            % EJ_library.assert.sizes() to derive nRecords and nSps since both
            % zvEpoch and zvDataCa are in the process of being transformed (LSF
            % subset, SWF-->CWF), and the values need to be derived and used in
            % the middle of that transformation.
            nRecords    = size(zvEpoch, 1);        % NOTE: After selecting records.
            nSps        = size(zvDataCa{1}, 2);    % SPS = Samples Per Snapshot
            assert(nSps >= 2)
            
            zvEpoch  = EJ_library.so.convert_N_to_1_SPR_Epoch(zvEpoch, nSps, ones(nRecords, 1)*samplFreqHz);
            for i = 1:nChannels
                zvData = zvDataCa{i}(bRecords, :);
                
                if removeMean(i)
                    % Remove mean of each snapshot separately.
                    % /YK 2020-10-13
                    %
                    % IMPLEMENTATION NOTE: Ignore NaN so that can handle
                    % varying-length TDS snapshots which pad the end with NaN/fill
                    % value.
                    zvData = zvData - repmat(mean(zvData, 2, 'omitnan'), [1, nSps]);
                end
                
                zvDataCa{i} = EJ_library.so.convert_N_to_1_SPR_redistribute(zvData);
            end
            
            
            
            % CALL INSTANCE METHOD
            obj.add_panel_time_series_CWF_general(...
                sprintf('%s %s time series', panelTagSignalsStr, lsfName), ...
                zvEpoch, [zvDataCa{:}], {}, {}, '[V]', lsfName, trLegend)
        end
        
        

        % Add one panel for spectrogram for ~CWF data.
        %
        % ARGUMENTS
        % =========
        % panelTag      :
        % yLabelNonUnit : y label without unit (unit is at the color bar;
        %                 Assumes "Ts" uses volt).
        %
        function add_panel_spectrogram_CWF(obj, panelTag, zvEpoch, zvData, zvSamplingFreqHz, yLabelNonUnit)
            assert(~obj.figureComplete)
            assert(nargin == 1+5)
            
            obj.add_panel_internal_vars(...
                @() (panel_spectrogram()), zvEpoch, 0, 1);
            
            Ts  = irf.ts_scalar(zvEpoch, zvData);
            
            
            
            function hAxes = panel_spectrogram()
                
                hAxes = irf_panel(panelTag);
                
                % SS = SubSequence
                [iSs1Array, iSs2Array, nSs] = EJ_library.utils.split_by_change(zvSamplingFreqHz);
                SpecrecCa = cell(nSs, 1);
                parfor jSs = 1:nSs    % PARFOR
                    
                    iSsArray = iSs1Array(jSs) : iSs2Array(jSs);
                    
                    SpecrecCa{jSs} = irf_powerfft(Ts(iSsArray), ...
                        solo.ql.summary_plot.N_SAMPLES_PER_SPECTRUM, ...
                        zvSamplingFreqHz(iSs1Array(jSs)));
                end
                
                Specrec = EJ_library.utils.merge_Specrec(SpecrecCa);
                Specrec.p_label = {'log_{10} [V^2/Hz]'};     % Replaces colorbarlabel
                irf_spectrogram(hAxes, Specrec);    % Replaces irf_plot
                
                set(hAxes, 'yscale','log')
                %caxis(hAxes, [-13, -4])
                % NOTE: Adding frequency unit on separate row.
                ylabel(hAxes, {yLabelNonUnit; 'f [Hz]'})
                
                colormap(solo.ql.summary_plot.COLORMAP)
                
                set(hAxes, 'YTick', [0.1, 1, 10, 100, 1e3, 1e4, 1e5])
            end
        end
        
        
        
        % Convenient wrapper around spectrum_panel.
        % Converts from zVar-like variables to what is actually used for plotting.
        %
        function add_panel_spectrogram_SWF_LSF(obj, ...
            panelTagSignalsStr, zvEpoch, zvData, ...
            zvSamplFreqHz, iLsf, trLegend, colLimits)
        
            assert(~obj.figureComplete)
            
            samplFreqHz = EJ_library.so.constants.LSF_HZ(iLsf);
            lsfName     = EJ_library.so.constants.LSF_NAME_ARRAY{iLsf};
            
            bRecords = (zvSamplFreqHz == samplFreqHz);
            
            zvEpoch = zvEpoch(bRecords, :);
            zvData  = zvData( bRecords, :);

            pcfc = @() (solo.ql.summary_plot.panel_spectrogram_snapshots(...
                sprintf('%s %s spectrogram', panelTagSignalsStr, lsfName), ...
                zvEpoch, ...
                zvData, ...
                samplFreqHz, ...
                lsfName, ...
                trLegend, ...
                colLimits));

            obj.add_panel_internal_vars(pcfc, zvEpoch, 0, 1);
        end
        


        % Method that should be called last, after that all panels have been
        % added.
        %
        %
        % ARGUMENTS
        % =========
        % plotTypeStr : Top title string, e.g. "LFR CWF".
        % filePath    : Path to file. Will only use the filename but permits
        %               submitting entire path for convenience.
        %
        function hAxesArray = finalize(obj, plotTypeStr, filePath)
            
            assert(~obj.figureComplete)
            
            nPanels = numel(obj.pcfcCa);
            assert(nPanels > 0, ...
                'Class is configured with zero panels. Can not handle this case.')
            
            irf_plot(nPanels, 'newfigure');
            
            hAxesArray = [];
            for i = 1:nPanels
                funcHandle = obj.pcfcCa{i}();
                hAxesArray(end+1) = funcHandle();
            end
           
            % For aligning MATLAB axes (taking color legends into account).
            irf_plot_axis_align(hAxesArray)
            % For aligning the content of the MATLAB axes.
            irf_zoom(hAxesArray, 'x', irf.tint([obj.minTt2000; obj.maxTt2000]))
            
            % Remove duplicate x labels.
            % Empirically: Must come after irf_zoom.
            % Only needed for HK?!!
            set(hAxesArray(1:end-1), 'XLabel', [])
            
            
            
            %=============================
            % Adjust the height of panels
            %=============================
            % 'Position' : [left bottom width height]. Size and location,
            % excluding a margin for the labels.
            positionCa = get(hAxesArray, 'Position');    % CA = Cell Array
            yPanelArray1      = cellfun(@(x) ([x(2)]), positionCa);
            % Panel height before distributing height segments. Assumes that
            % panels are adjacent to each other.
            heightPanelArray1 = cellfun(@(x) ([x(4)]), positionCa);
            
            heightPanelArray2 = EJ_library.utils.distribute_segments(...
                sum(heightPanelArray1), ...
                obj.panelHeightFixedSizeArray, ...
                obj.panelHeightWeightArray);
            yPanelArray2 = cumsum([heightPanelArray2(2:end); yPanelArray1(end)], 'reverse');
            
            for i = 1:numel(hAxesArray)
                position = positionCa{i};
                position(2) = yPanelArray2(i);
                position(4) = heightPanelArray2(i);
                set(hAxesArray(i), 'InnerPosition', position)
            end
            
            

            solo.ql.summary_plot.set_std_title(plotTypeStr, filePath, hAxesArray(1))
            
            obj.figureComplete = true;
        end



    end    % methods
        
    
    
    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access = private)
        
        
        
        function add_panel_internal_vars(obj, pcfc, zvTt2000, panelHeightFixedSize, panelHeightWeight)
            % PROPOSAL: Use for assert(~obj.figureComplete) (only place).
            %   PROPOSAL: Rename function to make this inclusion more natural: ~add_panel_PCFC
            
            obj.pcfcCa{end+1} = pcfc;
            obj.panelHeightFixedSizeArray(end+1, 1) = panelHeightFixedSize;
            obj.panelHeightWeightArray   (end+1, 1) = panelHeightWeight;
            
            obj.minTt2000 = min([zvTt2000; obj.minTt2000]);
            obj.maxTt2000 = max([zvTt2000; obj.maxTt2000]);
        end
        
        
        
    end    % methods
    
    
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access = private)



        % Set a standardized plot title (calls "title").
        %
        % ARGUMENTS
        % =========
        % plotTypeStr : Top string, e.g. "LFR CWF".
        % filePath    : Path to file. Will only use the filename but permits
        %               submitting entire path for convenience.
        %
        function set_std_title(plotTypeStr, filePath, hTopAxes)
            assert(isscalar(hTopAxes))
            
            [~, basename, suffix] = fileparts(filePath);
            filename = [basename, suffix];
            
            labelTimestamp = datestr(clock, 'yyyy-mm-dd HH:MM:SS');
            title(hTopAxes, {plotTypeStr, EJ_library.graph.escape_str(...
                sprintf('Plot time: %s, %s', labelTimestamp, filename))})
        end
        
        
        
        % Make spectrogram for time series consisting of snapshots (one sampling
        % frequency). Spectrogram for each snapshot is expanded in time for
        % easy-of-use.
        %
        %
        % ARGUMENTS
        % =========
        % tlLegend : Top-left  (TL) legend string.
        % trLegend : Top-right (TR) legend string.
        %
        % SHORTENINGS
        % ===========
        % SS  : SnapShot
        % SSS : SnapShot Spectrogram
        %
        function hAxes = panel_spectrogram_snapshots(...
                panelTag, zvEpoch, zvData, samplingFreqHz, tlLegend, trLegend, colorLimits)
            
            % NOTE: Multiple-row labels causes trouble for the time series'
            % ylabels.
            %
            % IMPLEMENTATION NOTE: Implemented to potentially be modified in the
            % future to handle TDS snapshots that vary in length.
            
            EJ_library.assert.sizes(colorLimits, [1,2])
            
            
            
            TsCa = solo.ql.summary_plot.snapshot_per_record_2_TSeries(...
                zvEpoch, zvData, samplingFreqHz);
            nTs  = numel(TsCa);
            
            hAxes = irf_panel(panelTag);
            
            %====================
            % Calculate spectras
            %====================
            % IMPLEMENTATION NOTE: irf_powerfft is the most time-consuming part
            % of this code.
            %
            % NOTE: Using for-->parfor speeds up plot_LFR_SWF by
            % 29.912231 s-->21.303145 s (irony). /2020-09-04
            %
            SpecrecCa = cell(nTs, 1);
            ssCenterEpochUnixArray = zeros(nTs, 1);
            parfor i = 1:nTs    % PARFOR
                Ts = TsCa{i};
                
                SpecrecCa{i} = irf_powerfft(Ts, ...
                    solo.ql.summary_plot.N_SAMPLES_PER_SPECTRUM, samplingFreqHz, ...
                    solo.ql.summary_plot.SPECTRUM_OVERLAP_PERCENT);
                
                % IMPLEMENTATION NOTE: Later needs the snapshot centers in the
                % same time system as Specrec.t (epoch Unix).
                ssCenterEpochUnixArray(i) = ...
                    (Ts.time.start.epochUnix + Ts.time.stop.epochUnix)/2;
            end
            sssMaxWidthSecArray = ...
                solo.ql.summary_plot.derive_max_spectrum_width(ssCenterEpochUnixArray);
            
            %====================================================================
            % Set the display locations of individual spectras (override
            % defaults). Separately stretch out the collection of spectras that
            % stems from every snapshot.
            % IMPLEMENTATION NOTE: This can not be done in the first loop in
            % order to derive the (minimum) snapshot time distance.
            %====================================================================
            for i = 1:numel(TsCa)
                bKeep(i) = ~isempty(SpecrecCa{i});
                if ~isempty(SpecrecCa{i})
                    %ssLengthSec = TsCa{i}.time.stop.epochUnix - TsCa{i}.time.start.epochUnix;
                    
                    sssWidthSec = sssMaxWidthSecArray(i) * solo.ql.summary_plot.SNAPSHOT_WIDTH_FRACTION;
                    
                    %SpecrecCa{i} = solo.ql.downsample_Specrec(SpecrecCa{i}, 10);    % TEST
                    
                    % Stretch out spectra (for given snapshot) in time to be
                    % ALMOST adjacent between snapshots.
                    % NOTE: Specrec.dt is not set by irf_powerfft so there is no
                    % default value that can be scaled up.
                    % NOTE: Uses original spectrum positions and re-positions
                    % them relative to snapshot center.
                    
                    % Number of timestamps, but also spectras (within snapshot).
                    nTime = numel(SpecrecCa{i}.t);
                    % Distance from SS center to center of first/last FFT.
                    distToSssEdgeT = sssWidthSec/2 - sssWidthSec/(2*nTime);
                    SpecrecCa{i}.t  = ssCenterEpochUnixArray(i) ...
                        + linspace(-distToSssEdgeT, distToSssEdgeT, nTime);
                    SpecrecCa{i}.dt = ones(nTime, 1) * sssWidthSec/(2*nTime);
                end
            end
            
            SpecrecCa(~bKeep) = [];
            Specrec = solo.ql.summary_plot.merge_specrec(SpecrecCa);
            
            Specrec.p_label = {'log_{10} [V^2/Hz]'};    % Replaces colorbarlabel
            irf_spectrogram(hAxes, Specrec);   % Replaces irf_plot
            
            set(hAxes, 'yscale','log')
            
            irf_legend(hAxes, tlLegend, solo.ql.summary_plot.LEGEND_TOP_LEFT_POSITION, 'color', 'k')
            irf_legend(hAxes, trLegend, solo.ql.summary_plot.LEGEND_TOP_RIGHT_POSITION)
            
            colormap(solo.ql.summary_plot.COLORMAP)
            
            caxis(hAxes, colorLimits)
            set(hAxes, 'YTick', [0.1, 1, 10, 100])
        end
        
        
        
        % Convert zVar-like variables for snapshots to cell array of TSeries.
        %
        % ARGUMENTS
        % =========
        % zvEpoch : Nx1 array. 
        % zvData  : NxM array. (iRecord, iSampleWithinSnapshot). 1 record=1 snapshot.
        % TsCa    : {iSnapshot} 1D cell array of TSeries.
        %           IMPLEMENTATION NOTE: Can not(?) be struct array since MATLAB
        %           confuses indexing a TSeries array (with brackets) with some special
        %           TSeries functionality for calling its code with brackets (calling
        %           TSeries' method "subsref").
        %
        % IMPLEMENTATION NOTE: Function is written to some day be easily extended to be
        % used for use with TDS's length-varying snapshots.
        %
        function TsCa = snapshot_per_record_2_TSeries(zvEpoch, zvData, samplingFreqHz)
            %
            % NOTE: No special treatment of snapshots with only NaN.
            
            assert(isscalar(samplingFreqHz))
            [nRecords, nSps] = EJ_library.assert.sizes(...
                zvEpoch, [-1], ...
                zvData,  [-1, -2]);
            %assert(size(zvEpoch, 1) == size(zvData, 1))   % Same number of records
            bicas.proc_utils.assert_zv_Epoch(zvEpoch)
            
            %nRecords = size(zvData, 1);
            %nSps     = size(zvData, 2);
            assert(nSps >= 2)
            
            % Relative timestamps inside CDF record/snapshot.
            epochRelArray = int64([0:(nSps-1)] * 1/samplingFreqHz * 1e9);
            
            TsCa = {};
            for i = 1:nRecords
                epochRecord = zvEpoch(i) + epochRelArray;
               
                TsCa{i}  = irf.ts_scalar(epochRecord, zvData(i, :));
            end
            
        end
        
        
        
        % For every snapshot, return the available width (in time; centered on
        % snapshot center) for displaying the snapshot spectrogram. Time offset
        % and unit unimportant. Argument and return values have same unit.
        %
        function sssMaxWidthArray = derive_max_spectrum_width(ssCenterArray)
            % Use distance to nearest snapshot for each snapshot separately.
            % NOTE: Should NOT be multiplied by two, since using entire distance.
            sssMaxWidthArray = min([Inf; diff(ssCenterArray(:))], [diff(ssCenterArray(:)); Inf]);
            
            % Use smallest distance between any two consecutive snapshots (one
            % global value for all snapshots). Sometimes yields too narrow
            % spectrograms.
            % Ex: solo_L2_rpw-lfr-surv-swf-e-cdag_20200228_V01.cdf
            %sssMaxWidthArray = min(diff(ssCenterArray)) * ones(size(ssCenterArray));
            
            % NOTE: Can not assume that both input and output have same size, only same
            % length.
            assert(numel(ssCenterArray) == numel(sssMaxWidthArray))
        end
        
        
        
        % Merge multiple instances of "specrec" structs as returned by
        % irf_powerfft, with identical frequencies.
        % NOTE: Optionally added fields must be added after merging.
        % NOTE: Cf EJ_library.utils.merge_Specrec which is more powerful but
        % which is unnecessary here since only merging spectras with the same
        % frequencies (underlying data uses the same sampling frequency).
        %
        % ARGUMENTS
        % =========
        % SpecrecCa : Cell array of "Specrec" structs as returned by irf_powerfft, but
        %             with .dt (column array) added to it.
        %             NOTE: Requires dt (column array of scalars).
        %             NOTE: Assumes that all specrec use the same frequencies.
        %             IMPLEMENTATION NOTE: Uses cell array instead of struct array to be
        %             able to handle (and ignore) the case specrec = [] which can be
        %             returned by irf_powerfft.
        %
        % RETURN VALUE
        % ============
        % Specrec   : Struct array that can be used by irf_spectrogram.
        %
        function Specrec = merge_specrec(SpecrecCa)
            % PROPOSAL: Assertion for frequencies.
            
            Specrec.f  = [];
            Specrec.p  = {[]};   % NOTE: Must 1x1 cell array. The array INSIDE the cell array is added to.
            Specrec.t  = [];
            Specrec.dt = [];
            
            for i = 1:numel(SpecrecCa)
                
                S = SpecrecCa{i};
                if ~isempty(S)
                    EJ_library.assert.struct(S, {'f', 'p', 't', 'dt'}, {});
                    assert(iscolumn(S.dt), 'S.dt is not a column.')
                    assert(numel(S.dt) == numel(S.t), 'Badly formatted SpecrecCa{%i}.', i)
                    
                    Specrec.f    = S.f;                       % NOTE: Not adding to array, but setting it in its entirety.
                    Specrec.p{1} = [Specrec.p{1}; S.p{1}];    % NOTE: Add to array inside cell array.
                    Specrec.t    = [Specrec.t;    S.t(:)];    % NOTE: Has to be column vector.
                    Specrec.dt   = [Specrec.dt;   S.dt(:)];   % NOTE: Has to be column vector.
                end
            end
            
            assert(issorted(Specrec.t))   % Not sure if good assertion.
        end
        
        
        
        function fade_color(hArray)
            
            for i = 1:numel(hArray)
                legendColor = get(hArray(i), 'Color');
                legendColor = 1 - solo.ql.summary_plot.C_FADE*(1-legendColor);
                
                set(hArray(i), 'Color', legendColor);   % Fade color (move toward white).
            end
        end



    end    % methods
    
    
    
end    % classdef
