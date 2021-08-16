%
% Functions (static methods) associated with bicas.proc.L1L2.cal using RCTs.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-08-16, by moving out functions from bicas.proc.L1L2.cal.
%
classdef cal_RCT   % < handle
    % PROPOSAL: Automatic test code.
    
    
    
    %###################
    %###################
    % PUBLIC PROPERTIES
    %###################
    %###################
    properties(GetAccess=public, Constant)
        
        % Constants relating to each different type of RCT
        % ------------------------------------------------
        % containers.Map: RCT Type ID --> Info about RCT type.
        % Its keys defines the set of RCT Type ID strings.
        RCT_TYPES_MAP = bicas.proc.L1L2.cal_RCT.init_RCT_TYPES_MAP();
        
    end



    %#####################
    %#####################
    % PRIVATE PROPERTIES
    %#####################
    %#####################
    properties(Access=private, Constant)
        
        % LL = Log Level
        READING_RCT_PATH_LL = 'info';
        RCT_DATA_LL         = 'debug';
        
    end



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Code to initialize hard-coded static constant.
        %
        function RctTypesMap = init_RCT_TYPES_MAP()
            RctTypesMap = containers.Map();
            
            % NOTE: Not TF.
            MODIFY_TDS_CWF_DATA_FUNC = @(S) (S);
            
            RctTypesMap('BIAS') = entry(...
                @bicas.RCT.read_BIAS_RCT, ...
                @bicas.proc.L1L2.cal_RCT.modify_BIAS_RCT_data, ...
                @bicas.proc.L1L2.cal_RCT.log_BIAS_RCT, ...
                'PROCESSING.RCT_REGEXP.BIAS');
            
            RctTypesMap('LFR') = entry(...
                @bicas.RCT.read_LFR_RCT, ...
                @bicas.proc.L1L2.cal_RCT.modify_LFR_RCT_data, ...
                @bicas.proc.L1L2.cal_RCT.log_LFR_RCTs, ...
                'PROCESSING.RCT_REGEXP.LFR');
            
            RctTypesMap('TDS-CWF') = entry(...
                @bicas.RCT.read_TDS_CWF_RCT, ...
                MODIFY_TDS_CWF_DATA_FUNC, ...
                @bicas.proc.L1L2.cal_RCT.log_TDS_CWF_RCTs, ...
                'PROCESSING.RCT_REGEXP.TDS-LFM-CWF');
            
            RctTypesMap('TDS-RSWF') = entry(...
                @bicas.RCT.read_TDS_RSWF_RCT, ...
                @bicas.proc.L1L2.cal_RCT.modify_TDS_RSWF_RCT_data, ...
                @bicas.proc.L1L2.cal_RCT.log_TDS_RSWF_RCTs, ...
                'PROCESSING.RCT_REGEXP.TDS-LFM-RSWF');
            
            %###################################################################
            function Entry = entry(readRctFunc, modifyRctFunc, logRctFunc, filenameRegexpSettingKey)
                Entry = struct(...
                    'readRctFunc',              readRctFunc, ...      % Pointer to function that reads one RCT.
                    'modifyRctFunc',            modifyRctFunc, ...    % Pointer to function that modifies data for one RCT.
                    'logRctFunc',               logRctFunc, ...       % Pointer to function that logs data for one RCT.
                    'filenameRegexpSettingKey', filenameRegexpSettingKey);
            end
            %###################################################################
        end
        
        
        
        % Determine the path to the RCT that should be used according to
        % algorithm specified in the documentation(?). If there are multiple
        % matching candidates, choose the latest one as indicated by the
        % filename.
        %
        %
        % IMPLEMENTATION NOTES
        % ====================
        % Useful to have this as separate functionality so that the chosen RCT
        % to use can be explicitly overridden via e.g. settings.
        %
        function path = find_RCT_regexp(rctDir, filenameRegexp, L)

            %=================================================
            % Find candidate files and select the correct one
            %=================================================
            dirObjectList = dir(rctDir);
            dirObjectList([dirObjectList.isdir]) = [];    % Eliminate directories.
            filenameList = {dirObjectList.name};
            % Eliminate non-matching filenames.
            filenameList(~EJ_library.str.regexpf(filenameList, filenameRegexp)) = [];
            
            % ASSERTION / WARNING
            if numel(filenameList) == 0
                % ERROR
                error('BICAS:CannotFindRegexMatchingRCT', ...
                    ['Can not find any calibration file that matches regular', ...
                    ' expression "%s" in directory "%s".'], ...
                    filenameRegexp, rctDir);
            end
            % CASE: There is at least one candidate file.
            
            filenameList = sort(filenameList);
            filename     = filenameList{end};
            path         = fullfile(rctDir, filename);
            
            if numel(filenameList) > 1
                % WARNING/INFO/NOTICE
                msg = sprintf(...
                    ['Found multiple calibration files matching regular', ...
                    ' expression "%s"\n in directory "%s".\n', ...
                     'Selecting the latest one as indicated by', ...
                     ' the filename: "%s".\n'], ...
                    filenameRegexp, rctDir, filename);
                for i = 1:numel(filenameList)
                    msg = [msg, sprintf('    %s\n', filenameList{i})];
                end
                L.log('debug', msg)
            end
            
            % IMPLEMENTATION NOTE: Not logging which calibration file is
            % selected, since this function is not supposed to actually load the
            % content.
        end



        % Load all non-BIAS RCTs (all types) using assumptions on filenames.
        %
        % NOTES
        % =====
        % NOTE: Can be useful for manual experimentation with calibration.
        % NOTE: Necessary when processing L1-->L2 (inofficially) since L1 does
        %       not have CALIBRATION_TABLE+CALIBRATION_TABLE_INDEX.
        % NOTE: Will only load ONE of each RCT type (no potential RCT time
        %       dependence as per global attribute CALIBRATION_TABLE) and
        %       requires user to not use CALIBRATION_TABLE_INDEX.
        %
        % IMPLEMENTATION NOTE: BICAS only needs one non-BIAS RCT type at a time.
        % However, it is useful to be able to initialize bicas.proc.L1L2.cal so
        % that it can simultanteously calibrate all kinds of data for debugging
        % purposes. Therefore loads ALL non-BIAS RCT types.
        %
        %
        % RETURN VALUE
        % ============
        % RctDataMap : containers.Map.
        %       One key per non-BIAS RCT type ID. Value = 1x1 cell array with
        %       RCT data.
        %       IMPLEMENTATION NOTE: Returns containers.Map to provide the same
        %       interface to bicas.proc.L1L2.cal constructor as
        %       bicas.proc.L1L2.cal_RCT.find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE.
        % 
        function RctDataMap = find_read_non_BIAS_RCTs_by_regexp(rctDir, SETTINGS, L)
            
            RctDataMap = containers.Map();
            
            for rctTypeId = {'LFR', 'TDS-CWF', 'TDS-RSWF'}
                
                settingKey     = bicas.proc.L1L2.cal_RCT.RCT_TYPES_MAP(...
                    rctTypeId{1}).filenameRegexpSettingKey;
                filenameRegexp = SETTINGS.get_fv(settingKey);
                filePath       = bicas.proc.L1L2.cal_RCT.find_RCT_regexp(rctDir, filenameRegexp, L);
                RctDataList    = {bicas.proc.L1L2.cal_RCT.read_RCT_modify_log(...
                    rctTypeId{1}, filePath, L)};
                
                % NOTE: Placing all non-BIAS RCT data inside 1x1 cell arrays so
                % that they are stored analogously with when using ga.
                % CALIBRATION_TABLE.
                RctDataMap(rctTypeId{1}) = RctDataList;
            end
        end



        % Load non-BIAS RCT(s) of ONE type (rctTypeId) using CDF global
        % attribute CALIBRATION_TABLE and zVars CALIBRATION_TABLE_INDEX and BW.
        %
        % IMPLEMENTATION NOTE
        % ===================
        % May load MULTIPLE RCTs (of the same RCT type) but will only load those
        % RCTs which are actually needed, as indicated by
        % zVariables CALIBRATION_TABLE_INDEX and BW. This is necessary since
        % CALIBRATION_TABLE may reference unnecessary RCTs of types not
        % recognized by BICAS (LFR's ROC-SGSE_CAL_RCT-LFR-VHF_V01.cdf
        % /2019-12-16), and which are therefore unreadable by BICAS (BICAS would
        % crash).
        %
        %
        % ARGUMENTS
        % =========
        % ga_CALIBRATION_TABLE
        %       1D cell array of strings. LFR/TDS RCT global attribute
        %       CALIBRATION_TABLE.
        % zv_CALIBRATION_TABLE_INDEX
        %       LFR/TDS BICAS input dataset zVariable CALIBRATION_TABLE_INDEX.
        % zv_BW
        %       Either
        %       (1) [] (as for TDS data), or
        %       (2) LFR input dataset zVariable BW.
        %
        %
        % RETURN VALUE
        % ============
        % RctDataMap
        %       containers.Map with
        %           keys   = non-BIAS RCT type ID.
        %           values = 1D cell.
        %       Non-empty indices {iRct} come from
        %       zv_CALIBRATION_TABLE_INDEX(i,1). Each element is the content of
        %       the corresponding RCT mentioned in ga_CALIBRATION_TABLE.
        %       IMPLEMENTATION NOTE: Returns containers.Map to provide the same
        %       interface to bicas.proc.L1L2.cal constructor as
        %       bicas.proc.L1L2.cal_RCT.find_read_non_BIAS_RCTs_by_regexp().
        %
        function RctDataMap = find_read_non_BIAS_RCTs_by_CALIBRATION_TABLE(...
                rctDir, rctTypeId, ...
                ga_CALIBRATION_TABLE, zv_CALIBRATION_TABLE_INDEX, zv_BW, L)
            
            % ASSERTION
            assert(iscell(ga_CALIBRATION_TABLE))
            
            if isempty(zv_BW)
                % ASSERTION
                nCt = EJ_library.assert.sizes(...
                    ga_CALIBRATION_TABLE,       [-1, 1], ...
                    zv_CALIBRATION_TABLE_INDEX, [-2, 2]);
                
                % CT = CALIBRATION_TABLE
                iCtArray = unique(zv_CALIBRATION_TABLE_INDEX(:, 1));
                
            else
                % ASSERTIONS
                nCt = EJ_library.assert.sizes(...
                    ga_CALIBRATION_TABLE,       [-1, 1], ...
                    zv_CALIBRATION_TABLE_INDEX, [-2, 2], ...
                    zv_BW,                      [-2, 1]);
                assert(all(ismember(zv_BW, [0,1])))
                
                iCtArray = unique(zv_CALIBRATION_TABLE_INDEX(logical(zv_BW), 1));
            end
            
            
            
            % Cell array of paths to RCTs of the same RCT type.
            RctDataList = cell(nCt, 1);
            
            % NOTE: Iterate over those entries in CALIBRATION_TABLE that should
            % be considered, NOT all indices. May therefore legitimately leave
            % some cells in cell array empty.
            for i = 1:numel(iCtArray)
                % NOTE: Cell array index is one greater than the stored value.
                j              = iCtArray(i) + 1;
                filePath       = fullfile(rctDir, ga_CALIBRATION_TABLE{j});
                RctDataList{j} = bicas.proc.L1L2.cal_RCT.read_RCT_modify_log(...
                    rctTypeId, filePath, L);
            end
            
            RctDataMap = containers.Map();
            RctDataMap(rctTypeId) = RctDataList;
        end


        
        function RctData = modify_BIAS_RCT_data(RctData)
            
            FtfRctSet = RctData.FtfSet;
            
            % Change name of field (sic!).
            % (There are many fields which are just kept untouched by this
            % function.)
            RctData = rmfield(RctData, 'FtfSet');
            RctData.FtfRctSet = FtfRctSet;
            
            % ASSERTIONS
            nTime = EJ_library.assert.sizes(...
                FtfRctSet.DcSingleAvpiv,   [-1, 1], ...
                FtfRctSet.DcDiffAvpiv,     [-1, 1], ...
                FtfRctSet.AcLowGainAvpiv,  [-1, 1], ...
                FtfRctSet.AcHighGainAvpiv, [-1, 1]);
            
            % NOTE: Derive ITFs.
            ItfSet = [];
            for iTf = 1:nTime
                % INVERT: FTF --> ITF
                
                % Temporary variables which are stored in the definitions of
                % anonymous functions later.
                % * Might speed up code by eliminating calls to method .inverse()
                % * Reduces size of individual expressions.
                TempItfDcSingleAvpiv   = FtfRctSet.DcSingleAvpiv{  iTf}.inverse();
                TempItfDcDiffAvpiv     = FtfRctSet.DcDiffAvpiv{    iTf}.inverse();
                TempItfAcLowGainAvpiv  = FtfRctSet.AcLowGainAvpiv{ iTf}.inverse();
                TempItfAcHighGainAvpiv = FtfRctSet.AcHighGainAvpiv{iTf}.inverse();
                
                ItfSet.dcSingleAvpiv{  iTf} = @(omegaRps) (TempItfDcSingleAvpiv.eval(omegaRps));
                ItfSet.dcDiffAvpiv{    iTf} = @(omegaRps) (TempItfDcDiffAvpiv.eval(omegaRps));
                ItfSet.acLowGainAvpiv{ iTf} = @(omegaRps) (TempItfAcLowGainAvpiv.eval(omegaRps));
                ItfSet.acHighGainAvpiv{iTf} = @(omegaRps) (TempItfAcHighGainAvpiv.eval(omegaRps));
            end
            
            RctData.ItfSet = ItfSet;
            
        end
        
        
            
        function RctData2 = modify_LFR_RCT_data(RctData1)
            
            FtfRctTpivCaCa = RctData1.FtfTpivTable;
            
            % Read LFR FTFs, derive ITFs and modify them.
            itfModifIvptCaCa = {};
            for iLsf = 1:numel(FtfRctTpivCaCa)
                
                itfModifIvptCaCa{end+1} = {};
                for iBlts = 1:numel(FtfRctTpivCaCa{iLsf})

                    % INVERT: tabulated FTF --> tabulated ITF
                    ItfIvpt = FtfRctTpivCaCa{iLsf}{iBlts}.inverse();
                    
                    % MODIFY tabulated ITF
                    ItfModifIvpt = bicas.proc.L1L2.cal_utils.extrapolate_tabulated_TF_to_zero_Hz(ItfIvpt);
                    
                    % MODIFY tabulated ITF --> Function TF
                    %
                    % NOTE: Can not blindly forbid extrapolation (beyond the
                    % extrapolation to 0 Hz already done above) by setting value
                    % outside table=NaN (which deliberately triggers error
                    % elsewhere). LFR's tabulated TFs do in theory cover
                    % frequencies up to the Nyquist frequency, but in practice,
                    % the actual sampling frequency varies slightly. This means
                    % that when the Nyquist frequency also varies slightly and
                    % sometimes it exceeds the tabulated frequencies.
                    % Ex: solo_L1R_rpw-lfr-surv-cwf-e-cdag_20201102_V01.cd
                    %
                    % 2020-11-06: LFR tables (RCT):
                    % F0=24576 Hz: f={  12--12288} [Hz]
                    % F1= 4096 Hz: f={0.01-- 2048} [Hz]
                    % F2=  256 Hz: f={0.01--  128} [Hz]
                    % F3=   16 Hz: f={0.01--    8} [Hz]
                    VALUE_OUTSIDE_TABLE = 0;
                    %VALUE_OUTSIDE_TABLE = NaN;   % Does not work. See comments above.
                    itfModifIvpt = @(omegaRps) (bicas.proc.L1L2.cal_utils.eval_tabulated_TF(...
                        ItfModifIvpt, omegaRps, VALUE_OUTSIDE_TABLE));
                    clear ItfModifIvpt
                    
                    itfModifIvptCaCa{iLsf}{iBlts} = itfModifIvpt;
                end
            end
            
            RctData2 = [];
            % NOTE: RctData.FtfRctTpivCaCa is still kept (for debugging).
            RctData2.FtfRctTpivCaCa   = FtfRctTpivCaCa;    % Just copied.
            RctData2.ItfModifIvptCaCa = itfModifIvptCaCa;
        end
        
        
        
        function RctData2 = modify_TDS_RSWF_RCT_data(RctData1)
            RctData2 = [];
            
            % Modify tabulated TDS-RSWF TFs.
            for iBlts = 1:numel(RctData1.ItfIvptList)
                % NOTE: Overwriting.
                
                ItfRctIvpt = RctData1.ItfIvptList{iBlts};
                
                % Store tabulated ITF EXACTLY AS THEY ARE in the RCT (before
                % modification).
                % NOTE: Struct field does not need to be
                % pre-initialized/pre-allocated.
                RctData2.ItfRctIvptCa{iBlts} = ItfRctIvpt;
                
                % MODIFY __tabulated__ ITF
                % (Does NOT wrap function handle in function handle.)
                ItfModifIvpt = bicas.proc.L1L2.cal_utils.extrapolate_tabulated_TF_to_zero_Hz(ItfRctIvpt);
                
                % MODIFY tabulated ITF --> function ITF
                %
                % NOTE: Use zero outside of tabulated frequencies (beyond
                % already made extrapolation). TDS-RSWF data requires
                % extrapolation.
                VALUE_OUTSIDE_TABLE = 0;
                itfModifIvpt = @(omegaRps) (bicas.proc.L1L2.cal_utils.eval_tabulated_TF(...
                    ItfModifIvpt, omegaRps, VALUE_OUTSIDE_TABLE));
                
                
                RctData2.itfModifIvptCa{iBlts} = itfModifIvpt;
                    
            end
            
        end
        


        % Log some indicative value(s) for a BIAS RCT.
        %
        % NOTE: Does not log file path. Caller is assumed to do that.
        function log_BIAS_RCT(RctData, L)
            
            % Logging parameters
            DC_FREQ_HZ       = [0];   % Single & diffs.
            AC_DIFF_FREQS_HZ = [0, 1000];
            LL               = bicas.proc.L1L2.cal_RCT.RCT_DATA_LL;
            
            %=====================
            % Iterate over EpochL
            %=====================
            for iEpochL = 1:numel(RctData.epochL)
                
                L.logf(LL, 'Below values are used for data beginning %s:', ...
                    EJ_library.cdf.TT2000_to_UTC_str(RctData.epochL(iEpochL)))
                
                % Log bias current calibration
                L.logf(LL, '    BIAS current offsets: %s [aampere]',         bicas.proc.L1L2.cal_utils.vector_string('% 10e', RctData.Current.offsetsAAmpere(iEpochL, :)))
                L.logf(LL, '    BIAS current gain   : %s [aampere/TM unit]', bicas.proc.L1L2.cal_utils.vector_string('% 10e', RctData.Current.gainsAapt(     iEpochL, :)))
                
                % Log transfer functions (frequency domain), selected frequencies.
                L.logf(LL, ...
                    '    Note: Not logging the exact RCT BIAS TFs (FTFs; RctData.FtfRctSet) since the inversion is trivial.')
                log_TF('    BIAS ITF DC single',          DC_FREQ_HZ,       RctData.ItfSet.dcSingleAvpiv)
                log_TF('    BIAS ITF DC diff',            DC_FREQ_HZ,       RctData.ItfSet.dcDiffAvpiv)
                log_TF('    BIAS ITF AC diff, low  gain', AC_DIFF_FREQS_HZ, RctData.ItfSet.acLowGainAvpiv)
                log_TF('    BIAS ITF AC diff, high gain', AC_DIFF_FREQS_HZ, RctData.ItfSet.acHighGainAvpiv)
            end
            
            %=====================
            % Iterate over EpochH
            %=====================
            % NOTE: Must work for multiple CDF records.
            dcDiffOffsetsAVolt = [...
                RctData.DcDiffOffsets.E12AVolt, ...
                RctData.DcDiffOffsets.E13AVolt, ...
                RctData.DcDiffOffsets.E23AVolt];
            EJ_library.assert.sizes(dcDiffOffsetsAVolt, [NaN, 3]);            
            for iEpochH = 1:numel(RctData.epochH)
                L.logf(LL, 'Below values are used for data beginning %s:', ...
                    EJ_library.cdf.TT2000_to_UTC_str(RctData.epochH(iEpochH)))
                
                L.logf(LL, '    BIAS DC single voltage offsets ( V1, V2, V3): %s [avolt]', ...
                    bicas.proc.L1L2.cal_utils.vector_string('%g', ...
                    RctData.dcSingleOffsetsAVolt(iEpochH, :)))
                L.logf(LL, '    BIAS DC diff   voltage offsets (E12,E13,E23): %s [avolt]', ...
                    bicas.proc.L1L2.cal_utils.vector_string('%g', ...
                    dcDiffOffsetsAVolt(iEpochH)))

            end
                
            %###################################################################
            % Nested utility function.
            % NOTE: Implicitly function of iEpochL, L, LL.
            function log_TF(name, freqArray, ItfList)
%                 bicas.proc.L1L2.cal_utils.log_TF_function_handle(...
%                     LL, name, 'avolt/ivolt', freqArray, ...
%                     @(omegaRps) (ItfList{iEpochL}.eval(omegaRps)), L);
                bicas.proc.L1L2.cal_utils.log_TF_function_handle(...
                    LL, name, 'avolt/ivolt', freqArray, ...
                    ItfList{iEpochL}, L);
            end
            %###################################################################
        end


        
        % Analogous to log_BIAS_RCT.
        function log_LFR_RCTs(RctData, L)
            % NOTE: Frequencies may go outside of tabulated data.
            %FREQ_HZ = [0, 1, 5];
            FREQ_HZ = [0, 1, 100];
            
            % CASE: This index corresponds to an actually loaded RCT (some are
            % intentionally empty).
            for iLsf = 1:4
                if iLsf ~= 4   nBltsMax = 5;
                else           nBltsMax = 3;
                end
                
                for iBlts = 1:nBltsMax
                    
                    itfNamePrefix = sprintf('LFR, F%i, BLTS/BIAS_%i', iLsf-1, iBlts);
                    
                    itfName = sprintf('%s FTF (as in RCT)', itfNamePrefix);
                    bicas.proc.L1L2.cal_utils.log_TF_tabulated(...
                        bicas.proc.L1L2.cal_RCT.RCT_DATA_LL, ...
                        itfName, ...
                        RctData.FtfRctTpivCaCa{iLsf}{iBlts}, ...
                        L);
                    
                    itfIvpt          = RctData.ItfModifIvptCaCa{iLsf}{iBlts};
                    itfName = sprintf('%s ITF (modif., interp.)', itfNamePrefix);
                    
                    bicas.proc.L1L2.cal_utils.log_TF_function_handle(...
                        bicas.proc.L1L2.cal_RCT.RCT_DATA_LL, ...
                        itfName, ...
                        'ivolt/TM unit', FREQ_HZ, itfIvpt, L)
                    
                end
            end    % for
            
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_CWF_RCTs(RctData, L)
            
            L.logf(bicas.proc.L1L2.cal_RCT.RCT_DATA_LL, ...
                'TDS CWF calibration factors: %s [ivolt/TM]', ...
                bicas.proc.L1L2.cal_utils.vector_string('%g', RctData.factorsIvpt));
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_RSWF_RCTs(RctData, L)
            % TODO: Log tabulated TFs.
            
            FREQ_HZ = 0;
            
            for iBlts = 1:3
                itfNamePrefix = sprintf('TDS RSWF, BLTS/BIAS_%i, ITF', iBlts);
                
                bicas.proc.L1L2.cal_utils.log_TF_tabulated(...
                    bicas.proc.L1L2.cal_RCT.RCT_DATA_LL, ...
                    sprintf('%s (as in RCT)', itfNamePrefix), ...
                    RctData.ItfRctIvptCa{iBlts}, ...
                    L);
                
%                 bicas.proc.L1L2.cal_utils.log_TF_tabulated(...
%                     bicas.proc.L1L2.cal_RCT.RCT_DATA_LL, ...
%                     sprintf('%s (modified)', itfNamePrefix), ...
%                     RctData.ItfModifIvptList{iBlts}, ...
%                     L);
                
                bicas.proc.L1L2.cal_utils.log_TF_function_handle(...
                    bicas.proc.L1L2.cal_RCT.RCT_DATA_LL, ...
                    sprintf('%s (modif., interp.)', itfNamePrefix), ...
                    'ivolt/TM unit', FREQ_HZ, RctData.itfModifIvptCa{iBlts}, L)
            end
        end
        
        
        
        % For a given RCT file, do the following operations, customized for the
        % type of RCT:
        % (1) read RCT file,
        % (2) modify the content as required (in practice extrapolate TFs), and
        % (3) log it.
        % Effectively wraps the different RCT-reading functions.
        % 
        %
        % IMPLEMENTATION NOTES
        % ====================
        % This method exists to
        % (1) run shared code that should be run when reading any RCT (logging,
        %     modifying data),
        % (2) separate the logging from the RCT-reading code, so that external
        %     code can read RCTs without BICAS.
        %
        %
        % ARGUMENTS
        % =========
        % rctTypeId : String constants representing pipeline and RCT to be read.
        %
        function RctData = read_RCT_modify_log(rctTypeId, filePath, L)
            
            L.logf(bicas.proc.L1L2.cal_RCT.READING_RCT_PATH_LL, ...
                'Reading RCT (rctTypeId=%s): "%s"', rctTypeId, filePath)
            
            readRctFunc   = bicas.proc.L1L2.cal_RCT.RCT_TYPES_MAP(rctTypeId).readRctFunc;
            modifyRctFunc = bicas.proc.L1L2.cal_RCT.RCT_TYPES_MAP(rctTypeId).modifyRctFunc;
            logRctFunc    = bicas.proc.L1L2.cal_RCT.RCT_TYPES_MAP(rctTypeId).logRctFunc;
            
            RctData = readRctFunc(filePath);
            RctData = modifyRctFunc(RctData);
            logRctFunc(RctData, L);
        end
        
        

    end    % methods(Static)

    
    
end
