%
% Functions (static methods) associated with bicas.proc.L1L2.cal_RCT using
% different types of RCTs (except for generic reading of RCTs), so that
% bicas.proc.L1L2.cal_RCT does not need to deal with them as special cases
% (almost).
%
%
% DESIGN INTENT
% =============
% This class handles the in-memory modification of calibration data read from
% RCTs (e.g. inverting FTFs) that bicas.proc.L1L2.cal.Cal needs so that generic
% RCT-reading code (bicas.RCT) does not need to.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-19, by moving out functions from
% bicas.proc.L1L2.cal_RCT.
%
classdef cal_RCT_types    
    % PROPOSAL: Class for init_RCT_TYPES_MAP:"Entry" structs.
    % PROPOSAL: Rename.
    %   CON: Does not conform to naming conventions.
    %   PRO: All code is about modifying loaded RCT data.
    % PROPOSAL: Merge with bicas.RCT.



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
        RCT_TYPES_MAP = bicas.proc.L1L2.cal_RCT_types.init_RCT_TYPES_MAP();

    end



    %#####################
    %#####################
    % PRIVATE PROPERTIES
    %#####################
    %#####################
    properties(Access=private, Constant)
        
        % LL = Log Level
        RCT_DATA_LL         = 'debug';
        
    end



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Code to initialize hard-coded static constant.
        %
        function RctTypesMap = init_RCT_TYPES_MAP()
            RctTypesMap = containers.Map();
            
            RctTypesMap('BIAS') = entry(...
                @bicas.RCT.read_BIAS_RCT, ...
                @bicas.proc.L1L2.cal_RCT_types.modify_BIAS_RCT_data, ...
                @bicas.proc.L1L2.cal_RCT_types.log_BIAS_RCT, ...
                'PROCESSING.RCT_REGEXP.BIAS');
            
            RctTypesMap('LFR') = entry(...
                @bicas.RCT.read_LFR_RCT, ...
                @bicas.proc.L1L2.cal_RCT_types.modify_LFR_RCT_data, ...
                @bicas.proc.L1L2.cal_RCT_types.log_LFR_RCTs, ...
                'PROCESSING.RCT_REGEXP.LFR');
            
            RctTypesMap('TDS-CWF') = entry(...
                @bicas.RCT.read_TDS_CWF_RCT, ...
                @bicas.proc.L1L2.cal_RCT_types.modify_TDS_CWF_RCT_data, ...
                @bicas.proc.L1L2.cal_RCT_types.log_TDS_CWF_RCTs, ...
                'PROCESSING.RCT_REGEXP.TDS-LFM-CWF');
            
            RctTypesMap('TDS-RSWF') = entry(...
                @bicas.RCT.read_TDS_RSWF_RCT, ...
                @bicas.proc.L1L2.cal_RCT_types.modify_TDS_RSWF_RCT_data, ...
                @bicas.proc.L1L2.cal_RCT_types.log_TDS_RSWF_RCTs, ...
                'PROCESSING.RCT_REGEXP.TDS-LFM-RSWF');
            
            %###################################################################
            function Entry = entry(...
                    readRctFunc, modifyRctFunc, logRctFunc, ...
                    filenameRegexpSettingKey)
                
                Entry = struct(...
                    ... % Pointer to function that reads one RCT.
                    'readRctFunc',              readRctFunc, ...      
                    ... % Pointer to function that modifies data for one RCT.
                    'modifyRctFunc',            modifyRctFunc, ...    
                    ... % Pointer to function that logs data for one RCT.
                    'logRctFunc',               logRctFunc, ... 
                    ... % Setting key to setting containing regex. for filename.
                    'filenameRegexpSettingKey', filenameRegexpSettingKey);
            end
            %###################################################################
        end
        
        
        
        function RctData = modify_BIAS_RCT_data(RctData)
            
            FtfRctSet = RctData.FtfSet;
            
            % Change name of field (sic!).
            % (There are many fields which are just kept untouched by this
            % function.)
            RctData = rmfield(RctData, 'FtfSet');
            RctData.FtfRctSet = FtfRctSet;
            
            % ASSERTIONS
            nTime = irf.assert.sizes(...
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
                    ItfModifIvpt = bicas.proc.L1L2.cal.utils.extrapolate_tabulated_TF_to_zero_Hz(ItfIvpt);
                    
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
                    itfModifIvpt = @(omegaRps) (bicas.proc.L1L2.cal.utils.eval_tabulated_TF(...
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
        
        
        
        % NOTE: Contains no TFs. Data is therefore trivial to use as it is in
        % the RCT.
        function RctData2 = modify_TDS_CWF_RCT_data(RctData1)
            RctData2 = RctData1;
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
                ItfModifIvpt = bicas.proc.L1L2.cal.utils.extrapolate_tabulated_TF_to_zero_Hz(ItfRctIvpt);
                
                % MODIFY tabulated ITF --> function ITF
                %
                % NOTE: Use zero outside of tabulated frequencies (beyond
                % already made extrapolation). TDS-RSWF data requires
                % extrapolation.
                VALUE_OUTSIDE_TABLE = 0;
                itfModifIvpt = @(omegaRps) (bicas.proc.L1L2.cal.utils.eval_tabulated_TF(...
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
            LL               = bicas.proc.L1L2.cal_RCT_types.RCT_DATA_LL;
            
            %=====================
            % Iterate over EpochL
            %=====================
            for iEpochL = 1:numel(RctData.epochL)
                
                L.logf(LL, 'Below values are used for data beginning %s:', ...
                    irf.cdf.TT2000_to_UTC_str(RctData.epochL(iEpochL)))
                
                % Log bias current calibration
                L.logf(LL, '    BIAS current offsets: %s [aampere]',         ...
                    bicas.proc.L1L2.cal.utils.vector_string(...
                        '% 10e', RctData.Current.offsetsAAmpere(iEpochL, :)))
                L.logf(LL, '    BIAS current gain   : %s [aampere/TM unit]', ...
                    bicas.proc.L1L2.cal.utils.vector_string(...
                        '% 10e', RctData.Current.gainsAapt(     iEpochL, :)))
                
                % Log transfer functions (frequency domain) at selected
                % frequencies.
                L.logf(LL, ...
                    ['    Note: Not logging the exact RCT BIAS TFs', ...
                    ' (FTFs; RctData.FtfRctSet) since the inversion is trivial.'])
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
            irf.assert.sizes(dcDiffOffsetsAVolt, [NaN, 3]);            
            for iEpochH = 1:numel(RctData.epochH)
                
                L.logf(LL, 'Below values are used for data beginning %s:', ...
                    irf.cdf.TT2000_to_UTC_str(RctData.epochH(iEpochH)))
                
                L.logf(LL, ...
                    '    BIAS DC single voltage offsets ( V1, V2, V3): %s [avolt]', ...
                    bicas.proc.L1L2.cal.utils.vector_string('%g', ...
                    RctData.dcSingleOffsetsAVolt(iEpochH, :)))
                L.logf(LL, ...
                    '    BIAS DC diff   voltage offsets (E12,E13,E23): %s [avolt]', ...
                    bicas.proc.L1L2.cal.utils.vector_string('%g', ...
                    dcDiffOffsetsAVolt(iEpochH)))
            end
                
            %###################################################################
            % Nested utility function.
            % NOTE: Implicitly function of iEpochL, L, LL.
            function log_TF(name, freqArray, ItfList)
                bicas.proc.L1L2.cal.utils.log_TF_function_handle(...
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
                    bicas.proc.L1L2.cal.utils.log_TF_tabulated(...
                        bicas.proc.L1L2.cal_RCT_types.RCT_DATA_LL, ...
                        itfName, ...
                        RctData.FtfRctTpivCaCa{iLsf}{iBlts}, ...
                        L);
                    
                    itfIvpt          = RctData.ItfModifIvptCaCa{iLsf}{iBlts};
                    itfName = sprintf('%s ITF (modif., interp.)', itfNamePrefix);
                    
                    bicas.proc.L1L2.cal.utils.log_TF_function_handle(...
                        bicas.proc.L1L2.cal_RCT_types.RCT_DATA_LL, ...
                        itfName, ...
                        'ivolt/TM unit', FREQ_HZ, itfIvpt, L)
                    
                end
            end    % for
            
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_CWF_RCTs(RctData, L)
            
            L.logf(bicas.proc.L1L2.cal_RCT_types.RCT_DATA_LL, ...
                'TDS CWF calibration factors: %s [ivolt/TM]', ...
                bicas.proc.L1L2.cal.utils.vector_string('%g', RctData.factorsIvpt));
        end
        
        
        
        % Analogous to log_BIAS_RCT.
        function log_TDS_RSWF_RCTs(RctData, L)
            % TODO: Log tabulated TFs.
            
            FREQ_HZ = 0;
            
            for iBlts = 1:3
                itfNamePrefix = sprintf('TDS RSWF, BLTS/BIAS_%i, ITF', iBlts);
                
                bicas.proc.L1L2.cal.utils.log_TF_tabulated(...
                    bicas.proc.L1L2.cal_RCT_types.RCT_DATA_LL, ...
                    sprintf('%s (as in RCT)', itfNamePrefix), ...
                    RctData.ItfRctIvptCa{iBlts}, ...
                    L);
                
                bicas.proc.L1L2.cal.utils.log_TF_function_handle(...
                    bicas.proc.L1L2.cal_RCT_types.RCT_DATA_LL, ...
                    sprintf('%s (modif., interp.)', itfNamePrefix), ...
                    'ivolt/TM unit', FREQ_HZ, RctData.itfModifIvptCa{iBlts}, L)
            end
        end

        
        
    end    % methods(Static, Access=private)



end
