%
% Functions (static methods) associated with bicas.proc.L1L2.cal.Cal finding,
% reading, and logging RCTs so that bicas.proc.L1L2.cal.Cal does not need to.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-16, by moving out functions from bicas.proc.L1L2.cal.Cal.
%
classdef cal_RCT    
    % PROPOSAL: Rename.
    %   CON: Does not conform to naming conventions.
    %   NOTE: Mostly about finding paths to RCTs to load. The exception is 
    %         function for actually loading an RCT.
    %
    % PROPOSAL: Normalize L1 & L1R by creating fake ga_CALIBRATION_TABLE,
    %           zv_CALIBRATION_TABLE_INDEX for L1.
    %   PRO: Can eliminate internal special cases in bicas.proc.L1L2.cal.Cal.
    %   NOTE: There is a function
    %         bicas.proc.L1L2.normalize_CALIBRATION_TABLE_INDEX() used by
    %         bicas.proc.L1L2.lfr.process_normalize_CDF() to produce NaN-valued
    %         zVar.
    %
    %   PROPOSAL: Have LFR&TDS:process_normalize_CDF() normalize
    %            CALIBRATION_TABLE_INDEX and CALIBRATION_TABLE by creating fake
    %            values (not NaN) corresponding to the actual values used wrt.
    %            RCTs (:,1) and maybe also index (:,2).
    %            produce_L1R_to_L2_LFR/TDS() then reads the values AFTER those
    %            functions and use them to initialize bicas.proc.L1L2.cal.Cal
    %            object.
    %       CON: Does not enable that much simplification. The only code that
    %            can be simplified is:
    %                (1) bicas.proc.L1L2.cal.Cal.calibrate_voltage_all()
    %                       selects between the CALIBRATION_TABLE_INDEX(:,1) and 0 depending on setting.
    %                (2) bicas.proc.pf.produce_L1R_to_L2_LFR/DS()
    %                       constructs bicas.proc.L1L2.cal.Cal object.
    %           NOTE: There is also some code associated with the settings:
    %               PROCESSING.L1R.LFR.USE_GA_CALIBRATION_TABLE_RCTS
    %               PROCESSING.L1R.LFR.USE_ZV_CALIBRATION_TABLE_INDEX2
    %               PROCESSING.L1R.TDS.CWF.USE_GA_CALIBRATION_TABLE_RCTS
    %               PROCESSING.L1R.TDS.RSWF.USE_GA_CALIBRATION_TABLE_RCTS
    %
    % PROPOSAL: Normalize TDS & LFR by creating a fake zv_BW for TDS.



    %#####################
    %#####################
    % PRIVATE PROPERTIES
    %#####################
    %#####################
    properties(Access=private, Constant)
        
        % LL = Log Level
        READING_RCT_PATH_LL = 'info';
        
    end



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



        % Load one RCT per selected RCT type using assumptions on filenames.
        %
        %
        % NOTES
        % =====
        % NOTE: Can be useful for manual experimentation with calibration of L1R
        %       (and L1) data.
        % NOTE: Necessary when processing L1-->L2 (unofficially) since L1 does
        %       not have CALIBRATION_TABLE+CALIBRATION_TABLE_INDEX.
        % NOTE: Will only load ONE of each RCT type (no potential RCT time
        %       dependence as per global attribute CALIBRATION_TABLE) and
        %       requires user to not use CALIBRATION_TABLE_INDEX.
        %
        % IMPLEMENTATION NOTE: BICAS only needs one non-BIAS RCT type at a time.
        % However, it is useful to be able to initialize bicas.proc.L1L2.cal.Cal so
        % that it can simultanteously calibrate all kinds of data for debugging
        % purposes. Therefore loads ALL non-BIAS RCT types.
        %
        %
        % RETURN VALUE
        % ============
        % RctDataMap
        %       containers.Map. Can be used for bicas.proc.L1L2.cal.Cal
        %       constructor even if there is no zVar CALIBRATION_TABLE.
        %       One key per specified RCT type ID in argument rctTypeIdCa.
        %       Exactly one RCT per RCT type.
        % 
        function RctDataMap = find_read_RCTs_by_regexp(...
                rctTypeIdCa, rctDir, SETTINGS, L)
            
            assert(iscell(rctTypeIdCa))
            
            RctDataMap = containers.Map();
            
            for i = 1:numel(rctTypeIdCa)
                rctTypeId = rctTypeIdCa{i};
                
                settingKey     = bicas.proc.L1L2.cal_RCT_types.RCT_TYPES_MAP(...
                    rctTypeId).filenameRegexpSettingKey;
                filenameRegexp = SETTINGS.get_fv(settingKey);
                filePath       = bicas.proc.L1L2.cal_RCT.find_RCT_regexp(...
                    rctDir, filenameRegexp, L);
                RctDataList    = {bicas.proc.L1L2.cal_RCT.read_RCT_modify_log(...
                    rctTypeId, filePath, L)};
                
                % NOTE: Placing all non-BIAS RCT data inside 1x1 cell arrays so
                % that they are stored analogously with when using ga.
                % CALIBRATION_TABLE.
                RctDataMap(rctTypeId) = RctDataList;
            end
        end



        % (1) Load one BIAS RCT by regular expression.
        % (2) Load one or multiple non-BIAS RCT(s) of the selected type
        % (rctTypeId) using CDF global attribute CALIBRATION_TABLE and ZVs
        % CALIBRATION_TABLE_INDEX and BW.
        %
        %
        % IMPLEMENTATION NOTE
        % ===================
        % May load MULTIPLE RCTs of the same RCT type, but will only load those
        % RCTs which are actually needed, as indicated by
        % zVariables CALIBRATION_TABLE_INDEX and BW. This is necessary since
        % CALIBRATION_TABLE may reference unnecessary RCTs of types not
        % recognized by BICAS (LFR's ROC-SGSE_CAL_RCT-LFR-VHF_V01.cdf
        % /2019-12-16), and which are therefore unreadable by BICAS (BICAS will
        % crash).
        %
        %
        % ARGUMENTS
        % =========
        % nonBiasRctTypeId
        %       RCT type string constant for a non-BIAS RCT.
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
        %       Returns containers.Map that can be used for bicas.proc.L1L2.cal.Cal
        %       constructor.
        %
        function RctDataMap = find_read_RCTs_by_regexp_and_CALIBRATION_TABLE(...
                nonBiasRctTypeId, rctDir, ...
                ga_CALIBRATION_TABLE, ...
                zv_CALIBRATION_TABLE_INDEX, ...
                zv_BW, SETTINGS, L)            
            
            BiasRctDataMap = bicas.proc.L1L2.cal_RCT.find_read_RCTs_by_regexp(...
                {'BIAS'}, rctDir, SETTINGS, L);
            
            RctDataList = bicas.proc.L1L2.cal_RCT.find_read_RCTs_by_CALIBRATION_TABLE(...
                nonBiasRctTypeId, rctDir, ...
                ga_CALIBRATION_TABLE, ...
                zv_CALIBRATION_TABLE_INDEX, ...
                zv_BW, L);
            
            RctDataMap                   = containers.Map();
            RctDataMap('BIAS')           = BiasRctDataMap('BIAS');
            RctDataMap(nonBiasRctTypeId) = RctDataList;            
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
        % --
        % NOTE: Only public due to automatic testing.
        %
        function path = find_RCT_regexp(rctDir, filenameRegexp, L)

            %=================================================
            % Find candidate files and select the correct one
            %=================================================
            dirObjectList = dir(rctDir);
            dirObjectList([dirObjectList.isdir]) = [];    % Eliminate directories.
            filenameList = {dirObjectList.name};
            % Eliminate non-matching filenames.
            filenameList(~irf.str.regexpf(filenameList, filenameRegexp)) = [];
            
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
        

        
    end    % methods(Static)

    
    
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)
        
        
        
        % ARGUMENTS
        % =========
        % zv_BW
        %       LFR L1/L1R zVar BW. If it does not exist (if processing TDS
        %       L1/L1R), then the caller should (!) create a fake one and submit
        %       it (normalize input).
        %
        function RctDataList = find_read_RCTs_by_CALIBRATION_TABLE(...
                nonBiasRctTypeId, rctDir, ...
                ga_CALIBRATION_TABLE, ...
                zv_CALIBRATION_TABLE_INDEX, ...
                zv_BW, L)
            
            % CT = glob.attr. CALIBRATION_TABLE
            
            % ASSERTION            
            assert(iscell(ga_CALIBRATION_TABLE))
            nCt = irf.assert.sizes(...
                ga_CALIBRATION_TABLE,       [-1, 1], ...
                zv_CALIBRATION_TABLE_INDEX, [-2, 2], ...
                zv_BW,                      [-2, 1]);
            assert(all(ismember(zv_BW, [0,1])))

            % Obtain indices into glob.attr. CALIBRATION_TABLE
            % ------------------------------------------------
            % NOTE: May exclude some in zv_CALIBRATION_TABLE_INDEX due to zv_BW.
            iCtArray = unique(zv_CALIBRATION_TABLE_INDEX(logical(zv_BW), 1));

            
            
            % Cell array of paths to RCTs of the same RCT type.
            RctDataList = cell(nCt, 1);
            
            % IMPLEMENTATION NOTE: Iterate over those entries in
            % CALIBRATION_TABLE that should be CONSIDERED, i.e. NOT all indices.
            % May therefore legitimately leave some cells in cell array empty.
            for i = 1:numel(iCtArray)
                % NOTE: Cell array index is one greater than the stored value.
                j              = iCtArray(i) + 1;
                filePath       = fullfile(rctDir, ga_CALIBRATION_TABLE{j});
                RctDataList{j} = bicas.proc.L1L2.cal_RCT.read_RCT_modify_log(...
                    nonBiasRctTypeId, filePath, L);
            end
            
        end
        
        
        
        % For a given RCT file, do the following operations, customized for the
        % type of RCT:
        %   (1) read RCT file,
        %   (2) modify the content as required (in practice extrapolate TFs),
        %       and
        %   (3) log the modified RCT content.
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
        % rctTypeId
        %       String constants representing pipeline and RCT to be read.
        %
        function RctData = read_RCT_modify_log(rctTypeId, filePath, L)
            
            L.logf(bicas.proc.L1L2.cal_RCT.READING_RCT_PATH_LL, ...
                'Reading RCT (rctTypeId=%s): "%s"', rctTypeId, filePath)
            
            RCT_TYPES_MAP = bicas.proc.L1L2.cal_RCT_types.RCT_TYPES_MAP;
            
            readRctFunc   = RCT_TYPES_MAP(rctTypeId).readRctFunc;
            modifyRctFunc = RCT_TYPES_MAP(rctTypeId).modifyRctFunc;
            logRctFunc    = RCT_TYPES_MAP(rctTypeId).logRctFunc;
            
            RctDataTemp = readRctFunc(filePath);
            RctData     = modifyRctFunc(RctDataTemp);
            logRctFunc(RctData, L);
        end
        
        

    end    % methods(Static, Access=private)

    
    
end
