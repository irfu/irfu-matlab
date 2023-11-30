%
% Functions (static methods) associated with bicas.proc.L1L2.cal.rct.findread using
% different types of RCTs (except for generic reading of RCTs), so that
% bicas.proc.L1L2.cal.rct.findread does not need to deal with them as special cases
% (almost).
%
%
% DESIGN INTENT
% =============
% This class handles the in-memory modification of calibration data read from
% RCTs (e.g. inverting FTFs) that bicas.proc.L1L2.cal.Cal needs so that generic
% RCT-reading code (bicas.proc.L1L2.cal.rct.fs) does not need to.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-19.
%
classdef typeproc
    % PROPOSAL: Class for init_RCT_TYPES_MAP:"Entry" structs.



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
        RCT_TYPES_MAP = bicas.proc.L1L2.cal.rct.typeproc.init_RCT_TYPES_MAP();
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
                @bicas.proc.L1L2.cal.rct.RctTypeBias.read_RCT, ...
                @bicas.proc.L1L2.cal.rct.RctTypeBias.modify_RCT_data, ...
                @bicas.proc.L1L2.cal.rct.RctTypeBias.log_RCT, ...
                bicas.proc.L1L2.cal.rct.RctTypeBias.filenameRegexpSettingKey);
            
            RctTypesMap('LFR') = entry(...
                @bicas.proc.L1L2.cal.rct.RctTypeLfr.read_RCT, ...
                @bicas.proc.L1L2.cal.rct.RctTypeLfr.modify_RCT_data, ...
                @bicas.proc.L1L2.cal.rct.RctTypeLfr.log_RCT, ...
                'PROCESSING.RCT_REGEXP.LFR');
            
            RctTypesMap('TDS-CWF') = entry(...
                @bicas.proc.L1L2.cal.rct.RctTypeTdsCwf.read_RCT, ...
                @bicas.proc.L1L2.cal.rct.RctTypeTdsCwf.modify_RCT_data, ...
                @bicas.proc.L1L2.cal.rct.RctTypeTdsCwf.log_RCT, ...
                'PROCESSING.RCT_REGEXP.TDS-LFM-CWF');
            
            RctTypesMap('TDS-RSWF') = entry(...
                @bicas.proc.L1L2.cal.rct.RctTypeTdsRswf.read_RCT, ...
                @bicas.proc.L1L2.cal.rct.RctTypeTdsRswf.modify_RCT_data, ...
                @bicas.proc.L1L2.cal.rct.RctTypeTdsRswf.log_RCT, ...
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
        
        
        

    end    % methods(Static, Access=private)



end
