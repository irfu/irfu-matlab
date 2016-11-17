% data_manager - Class that does the actual processing of datasets.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-10
%
%
%
% BASIC CONCEPT
% =============
% The basic concept is inspired by irfu-matlab's MMS data manager. The user gives the class the input data ("elementary
% input process data") it has, and asks for the output data ("elementary output process data") it wants.
%
% The class maintains an internal set of different "process data" variables (PDVs), a set of variables where each one is
% uniquely referenced by a unique string, a "process data ID" (PDID). These process data variables are all related to each
% other in a conceptual "web" (an acyclic directed graph) describing their dependencies on each other. Initially these
% process data variables are all empty.
%
%
% Dependencies
% -------------
% For a PDV "Y", there is a list of lists(!) of PDVs (or PDIDs) X_ij. To derive Y, you need to have available at least
% one X_ij, for every i. 
% 
%
% Execution/process
% -----------------
% One or a few PDVs are simply set at the outset. After that, the desired output PDV is requested. When a PDV
% "Y" is requested (and it has not already been derived and stored), the data manager tries to recursively request
% enough PDVs {X_ij} that could be used to later derive Y. After that those requests have been satisfied, Y is finally
% derived from the available PDVs {X_ij}.
%
%
% Example
% -------
% (Example assumes unique sets of required PDVs for simplicity. Information flows from left to right, i.e. PDVs on the
% left are used to derive the ones to the right. Each string represents a PDV.)
%    input_1 ---------------------------------------------- output_1 --
%                                                                      \
%                                                                       -- output_2
%
%    input_2 ------ intermediate_1 ------------------------ output_3
%                /                 \
%    input_3 ---|                   \
%                \                   \
%    input_4 -------------------------- intermediate_2 ---- output_4
%
%
% Advantages with architecture
% ----------------------------
% - Easy to implement CHANGING S/W modes (as defined by the RCS ICD) using this class, althoguh the class itself is
%   unaware of S/W modes.
%      Ex: Updates to the RPW pipeline, datasets.
%      Ex: Can use a different set of modes when validating using the official RCS validation software at LESIA (not all
%      datasets are available for input there, and one gets error when the S/W descriptor describes modes which require
%      dataset IDs not defined/existing there).
%      Ex: Implement inofficial modes for testing(?).
% - Easier to keep support for legacy datasets (older versions of CDFs).
% - Easier to maintain a varying set of test/debugging modes?
% - Indirect "caching".
% - Can reuse processing code?
% - PDVs can be derived from different sets of PDVs which is useful since many datasets (CDF files) are similar.
%
%
%
% DEFINITIONS OF TERMS, ABBREVIATIONS
% ===================================
% - Process data (PD)
%       In practice, one internal variable representing some form of data in some step of
%       the "data processing process".
% - Elementary input/output (EIn/EOut) process data
%       Process data that represents data that goes in and out of the data manager as a whole. The word "elementary" is
%       added to distinguish input/output from the input/output for a smaller operation, in particular when converting
%       to/from intermediate process data.
%       Elementary input (EIn)  process data :
%           Supplied BY the user. In practice this should correspond to the content of a CDF file.
%       Elementary output(EOut) process data :
%           Returned TO the user. In practice this should correspond to the content of a CDF file.
% - Intermediate process data
%       Process data that is not "elementary" input/output process data.
%       These only exist inside data_manager.
% - Process data ID (PDID)
%       A string that uniquely represents (refers to) a type of process data.
%       By convention, the PDIDs for elementary input/output PDIDs are a shortened version of
%       dataset_ID+skeleton_version_str, e.g. "L2S_LFR-SURV-CWF-E_V01". Hardcoded values of these occur throughout the
%       code as constants.
% - Processing function
%       Function that accepts PDs as arguments and from them derive another PD.
% - Pre-Demuxing-Calibration Data PreDCD)
%       Generic data format that can represent all forms of input datasets before demuxing and calibration.
%       Can use an arbitrary number of samples per record.
%       Consists of struct with fields:
%           .Epoch
%           .ACQUISITION_TIME
%           .DemuxerInput : struct with fields.
%               BIAS_1 to .BIAS_5  : NxM arrays, where M may be 1 (1 sample/record) or >1.
%           .freqHz                : Snapshot frequency in Hz. Unimportant for one sample/record data.
%           .DIFF_GAIN
%           .MUX_SET
%           QUALITY_FLAG
%           QUALITY_BITMASK
%           DELTA_PLUS_MINUS
%           (.SAMP_DTIME  ?)
%       Fields are "CDF-like": rows=records, all have same number of rows.
% - Post-Demuxing-Calibration Data (PostDCD)
%       Like PreDCD but with additional fields. Tries to capture a superset of the information that goes into any
%       dataset produced by BICAS.
%       Has extra fields:
%           .DemuxerOutput   : struct with fields.
%               V1, V2, V3,   V12, V13, V23,   V12_AC, V13_AC, V23_AC.
%           .IBIAS1
%           .IBIAS2
%           .IBIAS3
% 
% NOTE: So far (2016-11-11), every dataset ID is mapped (one-to-one) to a PDID. In principle, in the future, multiple
% elementary PDVs could each be associated with the one and same dataset ID, if the different datasets (CDF files) (with
% the same dataset ID) fulfill different roles. Example: Produce one output CDF covering the time intervals of data in
% input CDFs (all with the same dataset ID).
%
%
%
% CODE CONVENTIONS 
% ================
% - It is implicit that arrays/matrices representing CDF data, or "CDF-like" data, use the first MATLAB array index to
%   represent CDF records.
%
%
%
% IMPLEMENTATION NOTES
% ====================
% 1) As of now (2016-10-11), the "basic concept" is only partly taken advantage of, and is partly unnecessary.
% 2) This class is intended to:
%     - As far as possible, not be aware of S/W modes as defined in the RCS ICD.
%     - Not deal with writing CDFs. It does however read CDF files as a way of setting elementary input process data but
%       only since it is trivial.
%
classdef data_manager < handle     % Explicitly declare it as a handle class to avoid IDE warnings.
%#######################################################################################################################
% PROPOSAL: Use other class name that implices processing, and fits with defined term "process data".
%     "processing_manager"? "proc_manager"?
% --
% PROPOSAL: Derive DIFF_GAIN (from BIAS HK using time interpolation) in one code common to both LFR & TDS.
%   PROPOSAL: Function
%   PROPOSAL: In intermediate PDV?!
%   PRO: Uses flag for selecting interpolation time in one place.
% PROPOSAL: Derive HK_BIA_MODE_MUX_SET (from BIAS SCI or HK using time interpolation for HK) in one code common to both LFR & TDS.
%   PROPOSAL: Function
%   PROPOSAL: In intermediate PDV?!
%   PRO: Uses flag for selecting HK/SCI DIFF_GAIN in one place.
%   PRO: Uses flag for selecting interpolation time in one place.
%
% PROPOSAL: Move out calibration (not demuxing) from data_manager.
%   PROPOSAL: Reading of calibration files.
%   PROPOSAL: Function for calibrating with either constant factors and transfer functions. (Flag for choosing which.)
%       NOTE: Function needs enough information to split up data into sequences on which transfer functions can be applied.
%
% PROPOSAL: Use double for all numeric zVariables in the processing. Do not produce or require proper type, e.g. integers, in any
%           intermediate processing. Only convert to the proper data type/class when writing to CDF.
%   PRO: Variables can keep NaN to represent fill/pad value, also for "integers".
%   PRO: The knowledge of the dataset CDF formats is not spread out over the code.
%       Ex: Setting default values for PreDcd.QUALITY_FLAG, PreDcd.QUALITY_BITMASK, PreDcd.DELTA_PLUS_MINUS.
%       Ex: ACQUISITION_TIME.
%   CON: Less assertions can be made in utility functions.
%       Ex: dm_utils.ACQUISITION_TIME_*, dm_utils.tt2000_* functions.
%   CON: ROUNDING ERRORS. Can not be certain that values which are copied, are actually copied.
%   --
%   NOTE: Functions may in principle require integer math to work correctly.
%--
% PROPOSAL: Function for checking if list of elementary input PDIDs for a list of elementary output PDIDs can be satisfied.
%
% PROPOSAL: Functions for classifying PDIDs.
%   Ex: LFR/TDS, sample/record or snapshot/record, V01/V02 (for LFR variables changing name).
%   Ex: Number of samples per snapshot (useful for "virtual snapshots") in record (~size of record).
%   PRO: Could reduce/collect the usage of hardcoded PDIDs.
%   NOTE: Mostly/only useful for elementary input PDIDs?
%   PROPOSAL: Add field (sub-struct) "format_constants" in list of EIn PDIDs in constants?!!
%       NOTE: Already has version string in list of EIn PDIDs.
%--
% NOTE: Both BIAS HK and LFR SURV CWF contain MUX data (only LFR has one timestamp per snapshot). True also for other input datasets?
%#######################################################################################################################

    properties(Access=private)
        
        % Initialized by constructor.
        ProcessDataVariables = containers.Map('KeyType', 'char', 'ValueType', 'any');

        % Local constants
        INTERMEDIATE_PDIDs = {'PreDCD', 'PostDCD'};
        ALL_PDIDs          = [];  % Initialized by constructor.
        
    end
    
    %###################################################################################################################

    methods(Access=public)
        
        function obj = data_manager()
        % CONSTRUCTOR
        
            global CONSTANTS
        
            obj.ALL_PDIDs = {...
                CONSTANTS.INPUTS_PDIDS_LIST{:}, ...
                CONSTANTS.OUTPUTS_PDIDS_LIST{:}, ...
                obj.INTERMEDIATE_PDIDs{:}};
            
            obj.validate
        end

        
        
        function set_elementary_input_process_data(obj, pdid, processData)
        % Set elementary input process data directly as a struct.
        %
        % NOTE: This method is a useful interface for test code.
        
            global CONSTANTS        
            CONSTANTS.assert_EIn_PDID(pdid)            
            obj.set_process_data_variable(pdid, processData);
        end



        function processData = get_process_data_recursively(obj, pdid)
        % Get process data by deriving it from other process data recursively.
        %
        % Try in the following order
        % 1) If the process data is already available (stored), then return it.
        % 2) If there is a process function, try to derive process data from other process data recursively.
        % 3) Return nothing.
        %
        % RETURN VALUE
        % ============
        % processData : Empty if data was not already available, or could not be derived recursively.        
        %               ==> Empty is returnd for elementary input process data, or possibly for misconfigurations (error).
        

            % NOTE: This log message is particularly useful for following the recursive calls of this function.
            % sw_mode_ID comes first since that tends to make the log message values line up better.
            irf.log('n', sprintf('Begin function (pdid=%s)', pdid))

            processData = obj.get_process_data_variable(pdid);

            %============================================
            % Check if process data is already available
            %============================================
            if ~isempty(processData)
                % CASE: Process data is already available.
                irf.log('n', sprintf('End   function (pdid=%s) - PD already available', pdid))
                return   % NOTE: This provides caching so that the same process data are not derived twice.
            end

            % CASE: Process data is NOT already available.

            %=============================================================
            % Obtain processing function, and the PDs which are necessary
            %=============================================================
            [InputPdidsMap, processingFunc] = obj.get_processing_info(pdid);
            % ASSERTION
            if isempty(processingFunc)
                processData = [];
                irf.log('n', sprintf('End   function (pdid=%s) - PD can not be derived (is EIn)', pdid))
                return
            end

            %==============================================
            % Obtain the input process datas - RECURSIVELY
            %==============================================
            Inputs          = containers.Map();
            inputFieldsList = InputPdidsMap.keys();
            
            for iField = 1:length(inputFieldsList)
                inputField         = inputFieldsList{iField};
                inputFieldPdidList = InputPdidsMap(inputField);
                
                % See if can obtain any one PDV for the PDIDs listed (for this "input field").
                for iPdid = 1:length(inputFieldPdidList)
                
                    % NOTE: Should return error if can not derive data.
                    inputPd = obj.get_process_data_recursively(inputFieldPdidList{iPdid});   % NOTE: RECURSIVE CALL.
                    
                    if ~isempty(inputPd)
                        % CASE: Found an available PDV.
                        Input.pd   = inputPd;
                        Input.pdid = inputFieldPdidList{iPdid};
                        Inputs(inputField) = Input;
                        break;
                    end
                end
                
                if ~Inputs.isKey(inputField)
                    error('BICAS:data_manager:Assertion:SWModeProcessing', ...
                        'Can not derive necessary process data for pdid=%s, input field=%s', pdid, inputField)
                end

            end
            
            %========================================================
            % Derive the actual process data from other process data
            %========================================================
            irf.log('n', sprintf('Begin deriving PDID=%s using %s', pdid, func2str(processingFunc)))
            processData = processingFunc(Inputs);
            
            obj.set_process_data_variable(pdid, processData)
            
            % NOTE: This log message is useful for being able to follow the recursive calls of this function.
            irf.log('n', sprintf('End   function (pdid=%s) - PD was derived', pdid))
            
        end
        
        
    
        function SwModeInfo = get_extended_sw_mode_info(obj, cliParameter)
        % Return ~constants structure for a specific S/W mode as referenced by a CLI parameter.
        % 
        % This structure is automatically put together from other structures (constants) to avoid having to define too
        % many redundant constants.
        %
        % NOTE: This function takes what constants.SW_MODES_INFO_LIST returns and adds info about input and output datasets to it.
        % NOTE: The function takes the CLI_parameter as parameter, not the S/W mode ID!
        %
        % IMPLEMENTATION NOTE: The reason for that this function is located in data_manager instead of constants is the
        % call to "bicas.data_manager.get_elementary_input_PDIDs" and that constants.m should contain no reference to
        % data_manager.m (for initialization reasons at the very least). AMENDMENT 2016-11-16: Function no longer call
        % that function but might call a similar one in the future. Comment therefore still relevant.

            global CONSTANTS

            SwModeInfo = bicas.utils.select_structs(CONSTANTS.SW_MODES_INFO_LIST, 'CLI_parameter', {cliParameter});            
            SwModeInfo = SwModeInfo{1};

            % Collect all associated elementary input PDIDs.
            %input_PDIDs = obj.get_elementary_input_PDIDs(C_sw_mode.output_PDIDs, C_sw_mode.ID);

            try
                SwModeInfo.inputs = bicas.utils.select_structs(CONSTANTS.INPUTS_INFO_LIST,  'PDID', SwModeInfo.input_PDIDs);
            catch exception
                error('BICAS:Assertion:IllegalConfiguration', ...
                    'Can not identify all input PDIDs associated with S/W mode CLI parameter "%s".', cliParameter)
            end
            try
                SwModeInfo.outputs = bicas.utils.select_structs(CONSTANTS.OUTPUTS_INFO_LIST, 'PDID', SwModeInfo.output_PDIDs);
            catch exception
                error('BICAS:Assertion:IllegalConfiguration', ...
                    'Can not identify all output PDIDs associated with S/W mode CLI parameter "%s".', cliParameter)
            end
        end
        
        
        
    end   % methods: Instance, public
    
    %###################################################################################################################

    methods(Access=private)
        
        function validate(obj)
            
            global CONSTANTS
            
            % Initialize instance variable "process_data" with keys (with corresponding empty values).
            for pdid = obj.ALL_PDIDs
                obj.ProcessDataVariables(pdid{1}) = [];
            end
            
            % Validate
            %bicas.utils.assert_strings_unique(obj.ELEMENTARY_INPUT_PDIDs)
            %bicas.utils.assert_strings_unique(obj.ELEMENTARY_OUTPUT_PDIDs)
            bicas.utils.assert_strings_unique(obj.INTERMEDIATE_PDIDs)
            bicas.utils.assert_strings_unique(obj.ALL_PDIDs)
            
            
            
            %========================
            % Iterate over S/W modes
            %========================
            for i = 1:length(CONSTANTS.SW_MODES_INFO_LIST)
                
                % Get info for S/W mode.
                sw_mode_CLI_parameter = CONSTANTS.SW_MODES_INFO_LIST{i}.CLI_parameter;
                SwModeInfo = obj.get_extended_sw_mode_info(sw_mode_CLI_parameter);
                
                % Check unique inputs CLI_parameter.
                inputs_CLI_parameters = cellfun(@(s) ({s.CLI_parameter}), SwModeInfo.inputs);
                bicas.utils.assert_strings_unique(inputs_CLI_parameters)
                
                % Check unique outputs JSON_output_file_identifier.
                jsonOutputFileIdentifiers = cellfun(@(s) ({s.JSON_output_file_identifier}), SwModeInfo.outputs);

                bicas.utils.assert_strings_unique(jsonOutputFileIdentifiers)
            end
            
        end   % validate
        
        
        
        function assert_PDID(obj, pdid)
            if ~ischar(pdid) || ~ismember(pdid, obj.ALL_PDIDs)
                error('BICAS:data_manager:Assertion', 'There is no such PDID="\%s".', pdid)
            end
        end



        %==================================================================================
        % Return process data value.
        %
        % Does NOT try to fill process data variable with (derived) data if empty.
        % Can be used to find out whether process data has been set.
        %
        % ASSERTS: Valid PDID.
        % DOES NOT ASSERT: Process data has already been set.
        %==================================================================================
        function processData = get_process_data_variable(obj, pdid)
            obj.assert_PDID(pdid)
            
            processData = obj.ProcessDataVariables(pdid);
        end



        %==================================================================================
        % Set a process data variable.
        %
        % ASSERTS: Process data has NOT been set yet.
        % ASSERTS: Valid PDID.
        %==================================================================================
        function set_process_data_variable(obj, pdid, processData)
            % ASSERTIONS
            obj.assert_PDID(pdid)
            if ~isempty(obj.ProcessDataVariables(pdid))
                error('BICAS:data_manager:Assertion', 'There is already process data for the specified PDID="%s".', pdid);
            end
            
            obj.ProcessDataVariables(pdid) = processData;
        end

        
        
        % EXPERIMENTAL
        function [InputPdidsMap, processingFunc] = get_processing_info(obj, outputPdid)
        % For every PDID, return meta-information needed to derive the actual process data.
        %
        % HIGH-LEVEL DESCRIPTION
        % ======================
        % This function effectively does two things:
        % 1) (IMPORTANT) It defines/describes the dependencies between different PDIDs (an ~acyclic directed
        %    graph), by returning
        %    a) the "processing" function needed to produce the process data, and
        %    b) the PDIDs needed by the processing function.
        % 2) (LESS IMPORTANT) It implicitly decides how to informally distinguish the input process datas by giving them
        %    suitable struct field names which the processing functions recognize them with.
        %
        % The method itself is NOT recursive but it is designed so that other code can recurse over
        % the implicit "~acyclic graph" defined by it.
        % Other code can use this method to two things recursively:
        % (1) Derive process data.
        % (2) Derive all the necessary elementary input PDIDs needed to produce an elementary output
        %     PDID (without generating any process data). This is useful for automatically generating the required input
        %     dataset IDs for S/W modes (needed for (a) generating the S/W descriptor, and (b) generating requirements
        %     on CLI parameters).
        %
        % ARGUMENT AND RETURN VALUES
        % ==========================
        % pdid        : The PDID (string) for the PD that should be produced.
        % inputPdids  :
        %     Struct where every field is set to a cell array of PDIDs. In every such cell array, only one of its
        %     PDIDs/PDs is necessary for the processing function. The field names of the fields are "human-readable".
        %     NOTE: Elementary input PDIDs yield an empty struct (no field names).
        % processingFunc :
        %     Pointer to a function that can derive process data (for output_PDID) from other process data (for
        %     input_PDIDs).
        %         process_data = processing_func(inputs)
        %         ARGUMENTS:
        %             inputsMap : A containers.Map
        %                <keys>       : The same keys as InputPdidsMap (returned by get_processing_info).
        %                <values>     : A struct with fields .pd (Process data) and .pdid (PDID for .pd).
        %     Empty function is equivalent to that outputPdid is an EIn-PDID.
        
        % PROPOSAL: Change name of return value: InputPdidsLists? Something to imply a list of lists of PDIDs. *Map?
            
            obj.assert_PDID(outputPdid)
            
            % Assign value used for the case of elementary INPUT PDID (no input PDs <==> no fields).
            InputPdidsMap = containers.Map('KeyType', 'char', 'ValueType', 'any');
            
            

            % Assign value used for the case of there being no real processing function
            % -------------------------------------------------------------------------
            % IMPLEMENTATION NOTE: Choice of value for the case of there being no processing function:
            % 1) Using an empty value ==> The caller can immediately check whether it has received a processing function
            % and then immediately throw an assertion error if that was unexpected.
            % 2) Bad alternative: Use a nonsense function (one that always gives an assertion error, albeit a good error message)
            % ==> The (assertion) error will be triggered first when the attempting to call the nonexisting processing
            % function. ==> Error first very late. ==> Bad alternative.
            processingFunc = [];
            

            
            switch(outputPdid)
                %=====================================================
                % Elementary INPUT PDIDs : Return empty default value
                %=====================================================
                % NOTE: It is still useful to include cases which do nothing since it is a check on permitted values
                % (so that switch-otherwise can give error).
                % BIAS
                case 'HK_BIA_V01'
                % LFR
                case 'L2R_LFR-SBM1-CWF_V01'
                case 'L2R_LFR-SBM1-CWF_V02'
                case 'L2R_LFR-SBM2-CWF_V01'
                case 'L2R_LFR-SBM2-CWF_V02'
                case 'L2R_LFR-SURV-CWF_V01'
                case 'L2R_LFR-SURV-CWF_V02'
                case 'L2R_LFR-SURV-SWF_V01'
                case 'L2R_LFR-SURV-SWF_V02'
                % TDS
                case 'L2R_TDS-LFM-CWF_V01'                    
                case 'L2R_TDS-LFM-RSWF_V01'
                case 'L2R_TDS-LFM-RSWF_V02'
                    
                %====================
                % Intermediary PDIDs
                %====================
                
                case 'PreDCD'
                    InputPdidsMap('HK_cdf') = {'HK_BIA_V01'};
                    InputPdidsMap('SCI_cdf') = {...
                    'L2R_LFR-SBM1-CWF_V01', ...    % LFR    
                    'L2R_LFR-SBM1-CWF_V02', ...
                    'L2R_LFR-SBM2-CWF_V01', ...
                    'L2R_LFR-SBM2-CWF_V02', ...
                    'L2R_LFR-SURV-CWF_V01', ...
                    'L2R_LFR-SURV-CWF_V02', ...
                    'L2R_LFR-SURV-SWF_V01', ...
                    'L2R_LFR-SURV-SWF_V02', ...
                    'L2R_TDS-LFM-CWF_V01', ...     % TDS
                    'L2R_TDS-LFM-RSWF_V01', ...
                    'L2R_TDS-LFM-RSWF_V02'};
                    processingFunc = @bicas.data_manager.process_LFR_to_PreDCD;
                    
                case 'PostDCD'
                    InputPdidsMap('PreDCD') = {'PreDCD'};
                    processingFunc = @bicas.data_manager.process_demuxing_calibration;

                %=========================
                % Elementary OUTPUT PDIDs
                %=========================

                %-----
                % LFR
                %-----
                case {'L2S_LFR-SBM1-CWF-E_V02', ...
                      'L2S_LFR-SBM2-CWF-E_V02', ...
                      'L2S_LFR-SURV-CWF-E_V02', ...
                      'L2S_LFR-SURV-SWF-E_V02'}
                    InputPdidsMap('PostDCD') = {'PostDCD'};
                    processingFunc     = @(inputs) (bicas.data_manager.process_PostDCD_to_LFR(inputs, outputPdid));

                %-----
                % TDS
                %-----
                %case 'L2S_TDS-LFM-CWF-E_V02'
                %case 'L2S_TDS-LFM-RSWF-E_V02'                    

                %===========================================
                % OTHERWISE: Has no implementation for PDID
                %===========================================
                otherwise
                    % NOTE: The PDID can be valid without there being an implementation for it, i.e. the error should
                    % NOT be replaced with assert_PDID().
                    error('BICAS:data_manager:Assertion:OperationNotImplemented', ...
                        'This function has no implementation for this PDID, "%s".', outputPdid)
            end   % switch

        end   % get_processing_info
       
        
        
        
    end   % methods: Instance, private

    
    
    %###################################################################################################################

    
    
    methods(Static, Access=public)

        function PreDcd = process_LFR_to_PreDCD(InputsMap)
        % Processing function. Convert LFR CDF data (PDs) to PreDCD.
        %
        % Keeps number of samples/record. Treats 1 samples/record "length-one snapshots".
        
            global CONSTANTS
        
            % PROBLEM: Hardcoded CDF data types (MATLAB classes).

            SciPd   = InputsMap('SCI_cdf').pd;
            HkPd    = InputsMap('HK_cdf').pd;
            sciPdid = InputsMap('SCI_cdf').pdid;

            nRecords = size(SciPd.Epoch, 1);
            
            %=====================================================================
            % Handle differences between skeletons V01 and V02
            % ------------------------------------------------
            % LFR_V, LFR_E: zVars with different names (but identical meaning).
            % L1_REC_NUM  : Not defined in V01, but in V02 dataset skeletons.
            %=====================================================================
            switch(sciPdid)
                case {  'L2R_LFR-SBM1-CWF_V01', ...
                        'L2R_LFR-SBM2-CWF_V01', ...
                        'L2R_LFR-SURV-CWF_V01', ...
                        'L2R_LFR-SURV-SWF_V01'}
                    POTENTIAL  = SciPd.POTENTIAL;
                    ELECTRICAL = SciPd.ELECTRICAL;
                    L1_REC_NUM = NaN * zeros(nRecords, 1);   % Set to fill values.
                case {  'L2R_LFR-SBM1-CWF_V02', ...
                        'L2R_LFR-SBM2-CWF_V02', ...
                        'L2R_LFR-SURV-CWF_V02', ...
                        'L2R_LFR-SURV-SWF_V02'}
                    POTENTIAL  = SciPd.V;
                    ELECTRICAL = SciPd.E;
                    L1_REC_NUM = SciPd.L1_REC_NUM;
                otherwise
                    error('BICAS:data_manager:Assertion:SWModeProcessing:ConfigurationBug', ...
                        'Can not handle PDID="%s"', sciPdid)
            end
            
            %========================================================================================
            % Handle differences between datasets with and without zVAR FREQ:
            % LFR_FREQ: Corresponds to FREQ only defined in some LFR datasets.
            %========================================================================================
            switch(sciPdid)
                case {  'L2R_LFR-SBM1-CWF_V01', ...
                        'L2R_LFR-SBM1-CWF_V02'}
                    FREQ = ones(nRecords, 1) * 1;   % Always value "1".
                case {  'L2R_LFR-SBM2-CWF_V01', ...
                        'L2R_LFR-SBM2-CWF_V02'}
                    FREQ = ones(nRecords, 1) * 2;   % Always value "2".
                case {  'L2R_LFR-SURV-CWF_V01', ...
                        'L2R_LFR-SURV-CWF_V02', ...
                        'L2R_LFR-SURV-SWF_V01', ...
                        'L2R_LFR-SURV-SWF_V02'}
                    FREQ = SciPd.FREQ;
                otherwise
                    error('BICAS:data_manager:Assertion:ConfigurationBug', ...
                        'Can not handle PDID="%s"', sciPdid)
            end
            
            nSamplesPerRecord = size(POTENTIAL, 2);
            freqHz            = bicas.dm_utils.get_LFR_frequency( FREQ );   % NOTE: Needed also for 1 SPR.
            
            % Find the relevant value of zVariables R0, R1, R2, "R3".
            Rx = bicas.dm_utils.get_LFR_Rx( SciPd.R0, SciPd.R1, SciPd.R2, FREQ );   % NOTE: Function also handles the imaginary zVar "R3".
            
            PreDcd = [];
            PreDcd.Epoch = SciPd.Epoch;
            PreDcd.ACQUISITION_TIME = SciPd.ACQUISITION_TIME;
            PreDcd.DELTA_PLUS_MINUS = bicas.dm_utils.derive_DELTA_PLUS_MINUS(freqHz, nSamplesPerRecord);            
            PreDcd.freqHz           = freqHz;
            PreDcd.SAMP_DTIME       = bicas.dm_utils.derive_SAMP_DTIME(freqHz, nSamplesPerRecord);
            
            PreDcd.L1_REC_NUM = L1_REC_NUM;
            % Replace illegally empty data with fill values/NaN
            % -------------------------------------------------
            % IMPLEMENTATION NOTE: QUALITY_FLAG, QUALITY_BITMASK have been found empty in test data, but should have
            % attribute DEPEND_0 = "Epoch" ==> Should have same number of records as Epoch.
            % Can not save CDF with zVar with zero records (crashes when reading CDF). ==> Better create empty records.
            % Test data: MYSTERIOUS_SIGNAL_1_2016-04-15_Run2__7729147__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-SWF_7729147_CNE_V01.cdf
            PreDcd.QUALITY_FLAG    = SciPd.QUALITY_FLAG;
            PreDcd.QUALITY_BITMASK = SciPd.QUALITY_BITMASK;
            if isempty(PreDcd.QUALITY_FLAG)
                irf.log('w', 'QUALITY_FLAG from the SCI source dataset is empty. Filling with empty values.')
                PreDcd.QUALITY_FLAG    = NaN * zeros([nRecords, 1]);
            end
            if isempty(PreDcd.QUALITY_BITMASK)
                irf.log('w', 'QUALITY_BITMASK from the SCI source dataset is empty. Filling with empty values.')
                PreDcd.QUALITY_BITMASK = NaN * zeros([nRecords, 1]);
            end
            
            PreDcd.DemuxerInput        = [];
            PreDcd.DemuxerInput.BIAS_1 = POTENTIAL;
            PreDcd.DemuxerInput.BIAS_2 = bicas.dm_utils.filter_rows( ELECTRICAL(:,:,1), Rx==1 );
            PreDcd.DemuxerInput.BIAS_3 = bicas.dm_utils.filter_rows( ELECTRICAL(:,:,2), Rx==1 );
            PreDcd.DemuxerInput.BIAS_4 = bicas.dm_utils.filter_rows( ELECTRICAL(:,:,1), Rx==0 );
            PreDcd.DemuxerInput.BIAS_5 = bicas.dm_utils.filter_rows( ELECTRICAL(:,:,2), Rx==0 );
            
            

            % Define local convenience variables. AT = ACQUISITION_TIME
            hkAtTt2000  = bicas.dm_utils.ACQUISITION_TIME_to_tt2000(  HkPd.ACQUISITION_TIME );
            sciAtTt2000 = bicas.dm_utils.ACQUISITION_TIME_to_tt2000( SciPd.ACQUISITION_TIME );
            hkEpoch     = HkPd.Epoch;
            sciEpoch    = SciPd.Epoch;
            
            
            
            %=========================================================================================================
            % 1) Convert time to something linear in time that can be used for processing (not storing time to file).
            % 2) Effectively also chooses which time to use for the purpose of processing:
            %    ACQUISITION_TIME or Epoch.
            %=========================================================================================================
            if CONSTANTS.C.PROCESSING.USE_AQUISITION_TIME_FOR_HK_TIME_INTERPOLATION 
                hkInterpolationTimeTt2000  = hkAtTt2000;
                sciInterpolationTimeTt2000 = sciAtTt2000;
            else
                hkInterpolationTimeTt2000  = hkEpoch;
                sciInterpolationTimeTt2000 = sciEpoch;
            end



            %=========================================================================================================
            % Choose where to get MUX_SET from: LFR-SCI, or BIAS-HK
            % -----------------------------------------------------
            % NOTE: Only obtains one MUX_SET per record ==> Can not change MUX_SET in the middle of a record.
            %=========================================================================================================            
            PreDcd.MUX_SET = bicas.dm_utils.nearest_interpolate_float_records(...
                double(HkPd.HK_BIA_MODE_MUX_SET), hkInterpolationTimeTt2000, sciInterpolationTimeTt2000);   % Use BIAS HK.
            %PreDcd.MUX_SET = LFR_cdf.BIAS_MODE_MUX_SET;    % Use LFR SCI.



            %=========================================================================================================
            % Derive approximate DIFF_GAIN values for from BIAS HK
            %
            % NOTE: Not perfect handling of time when 1 snapshot/record, since one should ideally use time stamps
            % for every LFR _sample_.
            %=========================================================================================================
            PreDcd.DIFF_GAIN = bicas.dm_utils.nearest_interpolate_float_records(...
                double(HkPd.HK_BIA_DIFF_GAIN), hkInterpolationTimeTt2000, sciInterpolationTimeTt2000);


            
            % ASSERTIONS
            bicas.dm_utils.assert_unvaried_N_rows(PreDcd);
            bicas.dm_utils.assert_unvaried_N_rows(PreDcd.DemuxerInput);

            
            
            %==================================================================
            % Log time intervals to enable comparing available SCI and HK data
            %==================================================================
            bicas.dm_utils.log_tt2000_interval('HK  ACQUISITION_TIME', hkAtTt2000)
            bicas.dm_utils.log_tt2000_interval('SCI ACQUISITION_TIME', sciAtTt2000)
            bicas.dm_utils.log_tt2000_interval('HK  Epoch           ', hkEpoch)
            bicas.dm_utils.log_tt2000_interval('SCI Epoch           ', sciEpoch)
        end

        

        function assert_PreDCD(PreDcd)
            FIELDS = {'Epoch', 'ACQUISITION_TIME', 'DemuxerInput', 'freqHz', 'DIFF_GAIN', 'MUX_SET', 'QUALITY_FLAG', ...
                'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'L1_REC_NUM', 'SAMP_DTIME'};
            
            if ~isstruct(PreDcd) || ~isempty(setdiff(fieldnames(PreDcd), FIELDS))
                error('BICAS:data_manager:Assertion:SWModeProcessing', 'PDV structure is not on "PreDCD format".')
            end
            bicas.dm_utils.assert_unvaried_N_rows(PreDcd);
        end
        
        
        
        function assert_PostDCD(PostDcd)
            FIELDS = {'Epoch', 'ACQUISITION_TIME', 'DemuxerInput', 'freqHz', 'DIFF_GAIN', 'MUX_SET', 'QUALITY_FLAG', ...
                'QUALITY_BITMASK', 'DELTA_PLUS_MINUS', 'DemuxerOutput', 'IBIAS1', 'IBIAS2', 'IBIAS3', 'L1_REC_NUM', 'SAMP_DTIME'};
            
            if ~isstruct(PostDcd) || ~isempty(setdiff(fieldnames(PostDcd), FIELDS))
                error('BICAS:data_manager:Assertion:SWModeProcessing', 'PDV structure is not on "PostDCD format".')
            end
            bicas.dm_utils.assert_unvaried_N_rows(PostDcd);
        end
        
        

        function PostDcd = process_demuxing_calibration(InputsMap)
        % Processing function. Converts PreDCD to PostDCD, i.e. demux and calibrate data.
        
            PreDcd = InputsMap('PreDCD').pd;
            bicas.data_manager.assert_PreDCD(PreDcd);
                    
            % Log messages
            for f = fieldnames(PreDcd.DemuxerInput)'
                bicas.dm_utils.log_values_summary(f{1}, PreDcd.DemuxerInput.(f{1}));
            end
            
            PostDcd = PreDcd;
            
            % DEMUX
            PostDcd.DemuxerOutput = bicas.data_manager.simple_demultiplex(...
                PreDcd.DemuxerInput, PreDcd.MUX_SET, PreDcd.DIFF_GAIN);
            
            % Log messages
            for f = fieldnames(PostDcd.DemuxerOutput)'
                bicas.dm_utils.log_values_summary(f{1}, PostDcd.DemuxerOutput.(f{1}));
            end
            
            % BUG / TEMP: Set default values since the real values are not available.
            PostDcd.IBIAS1 = NaN * zeros(size(PostDcd.DemuxerOutput.V1));
            PostDcd.IBIAS2 = NaN * zeros(size(PostDcd.DemuxerOutput.V2));
            PostDcd.IBIAS3 = NaN * zeros(size(PostDcd.DemuxerOutput.V3));
            
            bicas.data_manager.assert_PostDCD(PostDcd)
        end
        

        
        function EOutPD = process_PostDCD_to_LFR(InputsMap, eoutPDID)
        % Processing function. Convert PostDCD to any one of several similar dataset PDs.
        
            PostDcd = InputsMap('PostDCD').pd;
            EOutPD = [];
            
            nSamplesPerRecord = size(PostDcd.DemuxerOutput.V1, 2);   % Samples per record.
            
            switch(eoutPDID)
                case  {'L2S_LFR-SBM1-CWF-E_V02', ...
                       'L2S_LFR-SBM2-CWF-E_V02', ...
                       'L2S_LFR-SURV-CWF-E_V02'}
                    
                    %=====================================================================
                    % Convert 1 snapshot/record --> 1 sample/record (if not already done)
                    %=====================================================================
                    EOutPD.Epoch = bicas.dm_utils.convert_N_to_1_SPR_Epoch( ...
                        PostDcd.Epoch, ...
                        nSamplesPerRecord, ...
                        PostDcd.freqHz  );
                    EOutPD.ACQUISITION_TIME = bicas.dm_utils.convert_N_to_1_SPR_ACQUISITION_TIME(...
                        PostDcd.ACQUISITION_TIME, ...
                        nSamplesPerRecord, ...
                        PostDcd.freqHz  );
                    
                    EOutPD.DELTA_PLUS_MINUS = bicas.dm_utils.convert_N_to_1_SPR_redistribute( PostDcd.DELTA_PLUS_MINUS );
                    EOutPD.L1_REC_NUM       = bicas.dm_utils.convert_N_to_1_SPR_repeat(       PostDcd.L1_REC_NUM,      nSamplesPerRecord);
                    EOutPD.QUALITY_FLAG     = bicas.dm_utils.convert_N_to_1_SPR_repeat(       PostDcd.QUALITY_FLAG,    nSamplesPerRecord);
                    EOutPD.QUALITY_BITMASK  = bicas.dm_utils.convert_N_to_1_SPR_repeat(       PostDcd.QUALITY_BITMASK, nSamplesPerRecord);
                    % F_SAMPLE, SAMP_DTIME: Omitting. Are not supposed to be present in BIAS CWF datasets.
                    
                    % Convert PostDcd.DemuxerOutput
                    for fn = fieldnames(PostDcd.DemuxerOutput)'
                        PostDcd.DemuxerOutput.(fn{1}) = bicas.dm_utils.convert_N_to_1_SPR_redistribute( ...
                            PostDcd.DemuxerOutput.(fn{1}) );
                    end
                    EOutPD.IBIAS1           = bicas.dm_utils.convert_N_to_1_SPR_redistribute( PostDcd.IBIAS1 );
                    EOutPD.IBIAS2           = bicas.dm_utils.convert_N_to_1_SPR_redistribute( PostDcd.IBIAS2 );
                    EOutPD.IBIAS3           = bicas.dm_utils.convert_N_to_1_SPR_redistribute( PostDcd.IBIAS3 );
                    EOutPD.V(:,:,1)         = PostDcd.DemuxerOutput.V1;
                    EOutPD.V(:,:,2)         = PostDcd.DemuxerOutput.V2;
                    EOutPD.V(:,:,3)         = PostDcd.DemuxerOutput.V3;
                    EOutPD.E(:,:,1)         = PostDcd.DemuxerOutput.V12;
                    EOutPD.E(:,:,2)         = PostDcd.DemuxerOutput.V13;
                    EOutPD.E(:,:,3)         = PostDcd.DemuxerOutput.V23;
                    EOutPD.EAC(:,:,1)       = PostDcd.DemuxerOutput.V12_AC;
                    EOutPD.EAC(:,:,2)       = PostDcd.DemuxerOutput.V13_AC;
                    EOutPD.EAC(:,:,3)       = PostDcd.DemuxerOutput.V23_AC;

                case  'L2S_LFR-SURV-SWF-E_V02'
                    % Check number of samples/record to see if one can just keep the samples as they are distributed on records.
                    if nSamplesPerRecord ~= 2048
                        error('BICAS:data_manager:Assertion:IllegalArgument', 'Number of samples per CDF record is not 2048.')
                    end
                    
                    EOutPD.Epoch            = PostDcd.Epoch;
                    EOutPD.ACQUISITION_TIME = PostDcd.ACQUISITION_TIME;                    
                    
                    EOutPD.DELTA_PLUS_MINUS = PostDcd.DELTA_PLUS_MINUS;
                    EOutPD.L1_REC_NUM       = PostDcd.L1_REC_NUM;
                    EOutPD.QUALITY_BITMASK  = PostDcd.QUALITY_BITMASK;
                    EOutPD.QUALITY_FLAG     = PostDcd.QUALITY_FLAG;
                    EOutPD.F_SAMPLE         = PostDcd.freqHz;
                    EOutPD.SAMP_DTIME       = PostDcd.SAMP_DTIME;

                    EOutPD.IBIAS1           = PostDcd.IBIAS1;
                    EOutPD.IBIAS2           = PostDcd.IBIAS2;
                    EOutPD.IBIAS3           = PostDcd.IBIAS3;
                    EOutPD.V(:,:,1)         = PostDcd.DemuxerOutput.V1;
                    EOutPD.V(:,:,2)         = PostDcd.DemuxerOutput.V2;
                    EOutPD.V(:,:,3)         = PostDcd.DemuxerOutput.V3;
                    EOutPD.E(:,:,1)         = PostDcd.DemuxerOutput.V12;
                    EOutPD.E(:,:,2)         = PostDcd.DemuxerOutput.V13;
                    EOutPD.E(:,:,3)         = PostDcd.DemuxerOutput.V23;
                    EOutPD.EAC(:,:,1)       = PostDcd.DemuxerOutput.V12_AC;
                    EOutPD.EAC(:,:,2)       = PostDcd.DemuxerOutput.V13_AC;
                    EOutPD.EAC(:,:,3)       = PostDcd.DemuxerOutput.V23_AC;
                    
                otherwise
                    error('BICAS:data_manager:Assertion:IllegalArgument', 'Function can not produce PDID=%s.', eoutPDID)
            end
            
            % ASSERTION            
            bicas.dm_utils.assert_unvaried_N_rows(EOutPD);
        end   % process_PostDCD_to_LFR



        function DemuxerOutput = simple_demultiplex(DemuxerInput, MUX_SET, DIFF_GAIN)
        % Wrapper around "simple_demultiplex_subsequence" to be able to handle multiple CDF records with changing
        % settings (mux_set, diff_gain).
        %
        % NOTE: NOT a processing function (does not derive a PDV).
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % mux_set   = Column vector. Numbers identifying the MUX/DEMUX mode. 
        % input     = Struct with fields BIAS_1 to BIAS_5.
        % diff_gain = Column vector. Gains for differential measurements. 0 = Low gain, 1 = High gain.
        %
        % NOTE: Can handle any arrays of any size as long as the sizes are consistent.
        
            bicas.dm_utils.assert_unvaried_N_rows(DemuxerInput)
            nRecords = length(MUX_SET);
            
            % Create empty structure to which new components can be added.
            DemuxerOutput = struct(...
                'V1',     [], 'V2',     [], 'V3',     [], ...
                'V12',    [], 'V23',    [], 'V13',    [], ...
                'V12_AC', [], 'V23_AC', [], 'V13_AC', []);
            
            iFirst = 1;    % First record in sequence of records with constant settings.
            while iFirst <= nRecords;
                
                % Find continuous sequence of records (i_first to i_last) having identical settings.
                iLast = bicas.dm_utils.find_last_same_sequence(iFirst, DIFF_GAIN, MUX_SET);
                MUX_SET_value   = MUX_SET  (iFirst);
                DIFF_GAIN_value = DIFF_GAIN(iFirst);
                irf.log('n', sprintf('Records %2i-%2i : Demultiplexing; MUX_SET=%-3i; DIFF_GAIN=%-3i', ...
                    iFirst, iLast, MUX_SET_value, DIFF_GAIN_value))    % "%-3" since value might be NaN.
                
                % Extract subsequence of records to "demux".
                demuxerInputSubseq = bicas.dm_utils.select_subset_from_struct(DemuxerInput, iFirst, iLast);
                
                %=================================================
                % CALL DEMUXER - See method/function for comments
                %=================================================
                demuxerOutputSubseq = bicas.data_manager.simple_demultiplex_subsequence(demuxerInputSubseq, MUX_SET_value, DIFF_GAIN_value);
                
                % Add demuxed sequence to the to-be complete set of records.
                DemuxerOutput = bicas.dm_utils.add_components_to_struct(DemuxerOutput, demuxerOutputSubseq);
                
                iFirst = iLast + 1;
                
            end   % while
            
        end   % simple_demultiplex


        
        function Output = simple_demultiplex_subsequence(Input, MUX_SET, DIFF_GAIN)
        % simple_demultiplex_subsequence   Demultiplex, with only constant factors for calibration (no transfer functions).
        %
        % This function implements Table 3 and Table 4 in "RPW-SYS-MEB-BIA-SPC-00001-IRF", iss1rev16.
        % Variable names are chosen according to these tables.
        %
        % NOTE: Conceptually, this function mixes demuxing with calibration which can (and most likely should) be separated.
        % - Demuxing is done on individual samples at a specific point in time.
        % - Calibration (with transfer functions) is made on a time series (presumably of one variable, but could be several).
        %
        % NOTE: This function can only handle one value for mux
        % NOTE: Function is intended for development/testing until there is proper code for using transfer functions.
        % NOTE: "input"/"output" refers to input/output for the function, which is (approximately) the opposite of
        % the physical signals in the BIAS hardware.
        %
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % Input     : Struct with fields BIAS_1 to BIAS_5.
        % MUX_SET   : Scalar number identifying the MUX/DEMUX mode.
        % DIFF_GAIN : Scalar gain for differential measurements. 0 = Low gain, 1 = High gain.
        % Output    : Struct with fields V1, V2, V3,   V12, V13, V23,   V12_AC, V13_AC, V23_AC.
        % 
        % NOTE: Will tolerate values of NaN for MUX_SET_value, DIFF_GAIN_value. The effect is NaN in the corresponding output values.
        %
        % NOTE: Can handle any arrays of any size as long as the sizes are consistent.

        
            %==========================================================================================================
            % QUESTION: How to structure the demuxing?
            % --
            % QUESTION: How split by record? How put together again? How do in a way which
            %           works for real transfer functions? How handle the many non-indexed outputs?
            % QUESTION: How handle changing values of diff_gain, mux_set, bias-dependent calibration offsets?
            % NOTE: LFR data can be either 1 sample/record or 1 snapshot/record.
            % PROPOSAL: Work with some subset of in- and out-values of each type?
            %   PROPOSAL: Work with exactly one value of each type?
            %       CON: Slow.
            %           CON: Only temporary implementation.
            %       PRO: Quick to implement.
            %   PROPOSAL: Work with only some arbitrary subset specified by array of indices.
            %   PROPOSAL: Work with only one row?
            %   PROPOSAL: Work with a continuous sequence of rows/records?
            %   PROPOSAL: Submit all values, and return structure. Only read and set subset specified by indices.
            %
            %
            % PROPOSAL: Could, maybe, be used for demuxing if the caller has already applied the
            % transfer function calibration on on the BIAS signals.
            % PROPOSAL: Validate with some "multiplexer" function?!
            % QUESTION: How handle overdetermined systems one gets when one probe fails?
            % QUESTION: Does it make sense to have BIAS values as cell array? Struct fields?!
            %   PRO: Needed for caller's for loop to split up by record.
            %
            % QUESTION: Is there some better way of implementing than giant switch statement?!
            %
            % PROPOSAL: Only work for scalar values of mux_set and diff_gain?
            % QUESTION: MUX mode 1-3 are overdetermined if we always have BIAS1-3?
            %           If so, how select what to calculate?! What if results disagree/are inconsistent? Check for it?
            %
            % PROPOSAL: Move translation diff_gain-->gamma to separate function (cf. dm_utils.get_LFR_frequency).
            %==========================================================================================================
            
            global CONSTANTS
            
            if numel(MUX_SET) ~= 1 || numel(DIFF_GAIN) ~= 1
                error('BICAS:data_manager:Assertion:IllegalArgument', 'Illegal argument value "mux_set" or "diff_gain". Must be scalars (not arrays).')
            end
            
            ALPHA = CONSTANTS.C.SIMPLE_DEMUXER.ALPHA;
            BETA  = CONSTANTS.C.SIMPLE_DEMUXER.BETA;
            switch(DIFF_GAIN)
                case 0    ; GAMMA = CONSTANTS.C.SIMPLE_DEMUXER.GAMMA_LOW_GAIN;
                case 1    ; GAMMA = CONSTANTS.C.SIMPLE_DEMUXER.GAMMA_HIGH_GAIN;
                otherwise
                    if isnan(DIFF_GAIN)
                        GAMMA = NaN;
                    else
                        error('BICAS:data_manager:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value "diff_gain"=%d.', DIFF_GAIN)                    
                    end
            end
            
            % Set default values which will remain for
            % variables which are not set by the demuxer.
            NaN_VALUES = ones(size(Input.BIAS_1)) * NaN;
            V1_LF     = NaN_VALUES;
            V2_LF     = NaN_VALUES;
            V3_LF     = NaN_VALUES;
            V12_LF    = NaN_VALUES;
            V13_LF    = NaN_VALUES;
            V23_LF    = NaN_VALUES;
            V12_LF_AC = NaN_VALUES;
            V13_LF_AC = NaN_VALUES;
            V23_LF_AC = NaN_VALUES;

            switch(MUX_SET)
                case 0   % "Standard operation" : We have all information.
                    
                    % Summarize what we have;
                    V1_DC  = Input.BIAS_1;
                    V12_DC = Input.BIAS_2;
                    V23_DC = Input.BIAS_3;
                    V12_AC = Input.BIAS_4;
                    V23_AC = Input.BIAS_5;
                    % Convert that which is trivial.
                    V1_LF     = V1_DC / ALPHA;
                    V12_LF    = V12_DC / BETA;
                    V23_LF    = V23_DC / BETA;
                    V12_LF_AC = V12_AC / GAMMA;
                    V23_LF_AC = V23_AC / GAMMA;
                    % Convert that which is less trivial.
                    V13_LF    = V12_LF    + V23_LF;
                    V2_LF     = V1_LF     - V12_LF;
                    V3_LF     = V2_LF     - V23_LF;                    
                    V13_LF_AC = V12_LF_AC + V23_LF_AC;
                    
                case 1   % Probe 1 fails
                    
                    V2_LF     = Input.BIAS_1 / ALPHA;
                    V3_LF     = Input.BIAS_2 / ALPHA;
                    V23_LF    = Input.BIAS_3 / BETA;
                    % input.BIAS_4 unavailable.
                    V23_LF_AC = Input.BIAS_5 / GAMMA;
                    
                case 2   % Probe 2 fails
                    
                    V1_LF     = Input.BIAS_1 / ALPHA;
                    V3_LF     = Input.BIAS_2 / ALPHA;
                    V13_LF    = Input.BIAS_3 / BETA;
                    V13_LF_AC = Input.BIAS_4 / GAMMA;
                    % input.BIAS_5 unavailable.
                    
                case 3   % Probe 3 fails
                    
                    V1_LF     = Input.BIAS_1 / ALPHA;
                    V2_LF     = Input.BIAS_2 / ALPHA;
                    V12_LF    = Input.BIAS_3 / BETA;
                    V12_LF_AC = Input.BIAS_4 / GAMMA;
                    % input.BIAS_5 unavailable.
                    
                case 4   % Calibration mode 0
                    
                    % Summarize what we have;
                    V1_DC  = Input.BIAS_1;
                    V2_DC  = Input.BIAS_2;
                    V3_DC  = Input.BIAS_3;
                    V12_AC = Input.BIAS_4;
                    V23_AC = Input.BIAS_5;
                    % Convert that which is trivial.
                    V1_LF     = V1_DC / ALPHA;
                    V2_LF     = V2_DC / ALPHA;
                    V3_LF     = V3_DC / ALPHA;
                    V12_LF_AC = V12_AC / GAMMA;
                    V23_LF_AC = V23_AC / GAMMA;
                    % Convert that which is less trivial.
                    V12_LF    = V1_LF     - V2_LF;
                    V13_LF    = V1_LF     - V3_LF;
                    V23_LF    = V2_LF     - V3_LF;
                    V13_LF_AC = V12_LF_AC + V23_LF_AC;
                    
                case {5,6,7}
                    error('BICAS:data_manager:Assertion:OperationNotImplemented', 'Not implemented for this value of mux_set yet.')
                    
                otherwise
                    if isnan(MUX_SET)
                        ;   % Do nothing. Allow the default values (NaN) to be returned.
                    else
                        error('BICAS:data_manager:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value for mux_set.')
                    end
            end   % switch
            
            % Create structure to return.
            Output = [];
            Output.V1     = V1_LF;
            Output.V2     = V2_LF;
            Output.V3     = V3_LF;
            Output.V12    = V12_LF;
            Output.V13    = V13_LF;
            Output.V23    = V23_LF;
            Output.V12_AC = V12_LF_AC;
            Output.V13_AC = V13_LF_AC;
            Output.V23_AC = V23_LF_AC;
            
        end  % simple_demultiplex_subsequence
        
        
        
    end   % methods: Static, Access=private

end
