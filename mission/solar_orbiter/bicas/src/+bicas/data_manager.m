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
% (Assuming unique sets of required PDVs. Information flows from left to right, i.e. PDVs on the left are used to derive
% the ones to the right. Each string represents a PDV.)
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
%           .demuxer_input : struct with fields.
%               BIAS_1 to .BIAS_5  : NxM arrays, where M may be 1 (1 sample/record) or >1.
%           .freq               : Snapshot frequency in Hz. Unimportant for one sample/record data.
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
%           .demuxer_output   : struct with fields.
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
%
% PROPOSAL: Better name for "process data" (all together or individual variable), process data type, elementary input/output.
%     PROPOSAL: "process data" (PD), "process data variables" = All process data in aggregate?
%     PROPOSAL: "process data variable" (PDV).
% --
% PROPOSAL: Move deriving DIFF_GAIN (from BIAS HK) from LFR & TDS code separately to combined code.
%       In intermediate PDID?
%
% PROPOSAL: Merge get_processing_info with definitions (constants) for mode and (elementary) input/output datasets?!!!
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
%       Ex: Setting default values for preDCD.QUALITY_FLAG, preDCD.QUALITY_BITMASK, preDCD.DELTA_PLUS_MINUS.
%       Ex: ACQUISITION_TIME.
%   CON: Less assertions can be made in utility functions.
%       Ex: dm_utils.ACQUISITION_TIME_*, dm_utils.tt2000_* functions.
%   CON: ROUNDING ERRORS. Can not be certain that values which are copied, are actually copied.
%   --
%   NOTE: Functions may in principle require integer math to work correctly.
%--
% PROPOSAL: Automatically generate lists of elementary input/output PDIDs, and all PDIDs from
% get_processing_info.
%    CON: Can not generate list of elementary INPUT PDIDs.
%    PROPOSAL: Check these against PDIDs in BICAS constants somehow.
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
        process_data_variables = containers.Map('KeyType', 'char', 'ValueType', 'any');

        INTERMEDIATE_PDIDs = {'PreDCD', 'PostDCD'};
        ALL_PDIDs          = [];  % Initialized by constructor.
        
    end
    
    %###################################################################################################################

    methods(Access=public)
        
        function obj = data_manager()
        % CONSTRUCTOR
        
            global CONSTANTS
        
            obj.ALL_PDIDs = {...
                CONSTANTS.EI_PDIDs{:}, ...
                CONSTANTS.EO_PDIDs{:}, ...
                obj.INTERMEDIATE_PDIDs{:}};
            
            obj.validate
        end

        
        
        function set_elementary_input_process_data(obj, PDID, process_data)
        % Set elementary input process data directly as a struct.
        %
        % NOTE: This method is a useful interface for test code.
        
            global CONSTANTS        
            CONSTANTS.assert_EI_PDID(PDID)            
            obj.set_process_data_variable(PDID, process_data);
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
            irf.log('n', sprintf('Arguments (pdid=%s)', pdid))

            processData = obj.get_process_data_variable(pdid);

            %============================================
            % Check if process data is already available
            %============================================
            if ~isempty(processData)
                % CASE: Process data is already available.
                return   % NOTE: This provides caching so that the same process data are not derived twice.
            end

            % CASE: Process data is NOT already available.

            %=============================================================
            % Obtain processing function, and the PDs which are necessary
            %=============================================================
            [inputPdidsStruct, processingFunc] = obj.get_processing_info(pdid);
            % ASSERTION
            if isempty(processingFunc)
                processData = [];
                return
            end

            %==============================================
            % Obtain the input process datas - RECURSIVELY
            %==============================================
            inputs = struct();
            inputFieldsList = fieldnames(inputPdidsStruct);
            
            for iField = 1:length(inputFieldsList)
                inputField = inputFieldsList{iField};
                pdidList   = inputPdidsStruct.(inputField);
                
                for iPdid = 1:length(pdidList)
                
                    % NOTE: Should return error if can not derive data.
                    inputPd = obj.get_process_data_recursively(pdidList{iPdid});    % NOTE: RECURSIVE CALL.
                    
                    if ~isempty(inputPd)
                        inputs.(inputField).PD   = inputPd;
                        inputs.(inputField).PDID = pdidList{iPdid};
                        break;
                    end
                end
                
                if ~isfield(inputs, inputField)
                    error('BICAS:data_manager:Assertion:SWModeProcessing', ...
                        'Can not derive necessary process data for pdid=%s, input field=%s', pdid, inputField)
                end
                    
            end
            
            %========================================================
            % Derive the actual process data from other process data
            %========================================================
            processData = processingFunc(inputs);
            
            obj.set_process_data_variable(pdid, processData)
            
            % NOTE: This log message is useful for being able to follow the recursive calls of this function.
            irf.log('n', sprintf('End function (pdid=%s)', pdid))
            
        end
        
        
    
        function C_sw_mode = get_C_sw_mode_full(obj, CLI_parameter)
        % Return ~constants structure for a specific S/W mode as referenced by a CLI parameter.
        % 
        % This structure is automatically put together from other structures (constants) to avoid having to define too
        % many redundant constants.
        %
        % NOTE: This function takes what constants.sw_modes returns and adds info about input and output datasets to it.
        % NOTE: The function takes the CLI_parameter as parameter, not the S/W mode ID!
        % IMPLEMENTATION NOTE: The reason for that this function is located in data_manager instead of constants is the
        % call to "bicas.data_manager.get_elementary_input_PDIDs" and that constants.m should contain no reference to
        % data_manager.m (for initialization reasons at the very least).

            global CONSTANTS

            C_sw_mode = bicas.utils.select_structs(CONSTANTS.sw_modes, 'CLI_parameter', {CLI_parameter});            
            C_sw_mode = C_sw_mode{1};

            % Collect all associated elementary input PDIDs.
            %input_PDIDs = obj.get_elementary_input_PDIDs(C_sw_mode.output_PDIDs, C_sw_mode.ID);

            try
                C_sw_mode.inputs = bicas.utils.select_structs(CONSTANTS.inputs,  'PDID', C_sw_mode.input_PDIDs);
            catch exception
                error('BICAS:Assertion:IllegalConfiguration', ...
                    'Can not identify all input PDIDs associated with S/W mode/CLI parameter "%s".', CLI_parameter)
            end
            try
                C_sw_mode.outputs = bicas.utils.select_structs(CONSTANTS.outputs, 'PDID', C_sw_mode.output_PDIDs);
            catch exception
                error('BICAS:Assertion:IllegalConfiguration', ...
                    'Can not identify all output PDIDs associated with S/W mode/CLI parameter "%s".', CLI_parameter)
            end
        end
        
        
        
    end   % methods: Instance, public
    
    %###################################################################################################################

    methods(Access=private)
        
        function validate(obj)
            
            global CONSTANTS
            
            % Initialize instance variable "process_data" with keys (with corresponding empty values).
            for PDID = obj.ALL_PDIDs
                obj.process_data_variables(PDID{1}) = [];
            end
            
            % Validate
            %bicas.utils.assert_strings_unique(obj.ELEMENTARY_INPUT_PDIDs)
            %bicas.utils.assert_strings_unique(obj.ELEMENTARY_OUTPUT_PDIDs)
            bicas.utils.assert_strings_unique(obj.INTERMEDIATE_PDIDs)
            bicas.utils.assert_strings_unique(obj.ALL_PDIDs)
            
            
            
            %========================
            % Iterate over S/W modes
            %========================
            for i = 1:length(CONSTANTS.sw_modes)
                
                % Get info for S/W mode.
                sw_mode_CLI_parameter = CONSTANTS.sw_modes{i}.CLI_parameter;
                C_sw_mode = obj.get_C_sw_mode_full(sw_mode_CLI_parameter);
                
                % Check unique inputs CLI_parameter.
                inputs_CLI_parameters = cellfun(@(s) ({s.CLI_parameter}), C_sw_mode.inputs);
                bicas.utils.assert_strings_unique(inputs_CLI_parameters)
                
                % Check unique outputs JSON_output_file_identifier.
                outputs_JSON_output_file_identifiers = cellfun(@(s) ({s.JSON_output_file_identifier}), C_sw_mode.outputs);
                bicas.utils.assert_strings_unique(outputs_JSON_output_file_identifiers)
            end
            
        end   % validate
        
        
        
        function assert_PDID(obj, PDID)
            if ~ischar(PDID) || ~ismember(PDID, obj.ALL_PDIDs)
                error('BICAS:data_manager:Assertion', 'There is no such PDID="\%s".', PDID)
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
        function process_data = get_process_data_variable(obj, PDID)
            obj.assert_PDID(PDID)
            
            process_data = obj.process_data_variables(PDID);
        end



        %==================================================================================
        % Set a process data variable.
        %
        % ASSERTS: Process data has NOT been set yet.
        % ASSERTS: Valid PDID.
        %==================================================================================
        function set_process_data_variable(obj, PDID, process_data)
            % ASSERTIONS
            obj.assert_PDID(PDID)
            if ~isempty(obj.process_data_variables(PDID))
                error('BICAS:data_manager:Assertion', 'There is already process data for the specified PDID="%s".', PDID);
            end
            
            obj.process_data_variables(PDID) = process_data;
        end

        
        
        % EXPERIMENTAL
        function [input_PDIDs, processing_func] = get_processing_info(obj, output_PDID)
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
        % PDID        : The PDID (string) for the PD that should be produced.
        % input_PDIDs :
        %     Struct where every field is set to a cell array of PDIDs. In every such cell array, only one of its
        %     PDIDs/PDs is necessary for the processing function. The field names of the fields are "human-readable".
        %     NOTE: Elementary input PDIDs yield an empty struct (no field names).
        % processing_func :
        %     Pointer to a function that can derive process data (for output_PDID) from other process data (for
        %     input_PDIDs).
        %         process_data = processing_func(inputs)
        %         ARGUMENTS:
        %             inputs : A struct 
        %                .<field>     : The same field names as input_PDIDs (returned by get_processing_info).
        %                     .PD     : Process data.
        %                     .PDID   : A single PDID (describing the sister field .PD).
        %     If there is no such function (i.e. output_PDID is an EIn-PDID), then empty.
            
            %global CONSTANTS         
            obj.assert_PDID(output_PDID)
            
            % Assign value used for the case of elementary INPUT PDID (no input PDs <==> no fields).
            input_PDIDs = struct();
            
            

            % Assign value used for the case of there being no real processing function
            % -------------------------------------------------------------------------
            % IMPLEMENTATION NOTE: Choice of value for the case of there being no processing function:
            % 1) Using an empty value ==> The caller can immediately check whether it has received a processing function
            % and then immediately throw an assertion error if that was unexpected.
            % 2) Bad alternative: Use a nonsense function (one that always gives an assertion error, albeit a good error message)
            % ==> The (assertion) error will be triggered first when the attempting to call the nonexisting processing
            % function. ==> Error first very late. ==> Bad alternative.
            processing_func = [];
            


            process_PostDCD_to_specific_EO_PD = @(inputs) (bicas.data_manager.process_PostDCD_to_EO_PD(inputs, output_PDID));
            %use_generic_processing_func = 0;   % Default value;
            switch(output_PDID)
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
                    input_PDIDs.HK_cdf = {'HK_BIA_V01'};
                    input_PDIDs.SCI_cdf = {...
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
                    processing_func = @bicas.data_manager.process_LFR_to_PreDCD;
                    
                case 'PostDCD'
                    input_PDIDs.preDCD = {'PreDCD'};
                    processing_func = @bicas.data_manager.demux_calib;

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
                    input_PDIDs.postDCD = {'PostDCD'};
                    processing_func     = process_PostDCD_to_specific_EO_PD;

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
                        'This function has no implementation for this PDID, "%s".', output_PDID)
            end   % switch

        end   % get_processing_info
       
        
        
        
    end   % methods: Instance, private

    
    
    %###################################################################################################################

    
    
    %methods(Static, Access=private)
    methods(Static, Access=public)

        function preDCD = process_LFR_to_PreDCD(inputs)
        % Processing function. Convert LFR CDF data (PDs) to PreDCD.
        
            % PROBLEM: Hardcoded CDF data types (MATLAB classes).

            SCI = inputs.SCI_cdf.PD;
            HK  = inputs.HK_cdf.PD;

            %===============================================================
            % Define variables LFR_V, LFR_E corresponding to zVars which
            % with different names (but not meaning) in different datasets.
            %===============================================================
            switch(inputs.SCI_cdf.PDID)
                case {  'L2R_LFR-SBM1-CWF_V01', ...
                        'L2R_LFR-SBM2-CWF_V01', ...
                        'L2R_LFR-SURV-CWF_V01', ...
                        'L2R_LFR-SURV-SWF_V01'}
                    LFR_V = SCI.POTENTIAL;
                    LFR_E = SCI.ELECTRICAL;
                case {  'L2R_LFR-SBM1-CWF_V02', ...
                        'L2R_LFR-SBM2-CWF_V02', ...
                        'L2R_LFR-SURV-CWF_V02', ...
                        'L2R_LFR-SURV-SWF_V02'}
                    LFR_V = SCI.V;
                    LFR_E = SCI.E;
                    %LFR_SAMP_DTIME = SCI.SAMP_DTIME;
                otherwise
                    error('BICAS:data_manager:Assertion:SWModeProcessing:ConfigurationBug', ...
                        'Can not handle input_PDIDs.SCI_cdf="%s"', input_PDIDs.SCI_cdf)
            end
            
            %========================================================================================
            % Define a variable LFR_FREQ (for all datasets) with the same meaning as LFR zVar FREQ
            % (which exists only for some LFR datasets).
            %========================================================================================
            N_records = size(LFR_V, 1);
            switch(inputs.SCI_cdf.PDID)
                case {  'L2R_LFR-SBM1-CWF_V01', ...
                        'L2R_LFR-SBM1-CWF_V02'}
                    LFR_FREQ = ones(N_records, 1) * 1;   % Always value "1".
                case {  'L2R_LFR-SBM2-CWF_V01', ...
                        'L2R_LFR-SBM2-CWF_V02'}
                    LFR_FREQ = ones(N_records, 1) * 2;   % Always value "2".
                case {  'L2R_LFR-SURV-CWF_V01', ...
                        'L2R_LFR-SURV-CWF_V02', ...
                        'L2R_LFR-SURV-SWF_V01', ...
                        'L2R_LFR-SURV-SWF_V02'}
                    LFR_FREQ = SCI.FREQ;
                otherwise
                    error('BICAS:data_manager:Assertion:ConfigurationBug', ...
                        'Can not handle inputs.SCI_cdf.PDID="%s"', inputs.SCI_cdf.PDID)
            end
            
            Rx = bicas.dm_utils.get_LFR_Rx( SCI.R0, SCI.R1, SCI.R2, LFR_FREQ );   % NOTE: Function can handles "R3".
            
            preDCD.freq = bicas.dm_utils.get_LFR_frequency( LFR_FREQ );
            
            preDCD.demuxer_input        = [];
            preDCD.demuxer_input.BIAS_1 = LFR_V;
            preDCD.demuxer_input.BIAS_2 = bicas.dm_utils.filter_rows( LFR_E(:,:,1), Rx==1 );
            preDCD.demuxer_input.BIAS_3 = bicas.dm_utils.filter_rows( LFR_E(:,:,2), Rx==1 );
            preDCD.demuxer_input.BIAS_4 = bicas.dm_utils.filter_rows( LFR_E(:,:,1), Rx==0 );
            preDCD.demuxer_input.BIAS_5 = bicas.dm_utils.filter_rows( LFR_E(:,:,2), Rx==0 );
            preDCD.Epoch = SCI.Epoch;
            preDCD.ACQUISITION_TIME = SCI.ACQUISITION_TIME;
            
            % IMPLEMENTATION NOTE: QUALITY_FLAG, QUALITY_BITMASK have been found empty in test data, but should have
            % attribute DEPEND_0 = "Epoch" ==> Should have same number of records as Epoch.
            % Can not save CDF with zVar with zero records (crashes when reading CDF). ==> Better create empty records.
            % Test data: MYSTERIOUS_SIGNAL_1_2016-04-15_Run2__7729147__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-SWF_7729147_CNE_V01.cdf
            if isempty(SCI.QUALITY_FLAG)
                preDCD.QUALITY_FLAG    = cast(  zeros(size(preDCD.Epoch)),  bicas.utils.convert_CDF_type_to_MATLAB_class('CDF_UINT1', 'Only CDF data types')  );
            else
                preDCD.QUALITY_FLAG    = SCI.QUALITY_FLAG;
            end
            if isempty(SCI.QUALITY_BITMASK)
                preDCD.QUALITY_BITMASK = cast(  zeros(size(preDCD.Epoch)),  bicas.utils.convert_CDF_type_to_MATLAB_class('CDF_UINT1', 'Only CDF data types')  );
            else
                preDCD.QUALITY_BITMASK = SCI.QUALITY_BITMASK;
            end
                
            % BUG / TEMP: Set default values.
            preDCD.DELTA_PLUS_MINUS = cast(  zeros(size(preDCD.Epoch)),  bicas.utils.convert_CDF_type_to_MATLAB_class('CDF_INT8',  'Only CDF data types')  );

            

            %=========================================================================================================
            % 1) Convert time to something linear in time that can be used for processing (not storing time to file).
            % 2) Effectively also chooses which time to use for the purpose of processing:
            %    ACQUISITION_TIME or Epoch.
            %=========================================================================================================
            t_HK_ACQUISITION_TIME  = bicas.dm_utils.ACQUISITION_TIME_to_tt2000(  HK.ACQUISITION_TIME );
            t_SCI_ACQUISITION_TIME = bicas.dm_utils.ACQUISITION_TIME_to_tt2000( SCI.ACQUISITION_TIME );
            t_HK_Epoch  = HK.Epoch;
            t_SCI_Epoch = SCI.Epoch;
            t_HK        = t_HK_ACQUISITION_TIME;
            t_SCI       = t_SCI_ACQUISITION_TIME;



            %=========================================================================================================
            % Choose where to get MUX_SET from: LFR-SCI, or BIAS-HK
            % -----------------------------------------------------
            % NOTE: Only obtains one MUX_SET per record ==> Can not change MUX_SET in the middle of a record.
            %=========================================================================================================
            preDCD.MUX_SET = interp1(double(t_HK), double(HK.HK_BIA_MODE_MUX_SET), double(t_SCI), 'nearest', NaN);   % Use BIAS HK.
            %std_data.MUX_SET = LFR_cdf.BIAS_MODE_MUX_SET;    % Use LFR SCI.



            %=========================================================================================================
            % Derive approximate DIFF_GAIN values for from BIAS HK
            %
            % NOTE: Not perfect handling of time when 1 snapshot/record, since one should ideally use time stamps
            % for every LFR _sample_.
            %=========================================================================================================
            preDCD.DIFF_GAIN = interp1(double(t_HK), double(HK.HK_BIA_DIFF_GAIN), double(t_SCI), 'nearest', NaN);


            
            % ASSERTIONS
            bicas.dm_utils.assert_unvaried_N_rows(preDCD);
            bicas.dm_utils.assert_unvaried_N_rows(preDCD.demuxer_input);

            
            
            %====================
            % Misc. log messages
            %====================
            % PROPOSAL: PDID, dataset ID.
            % PROPOSAL: Move to demux_calib.
            %    NOTE: Must have t_* variables.
            bicas.dm_utils.log_unique_values_all('LFR_FREQ',  LFR_FREQ);
            bicas.dm_utils.log_unique_values_all('Rx',        Rx);
            bicas.dm_utils.log_unique_values_all('DIFF_GAIN', preDCD.DIFF_GAIN);
            bicas.dm_utils.log_unique_values_all('MUX_SET',   preDCD.MUX_SET);
            
            % PROPOSAL: Move to reading of CDF, setting EIn PD.
            % PROPOSAL: Move to demux_calib.
            bicas.dm_utils.log_tt2000_interval('HK  ACQUISITION_TIME', t_HK_ACQUISITION_TIME)
            bicas.dm_utils.log_tt2000_interval('SCI ACQUISITION_TIME', t_SCI_ACQUISITION_TIME)
            bicas.dm_utils.log_tt2000_interval('HK  Epoch           ', t_HK_Epoch)
            bicas.dm_utils.log_tt2000_interval('SCI Epoch           ', t_SCI_Epoch)
        end



        function postDCD = demux_calib(inputs)
        % Processing function. Convert PreDCD to PostDCD, i.e. demux and calibrate data.
        
            preDCD = inputs.preDCD.PD;
                    
            % Log messages
            for f = fieldnames(preDCD.demuxer_input)'
                bicas.dm_utils.log_unique_values_summary(f{1}, preDCD.demuxer_input.(f{1}));
            end
            
            % ASSERTION
            bicas.dm_utils.assert_unvaried_N_rows(preDCD);
            
            postDCD = preDCD;
            
            % DEMUX
            postDCD.demuxer_output = bicas.data_manager.simple_demultiplex(...
                preDCD.demuxer_input, preDCD.MUX_SET, preDCD.DIFF_GAIN);
            
            % Log messages
            for f = fieldnames(postDCD.demuxer_output)'
                bicas.dm_utils.log_unique_values_summary(f{1}, postDCD.demuxer_output.(f{1}));
            end
            
            % BUG / TEMP: Set default values since the real values are not available.
            postDCD.IBIAS1 = NaN * zeros(size(postDCD.demuxer_output.V1));
            postDCD.IBIAS2 = NaN * zeros(size(postDCD.demuxer_output.V2));
            postDCD.IBIAS3 = NaN * zeros(size(postDCD.demuxer_output.V3));
        end
        

        
        function EO_PD = process_PostDCD_to_EO_PD(inputs, EO_PDID)
        % Processing function. Convert PostDCD to any one of several similar dataset PDs.
        
            postDCD = inputs.postDCD.PD;
            EO_PD = [];
            
            N_smpls_rec = size(postDCD.demuxer_output.V1, 2);   % Samples per record.
            
            switch(EO_PDID)
                case  {'L2S_LFR-SBM1-CWF-E_V02', ...
                       'L2S_LFR-SBM2-CWF-E_V02'}                        
                    error('BICAS:data_manager:OperationNotImplemented', 'Can not produce this EOut PDID.')
                    
                case  'L2S_LFR-SURV-CWF-E_V02'
                    
                    %=====================================================================
                    % Convert 1 snapshot/record --> 1 sample/record (if not already done)
                    %=====================================================================
                    EO_PD.ACQUISITION_TIME = bicas.dm_utils.ACQUISITION_TIME___expand_to_sequences(...
                        postDCD.ACQUISITION_TIME, ...
                        N_smpls_rec, ...
                        postDCD.freq  );
                    EO_PD.Epoch = bicas.dm_utils.tt2000___expand_to_sequences( ...
                        postDCD.Epoch, ...
                        N_smpls_rec, ...
                        postDCD.freq  );
                    for fn = fieldnames(postDCD.demuxer_output)'
                        postDCD.demuxer_output.(fn{1}) = bicas.dm_utils.reshape_to_1_sample_per_record( ...
                            postDCD.demuxer_output.(fn{1}) );
                    end
                    EO_PD.V   = [postDCD.demuxer_output.V1,     postDCD.demuxer_output.V2,     postDCD.demuxer_output.V3];
                    EO_PD.E   = [postDCD.demuxer_output.V12,    postDCD.demuxer_output.V13,    postDCD.demuxer_output.V23];
                    EO_PD.EAC = [postDCD.demuxer_output.V12_AC, postDCD.demuxer_output.V13_AC, postDCD.demuxer_output.V23_AC];
                    EO_PD.IBIAS1 = bicas.dm_utils.reshape_to_1_sample_per_record( postDCD.IBIAS1 );
                    EO_PD.IBIAS2 = bicas.dm_utils.reshape_to_1_sample_per_record( postDCD.IBIAS2 );
                    EO_PD.IBIAS3 = bicas.dm_utils.reshape_to_1_sample_per_record( postDCD.IBIAS3 );

                case  'L2S_LFR-SURV-SWF-E_V02'
                    % Check number of samples/record to see if one can just keep the samples as they are distributed on records.
                    if N_smpls_rec ~= 2048
                        error('BICAS:data_manager:Assertion:IllegalArgument', 'Number of samples per CDF record is not 2048.')
                    end
                    
                    EO_PD.Epoch            = postDCD.Epoch;
                    EO_PD.ACQUISITION_TIME = postDCD.ACQUISITION_TIME;                    
                    EO_PD.F_SAMPLE         = postDCD.freq;
                    EO_PD.V(:,:,1)   = postDCD.demuxer_output.V1;
                    EO_PD.V(:,:,2)   = postDCD.demuxer_output.V2;
                    EO_PD.V(:,:,3)   = postDCD.demuxer_output.V3;
                    EO_PD.E(:,:,1)   = postDCD.demuxer_output.V12;
                    EO_PD.E(:,:,2)   = postDCD.demuxer_output.V13;
                    EO_PD.E(:,:,3)   = postDCD.demuxer_output.V23;
                    EO_PD.EAC(:,:,1) = postDCD.demuxer_output.V12_AC;
                    EO_PD.EAC(:,:,2) = postDCD.demuxer_output.V13_AC;
                    EO_PD.EAC(:,:,3) = postDCD.demuxer_output.V23_AC;
                    EO_PD.IBIAS1 = postDCD.IBIAS1;
                    EO_PD.IBIAS2 = postDCD.IBIAS2;
                    EO_PD.IBIAS3 = postDCD.IBIAS3;
                    
                case 'L2S_TDS-LFM-CWF-E_V02'
                    error('BICAS:data_manager:OperationNotImplemented', 'Can not produce this EOut PDID.')
                case 'L2S_TDS-LFM-RSWF-E_V02'
                    error('BICAS:data_manager:OperationNotImplemented', 'Can not produce this EOut PDID.')
                otherwise
                    error('BICAS:data_manager:Assertion:IllegalArgument', 'Function can not produce this EOut PDID.')
            end
            
            % BUG: How handle?!!!
            % NOTE: Can not change number of records when just copying information!!!
            EO_PD.QUALITY_FLAG     = postDCD.QUALITY_FLAG;
            EO_PD.QUALITY_BITMASK  = postDCD.QUALITY_BITMASK;
            
            EO_PD.DELTA_PLUS_MINUS = postDCD.DELTA_PLUS_MINUS;

        end   % process_PostDCD_to_EO_PD



        function demuxer_output = simple_demultiplex(demuxer_input, MUX_SET, DIFF_GAIN)
        % Wrapper around "simple_demultiplex_subsequence" to be able to handle multiple CDF records with changing
        % settings (mux_set, diff_gain).
        %
        % ARGUMENTS AND RETURN VALUE
        % ==========================
        % mux_set   = Column vector. Numbers identifying the MUX/DEMUX mode. 
        % input     = Struct with fields BIAS_1 to BIAS_5.
        % diff_gain = Column vector. Gains for differential measurements. 0 = Low gain, 1 = High gain.
        %
        % NOTE: Can handle any arrays of any size as long as the sizes are consistent.
        
            bicas.dm_utils.assert_unvaried_N_rows(demuxer_input)
            N_records = length(MUX_SET);
            
            % Create empty structure to which new components can be added.
            demuxer_output = struct(...
                'V1',     [], 'V2',     [], 'V3',     [], ...
                'V12',    [], 'V23',    [], 'V13',    [], ...
                'V12_AC', [], 'V23_AC', [], 'V13_AC', []);
            
            i_first = 1;    % First record in sequence of records with constant settings.
            while i_first <= N_records;
                
                % Find continuous sequence of records (i_first to i_last) having identical settings.
                i_last = bicas.dm_utils.find_last_same_sequence(i_first, DIFF_GAIN, MUX_SET);
                mux_set   = MUX_SET  (i_first);
                diff_gain = DIFF_GAIN(i_first);
                irf.log('n', sprintf('Records %2i-%2i : Demultiplexing; MUX_SET=%-3i; DIFF_GAIN=%-3i', i_first, i_last, mux_set, diff_gain))  % "%-3" since value might be NaN.
                
                % Extract subsequence of records to "demux".
                demuxer_input_subseq = bicas.dm_utils.select_subset_from_struct(demuxer_input, i_first, i_last);
                
                %=================================================
                % CALL DEMUXER - See method/function for comments
                %=================================================
                demuxer_output_subseq = bicas.data_manager.simple_demultiplex_subsequence(demuxer_input_subseq, mux_set, diff_gain);
                
                % Add demuxed sequence to the to-be complete set of records.
                demuxer_output = bicas.dm_utils.add_components_to_struct(demuxer_output, demuxer_output_subseq);
                
                i_first = i_last + 1;
                
            end   % while
            
        end   % simple_demultiplex


        
        function output = simple_demultiplex_subsequence(input, mux_set, diff_gain)
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
        % input     : Struct with fields BIAS_1 to BIAS_5.
        % mux_set   : Scalar number identifying the MUX/DEMUX mode.
        % diff_gain : Gain for differential measurements. 0 = Low gain, 1 = High gain.
        % output    : Struct with fields V1, V2, V3,   V12, V13, V23,   V12_AC, V13_AC, V23_AC.
        % 
        % NOTE: Will tolerate values of NaN for mux_set, diff_gain. The effect is NaN in the corresponding output values.
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
            
            if numel(mux_set) ~= 1 || numel(diff_gain) ~= 1
                error('BICAS:data_manager:Assertion:IllegalArgument', 'Illegal argument value "mux_set" or "diff_gain". Must be scalars (not arrays).')
            end
            
            ALPHA = CONSTANTS.C.approximate_demuxer.alpha;
            BETA  = CONSTANTS.C.approximate_demuxer.beta;
            switch(diff_gain)
                case 0    ; GAMMA = CONSTANTS.C.approximate_demuxer.gamma_lg;
                case 1    ; GAMMA = CONSTANTS.C.approximate_demuxer.gamma_hg;
                otherwise
                    if isnan(diff_gain)
                        GAMMA = NaN;
                    else
                        error('BICAS:data_manager:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value "diff_gain"=%d.', diff_gain)                    
                    end
            end
            
            % Set default values which will remain for
            % variables which are not set by the demuxer.
            NaN_VALUES = ones(size(input.BIAS_1)) * NaN;
            V1_LF     = NaN_VALUES;
            V2_LF     = NaN_VALUES;
            V3_LF     = NaN_VALUES;
            V12_LF    = NaN_VALUES;
            V13_LF    = NaN_VALUES;
            V23_LF    = NaN_VALUES;
            V12_LF_AC = NaN_VALUES;
            V13_LF_AC = NaN_VALUES;
            V23_LF_AC = NaN_VALUES;

            switch(mux_set)
                case 0   % "Standard operation" : We have all information.
                    
                    % Summarize what we have;
                    V1_DC  = input.BIAS_1;
                    V12_DC = input.BIAS_2;
                    V23_DC = input.BIAS_3;
                    V12_AC = input.BIAS_4;
                    V23_AC = input.BIAS_5;
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
                    
                    V2_LF     = input.BIAS_1 / ALPHA;
                    V3_LF     = input.BIAS_2 / ALPHA;
                    V23_LF    = input.BIAS_3 / BETA;
                    % input.BIAS_4 unavailable.
                    V23_LF_AC = input.BIAS_5 / GAMMA;
                    
                case 2   % Probe 2 fails
                    
                    V1_LF     = input.BIAS_1 / ALPHA;
                    V3_LF     = input.BIAS_2 / ALPHA;
                    V13_LF    = input.BIAS_3 / BETA;
                    V13_LF_AC = input.BIAS_4 / GAMMA;
                    % input.BIAS_5 unavailable.
                    
                case 3   % Probe 3 fails
                    
                    V1_LF     = input.BIAS_1 / ALPHA;
                    V2_LF     = input.BIAS_2 / ALPHA;
                    V12_LF    = input.BIAS_3 / BETA;
                    V12_LF_AC = input.BIAS_4 / GAMMA;
                    % input.BIAS_5 unavailable.
                    
                case 4   % Calibration mode 0
                    
                    % Summarize what we have;
                    V1_DC  = input.BIAS_1;
                    V2_DC  = input.BIAS_2;
                    V3_DC  = input.BIAS_3;
                    V12_AC = input.BIAS_4;
                    V23_AC = input.BIAS_5;
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
                    if isnan(mux_set)
                        ;   % Do nothing. Allow the default values (NaN) to be returned.
                    else
                        error('BICAS:data_manager:Assertion:IllegalArgument:DatasetFormat', 'Illegal argument value for mux_set.')
                    end
            end   % switch
            
            % Create structure to return.
            output = [];
            output.V1     = V1_LF;
            output.V2     = V2_LF;
            output.V3     = V3_LF;
            output.V12    = V12_LF;
            output.V13    = V13_LF;
            output.V23    = V23_LF;
            output.V12_AC = V12_LF_AC;
            output.V13_AC = V13_LF_AC;
            output.V23_AC = V23_LF_AC;
            
        end  % simple_demultiplex_subsequence
        
        
        
    end   % methods: Static, Access=private

end
