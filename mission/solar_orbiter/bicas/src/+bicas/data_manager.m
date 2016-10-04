% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-10
%
% "Data manager". Does the actual processing of datasets. Inspired by irfu-matlab's MMS data manager.
%
%
%
% BASIC CONCEPT
% The user gives the class the input data ("elementary input process data") it has, and asks for the output data
% ("elementary output process data") it wants.
%
% The class maintains an internal list of different "process data" variables, a list of variables where each one is
% uniquely referenced by a unique string, a "process data type". These process data variables are all related to each
% other in a conceptual "web" (an acyclic directed graph) describing their dependencies on each other. Initially these
% process data variables are all empty.
%
% There are three forms of process data:
%   (1) elementary input process data  : Supplied BY the user. In practise this should correspond to the content of a CDF file.
%   (2) elementary output process data : Returned TO the user. In practise this should correspond to the content of a CDF file.
%                                        Derived from a fixed set of other process data.
%   (3) intermediate process data      : Derived from a fixed set of other process data. These only exist inside data_manager.
%
% When a process data variable Y is requested and it is not already derived, the other process data variables X_1, ...,
% X_n required for it are requested. Then Y is derived from the X_1, ..., X_n. If a process data variable X_i is not
% already derived, it will be derived the same way, recursively. If a variable can not be derived, because of error, or
% because it can only be supplied by the user and the user has not, the process fails.
%
% Example: (Information flows from left to right. Each string represents a process data variable.)
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
% NOTE: Advantages with architecture.
% -- Easy to implement _shifting_ S/W modes (as defined by the RCS ICD) using this class, althoguh the class itself is
%    unaware of S/W modes.
%      Ex: Updates to the RPW pipeline, datasets.
%      Ex: Can use a different set of modes when validating using the official RCS validation software at LESIA (not all
%      datasets are available for input there, and one gets error when the S/W descriptor describes modes which require
%      non-existing dataset IDs).
%      Ex: Implement inofficial modes for testing(?).
% -- Easier to keep support for legacy datasets (older versions of CDFs)?
% -- Easier to maintain a varying set of test/debugging modes?
% -- Indirect "caching".
% -- Can reuse processing code?
% -- Other code can automatically derive which datasets are required to derive which. ==> Can automatically derive (1) the
%    S/W modes in the S/W descriptor, (2) the required CLI parameters.
%
%
%
% DEFINITIONS OF TERMS
% -- process data
%        In practise, one internal instance variable representing some form of data in some step of
%        the "data processing process".
% -- elementary input/output process data
%        Process data that represents data that goes in and out of the data manager as a whole. The word "elementary" is
%        added to distinguish input/output from the input/output for a smaller operation.
% -- intermediate process data
%        Process data that is not elementary.
% -- process data type
%        A string that uniquely represents (refers to) a type of process data.
%        By convention, the process data types for elementary input/output process data types are a shortened version of
%        dataset_ID+skeleton_version_str, e.g. "L2S_LFR-SURV-CWF-E_V01". Hardcoded values of these occur throughout the
%        code as constants.
%
% CONVENTIONS: 
% -- It is implicit that arrays representing CDF data, or "CDF-like" data, use the first MATLAB array index to represent
%    CDF records.
%
% IMPLEMENTATION NOTE: This class is intended to:
% - Not be aware of S/W modes as defined in the RCS ICD.
% - Not deal with writing CDFs. It does however read CDF files as a way of setting elementary input process data since
%   it is trivial.
%
classdef data_manager
%#######################################################################################################################
% PROPOSAL: Use other class name that implices processing, and fits with defined term "process data".
%     "processing_manager"? "proc_manager"?
%
% PROPOSAL: Better name for "process data" (all together or individual variable), process data type, elementary input/output.
%     PROPOSAL: "process data", "process data variables" = All process data in aggregate?
%     PROPOSAL: "process data variable" (PDV)
%     PROPOSAL: "process data ID" (PDID).
%
% PROPOSAL: Different name for process data variable: data, datas, data_types, internal_data,
%     data_network, data_chain, chain_data, process_data, processing_data.
%     A one-letter name? Shortening? PD and PDT?
%
% QUESTION: How handle master CDFs?
%   NOTE: Code needs the fill/pad values for the final CDF format.
%         (1) Can not work entirely with NaN internally since some MATLAB classes (non-floating point) do not have NaN.
%         (2) In principle, there is only one NaN which can not both represent the pad value and the fill value.
%   NOTE: There is ALWAYS exactly one master CDF per elementary output.
%   NOTE: Code might want to use the master CDF zVariable sizes as input.
%   --
%   PROPOSAL: No master CDFs (no data, no pad/fill values) inside data_manager. Convert NaN to pad or fill value outside
%             before writing to file.
%       CON: Some MATLAB classes (non-floating point) type do not have NaN.
%   PROPOSAL: Supply the master CDF with the call to get_process_data_recursively.
%       CON: Won't work if it in turn requires other process data corresponding to another (output) CDF which requires
%            its own master CDF.
%   PROPOSAL: User gives the master CDF to the data_manager which stores it in instance variable (containers.Map for output process data).
%   PROPOSAL: data_manager has access to CONSTANTS/function which retrieves master CDFs.
%
%
% QUESTION: Should this class AT ALL handle reading and writing CDF _FILES_? Should that not be done in some separate layer?!
%    NOTE: Concerns input, output, and MASTER CDF files.
%    NOTE: Need to handle conversion to/from fill/pad values which more closely tied to the CDF files.
%    PROPOSAL: Only work with input/output data types that correspond tightly to CDF files.
%       CON: Would work against the idea of combining code with the S/W descriptor info.
%          CON: Not really. How?
%       CON: Might need access to fill/pad values from the CDF templates for the final CDF files.
%          CON: Could make conversions NaN-->pad/fill value outside of data_manager.
%
% PROPOSAL: Implement master CDFs as internal process data types?!!
%    PRO: Does actually need the fill/pad values.
%    CON: Elementary inputs should fit with S/W description and CLI parameters. Master CDFs will mess that up, or BICAS
%         will have to distinguish between different types of elementary input process data.
%    CON: Can just as well have the data_manager read master files itself.
%
% QUESTION: Is dataobj a handle class?! How (deep) copy? Needed?! If so, then matters for setting elementary input
% processing data directly.
%
% PROPOSAL: Functions assert_is_elem_input_PDT, assert_is_elem_output_PDT, assert_is_PT.
%
% PROPOSAL: Automatically generate lists of elementary input/output process data types, and all process data types from
% get_processing_info.
%    CON: Can not generate list of elementary INPUT process data types.
%    PROPOSAL: Check these against process data types in BICAS constants somehow.
%
% PROPOSAL: Merge get_processing_info with definitions (constants) for mode and (elementary) input/output datasets?!!!
%
%
%
% TODO: General functions for handling CDFs:
%    -validating input/output CDFs: dataset IDs, (same variable sizes?, same nbr of
%     records), obtain number of records(?!!).
%    -setting Parents+Parent_version,
%    -reading, finding, validating/checking (Skeleton_version) master CDF
%
% TODO: General functions for handling data.
%    -Convert fill/pad values <---> NaN.
%
% ~BUG: Documentation (beginning of file) does not mention dependencies varying due to sw_modes.
%
% NOTE: Both BIAS HK and LFR SURV CWF contain MUX data (only LFR has one timestamp per snapshot). True also for other input datasets?
%#######################################################################################################################

    properties(Access=private)
        % NOTE: The code does not need any instance variable "ALL_PROCESS_DATA_TYPES" since the keys in "process_data"
        % (initialized in the constructor) have the same function.
        
        
        
        % NOTE: Only set the map's keys here. This makes sure that no disallowed keys are ever used.
        % For data corresponding to CDFs: Use dataset IDs.
        % Initialized by constructor.
        process_data = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        ELEMENTARY_INPUT_PROCESS_DATA_TYPES  = {'L2R_LFR-SURV-CWF_V01',   'L2R_LFR-SURV-SWF_V01',  'HK_BIA_V01', 'L2R_TEST_V99'};
        ELEMENTARY_OUTPUT_PROCESS_DATA_TYPES = {'L2S_LFR-SURV-CWF-E_V01', 'L2S_LFR-SURV-SWF-E_V01',              'L2S_TEST_V99'};
        INTERMEDIATE_PROCESS_DATA_TYPES      = {};
        
    end
    
    %###################################################################################################################

    methods(Access=public)
        
        %=============
        % Constructor
        %=============
        function obj = data_manager()
            
            % Initialize instance variable "ALL_PROCESS_DATA_TYPES".
            ALL_PROCESS_DATA_TYPES = {...
                obj.ELEMENTARY_INPUT_PROCESS_DATA_TYPES{:}, ...
                obj.ELEMENTARY_OUTPUT_PROCESS_DATA_TYPES{:}, ...
                obj.INTERMEDIATE_PROCESS_DATA_TYPES{:}};
            
            % Initialize instance variable "process_data".
            for process_data_type = ALL_PROCESS_DATA_TYPES
                obj.process_data(process_data_type{1}) = [];
            end
            
            % Validate
            assert_strings_unique(obj.ELEMENTARY_INPUT_PROCESS_DATA_TYPES)
            assert_strings_unique(obj.ELEMENTARY_OUTPUT_PROCESS_DATA_TYPES)
            assert_strings_unique(obj.INTERMEDIATE_PROCESS_DATA_TYPES)
            assert_strings_unique(    ALL_PROCESS_DATA_TYPES)
        end

        
        
        %==================================================================================
        % Set elementary input process data via CDF file.
        % NOTE: Not preventing from setting data types/properties which are not CDF files.
        %==================================================================================
        function set_elementary_input_CDF(obj, process_data_type, file_path)
            % PROPOSAL: "Validate" the read CDF file here.
            
            irf.log('n', sprintf('process_data_type=%s: file_path=%s', process_data_type, file_path))    % NOTE: irf.log adds the method name.
            global ERROR_CODES
            
            % Assertion
            if ~ismember(process_data_type, obj.ELEMENTARY_INPUT_PROCESS_DATA_TYPES)
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Data type is not an elementary input data type.')
            end
            
            obj.set_process_data_variable(process_data_type, dataobj(file_path));
        end



        %=============================================================================================================
        % Get process data. If it is already available (stored), return it, if not, derive it from other process data
        % recursively.
        %
        % Should (indirectly) return error if can not return process data.
        %=============================================================================================================
        function process_data = get_process_data_recursively(obj, process_data_type, sw_mode_ID)
            
            global ERROR_CODES
            
            %==================================================================================
            % If process data already exists - Return it and exit
            % 
            % NOTE: This provides caching so that the same process data are not derived twice.
            %==================================================================================
            process_data = obj.get_process_data_variable(process_data_type);
            if ~isempty(process_data)
                return   % Return with already available process data.
            end
            
            [input_process_data_types, processing_func] = bicas.data_manager.get_processing_info(process_data_type, sw_mode_ID);            
            % Assertion
            if isempty(processing_func)
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Received no processing function necessary for deriving process data (type "%s").', process_data_type)
            end
            
            %==============================================
            % Obtain the input process datas - RECURSIVELY
            %==============================================
            input_process_datas = struct();
            fn_list = fieldnames(input_process_data_types);
            for i = 1:length(fn_list)
                fn = fn_list{i};
                
                % NOTE: Should return error if can not derive data.
                input_process_datas.(fn) = obj.get_process_data_recursively(...
                    input_process_data_types.(fn), ...
                    sw_mode_ID);                         % NOTE: RECURSIVE CALL.
            end
            
            %=============================================
            % Derive process data from other process data
            %=============================================
            process_data = processing_func(input_process_datas);

        end
    
    end   % methods: Instance, public
    
    %###################################################################################################################

    methods(Access=private)
        
        %==================================================================================
        % Return process data value.
        % Does NOT try to fill it with (derived) data if empty.
        % Can be used to find out whether process data has been set.
        %
        % ASSERTS: Valid process_data_type.
        % DOES NOT ASSERT: Process data has been set.
        %==================================================================================
        function process_data = get_process_data_variable(obj, process_data_type)
            global ERROR_CODES
            
            if ~obj.process_data.isKey(process_data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is no such process data type, "%s".', process_data_type);
            end
            process_data = obj.process_data(process_data_type);
        end



        %==================================================================================
        % Set a process data variable.
        %
        % ASSERTS: Process data has NOT been set yet.
        % ASSERTS: Valid process_data_type.
        %==================================================================================
        function set_process_data_variable(obj, process_data_type, process_data)
            global ERROR_CODES
            
            if ~obj.process_data.isKey(process_data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is no such process data type, "%s".', process_data_type);
            elseif ~isempty(obj.process_data(process_data_type))
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is already process data for the specified data type "%s".', process_data_type);
            end
            obj.process_data(process_data_type) = process_data;
        end
        
    end   % methods: Instance, private
    
    %###################################################################################################################
    
    methods(Static, Access=public)
    
        %===============================================================================================================
        % For every process data type (combined with an appropriate S/W mode), return meta-information needed to
        % derive the process data.
        %
        % This function effectively
        % 1) (IMPORTANT) defines/describes the dependencies between different process data types (an acyclic directed
        %    graph, albeit dependent on S/W modes), and
        % 2) (LESS IMPORTANT) implicitly decides how to informally distinguish the input process datas by giving them
        %    suitable struct field names which the processing functions recognize them with.
        %
        % The method itself is NOT recursive but it is designed so that other code can recurse over
        % the implicit "acyclic graph" defined by it.
        % It can be used to recursively:
        % (1) derive all the necessary elementary input process data types needed to produce an elementary output
        %     process data type.
        %     This is useful for automatically generating the required input datasets for S/W modes (for the S/W
        %     descriptor, for required CLI parameters).
        % (2) derive process data recursively.
        %
        % input_process_data_types :
        %     Struct with fields set to the necessary process data types.
        %     Elementary input process data types yield empty structs.
        % processing_func :
        %     Pointer to function that can derive the process data.
        %         process_data = processing_func(input_process_datas)
        %     The function accepts exactly one argument: a struct analogous to input_process_data_types but with fields set to
        %     the corresponding process data instead.
        % 
        % NOTE: A reasonable implementation does not check for the validity of sw_mode_ID in all cases.
        %===============================================================================================================
        function [input_process_data_types, processing_func] = get_processing_info(process_data_type, sw_mode_ID)
            % PROPOSAL: Warning/error for not assigning a processing function.
            % PROPOSAL: Should check process_data_type in addition to switch statement.
            % PROPOSAL: Reverse the switch statements for S/W mode and process data types. S/W mode outermost.
            %    CON: Would in total give more checks (switch), longer code.(?) Can not omit the check for cases with
            %    only one mode.
            % PROBLEM: Do not want to include to much dependence on modes.
            
            global ERROR_CODES
            
            % Default values. These are returned to the caller if not amended to or overwritten.
            % NOTE: "func" stores the variable values which are not its argument (process_data_type, sw_mode_ID).
            inputs = struct();
            func   = @(input_process_datas) (...
                errorp(...
                    ERROR_CODES.SW_MODE_PROCESSING_ERROR, ...
                    'The processing function for process data type "%s", S/W mode "%s" has not yet been implemented.', ...
                    process_data_type, sw_mode_ID) ...
                );

            
            
            switch(process_data_type)
                %=====================================
                % Elementary INPUT process data types
                %=====================================
                case 'HK_BIA_V01'

                % LFR
                case 'L2R_LFR-SBM1-CWF_V01'
                case 'L2R_LFR-SBM2-CWF_V01'
                case 'L2R_LFR-SURV-CWF_V01'
                case 'L2R_LFR-SURV-SWF_V01'
                    
                % TDS
                case 'L2R_TDS-LFM-CWF_V01'                    
                case 'L2R_TDS-LFM-RSWF_V01'
                case 'L2R_TDS-LFM-RSWF_V02'

                % TEST
                %case 'L2R_TEST_V99';

                %======================================
                % Elementary OUTPUT process data types
                %======================================
                
                %-----
                % LFR
                %-----                
                case 'L2S_LFR-SBM1-CWF-E_V01'
                case 'L2S_LFR-SBM2-CWF-E_V01'
                case 'L2S_LFR-SURV-CWF-E_V01'
                    inputs.HK_cdf  = 'HK_BIA_V01';
                    inputs.SCI_cdf = 'L2R_LFR-SURV-CWF_V01';
                    func = @bicas.data_manager.process_LFR_CWF_CDF_to_BIAS_CDF;
                case 'L2S_LFR-SURV-SWF-E_V01'
                    inputs.HK_cdf  = 'HK_BIA_V01';
                    inputs.SCI_cdf = 'L2R_LFR-SURV-SWF_V01';

                %-----
                % TDS
                %-----
                case 'L2S_TDS-LFM-CWF-E_V01'
                   inputs.HK_cdf  = 'HK_BIA_V01';
                   inputs.SCI_cdf = 'L2R_TDS-LFM-CWF_V01';
                            
                case 'L2S_TDS-LFM-RSWF-E_V01'
                    inputs.HK_cdf  = 'HK_BIA_V01';
                    switch(sw_mode_ID)
                        case 'TDS-LFM-RSWF-E_V01-V01' ; inputs.SCI_cdf = 'L2R_TDS-LFM-RSWF_V01';
                        case 'TDS-LFM-RSWF-E_V02-V01' ; inputs.SCI_cdf = 'L2R_TDS-LFM-RSWF_V02';
                        otherwise ; error_bad_sw_mode(process_data_type, sw_mode_ID)
                    end

                %---------------------------------
                % Intermediate process data types
                %---------------------------------

                case 'input_waveform_data_std_format'
                    % Represents the input waveform data, without any HK (incl. LFR's BIAS_MODE_MUX_SET).
                    % Name? raw?
                    switch(sw_mode_ID)
                        case 'TDS-LFM-CWF-E_V01-V01'  ; inputs.SCI_cdf = 'L2R_TDS-LFM-CWF_V01';
                        case 'TDS-LFM-RSWF-E_V01-V01' ; inputs.SCI_cdf = 'L2R_TDS-LFM-RSWF_V01';
                        case 'TDS-LFM-RSWF-E_V02-V01' ; inputs.SCI_cdf = 'L2R_TDS-LFM-RSWF_V02';
                        otherwise ; error_bad_sw_mode(process_data_type, sw_mode_ID)
                    end
                
                %case 'HK_std_format'
                %    % Could also take information from LFR since their datasets also contain (some) HK in BIAS_MODE_MUX_SET.
                %    inputs.HK = 'HK_BIA_V01';
                
                %------
                % TEST
                %------
                
                %case 'L2S_TEST_V99'
                %    inputs.SCI2  = 'L2R_TEST_V99';
                    
                %============================================
                % ERROR: Can not recognize process data type
                %============================================                    
                otherwise
                    errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Can not produce this process data type, "%s".', process_data_type)
            end   % switch
            
            
            
            % Change name before returning values.
            input_process_data_types = inputs;
            processing_func = func;
            
            
            
            %-----------------------------------------------------------------------------------------------------------
            function error_bad_sw_mode(process_data_type, sw_mode_ID)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Can not interpret S/W mode ID (%s) for this process data type (%s).', sw_mode_ID, process_data_type)
            end
            %-----------------------------------------------------------------------------------------------------------
        end
        
        
        
        %===============================================================================================================
        % Collect all the elementary input process data types needed to produce a given process data type.
        %
        % Can be used to determine which process data types are needed for a S/W mode. Recursive.
        %
        % elementary_input_process_data_types :
        %   Cell list of strings. Never contains any duplicates. If process_data_type is an elementary data type, then it
        %   returns a cell array containing only "process_data_type" itself.
        %===============================================================================================================
        function elementary_input_process_data_types = get_elementary_input_process_data_types(process_data_type, sw_mode_ID)
            % PROPOSAL: Change implementation. Iterate over list instead.
            %     PRO: Better way of handling duplicates.
            % T = {start_type}
            % %  T_handled = {}  // Use?!!
            % T_f = {}
            % while (T not empty)
            %     T_next = {}
            %     for t in T
            %         if t elementary input
            %             Put t in T_f
            %         else
            %             Add t's input types in T_next
            %         end
            %     end
            %     Remove doubles from T_next
            %     T = T_next
            % end
            
            [input_process_data_types, ~] = bicas.data_manager.get_processing_info(process_data_type, sw_mode_ID);
            
            elem_input_types = {};
            
            if isempty(fieldnames(input_process_data_types))
                % CASE: process_data_type is an ELEMENTARY INPUT process data type.
                elem_input_types = {process_data_type};
            else
                % CASE: process_data_type is NOT an NON-ELEMENTARY INPUT process data type.
                fn_list = fieldnames(input_process_data_types);
                for i = 1:length(fn_list)
                    fn = fn_list{i};
                
                    new_types = bicas.data_manager.get_elementary_input_process_data_types(...
                        input_process_data_types.(fn), ...
                        sw_mode_ID);                        % NOTE: RECURSIVE CALL
                    
                    elem_input_types = {elem_input_types{:}, new_types{:}};
                end
            end
            
            elementary_input_process_data_types = unique(elem_input_types);
        end
        
        
        
    end  % methods: Static, public
    
    %###################################################################################################################
        
    methods(Static, Access=private)

        function output = process_LFR_CWF_CDF_to_BIAS_CDF(input_process_datas)
            % INCOMPLETE/INCORRECT IMPLEMENTATION
            % Does not take LFR's R0/R1/R2 into account. Current implementation probably based on misinterpretation of how the dataset format works.
            %
            % PROPOSAL: Move assignment of base variables into record loop? Its own record loop?
            % PROPOSAL: Convert everything to/from 1 sample/record and implement core algorithms as 1 sample/record.
            %    NOTE: R0, R1.
            %    NOTE: Convert data from HK time stamps to SAMPLE time stamps (not snapshot or record time stamps).
            %    CON: Conversion back to waveforms, assuming that some values are constant over entire snapshot?
            %
            % Might be relevant:
            %    True/false depending on available BIAS_1/2/3? Should be in
            %    agreement with MUX_SET?
            %    input_BIAS_HK.BIAS_MODE_BIAS1_ENABLED; ...
            %    input_LFR.data.R0; ...

            BIAS_HK_cdf = input_process_datas.HK_cdf.data;
            LFR_cdf     = input_process_datas.SCI_cdf.data;

            N_LFR_records          = size(LFR_cdf.POTENTIAL.data, 1);
            N_samples_per_snapshot = size(LFR_cdf.POTENTIAL.data, 2);
            
            % Change to standard names.
            %R0 = LFR_cdf.R0.data;
            %R1 = LFR_cdf.R1.data;
            %R2 = LFR_cdf.R2.data;
            demuxer_input = [];
            demuxer_input.BIAS_1 = LFR_cdf.POTENTIAL.data;
            demuxer_input.BIAS_2 = LFR_cdf.ELECTRICAL.data(:,:,1);                 % BUG: Can not assume that these values should be used. Depends on R0/R1/R2.
            demuxer_input.BIAS_3 = LFR_cdf.ELECTRICAL.data(:,:,2);                 % BUG: Can not assume that these values should be used. Depends on R0/R1/R2.
            demuxer_input.BIAS_4 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;    % No data.
            demuxer_input.BIAS_5 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;    % No data.
            %================================================================================
            % Find sequences of records (1 snapshot/record) with constant settings incl. R0/R1/R2. - UNFINISHED; BASED ON FALSE PREMISES
            %================================================================================
            %                     i_first = 1;    % First record in sequence of records with constant settings.
            %                     while i_first <= N_LFR_records;
            %
            %                         % Find sequence of records having identical settings.
            %                         i_last = bicas.data_manager.find_last_same_sequence(...
            %                             i_first, ...
            %                             LFR_cdf.R0.data, ...
            %                             LFR_cdf.R1.data, ...
            %                             LFR_cdf.R2.data);   % BUG: Only one R0, R1, R2 needs to be checked for.
            %
            %                     	R0 = LFR_cdf.R0.data(i_first);
            %                     	R1 = LFR_cdf.R1.data(i_first);
            %                     	R1 = LFR_cdf.R1.data(i_first);
            %
            %                         demuxer_input.BIAS_2 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;    % No data.
            %                         demuxer_input.BIAS_3 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;    % No data.
            %                         demuxer_input.BIAS_4 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;    % No data.
            %                         demuxer_input.BIAS_5 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;    % No data.
            %                         if R0==1 && R1==0
            %                             demuxer_input.BIAS_2 = LFR_cdf.ELECTRICAL.data(i_first:i_last,:,1);
            %                             demuxer_input.BIAS_3 = LFR_cdf.ELECTRICAL.data(i_first:i_last,:,2);
            %                         elseif R0==0 && R1==1
            %                             demuxer_input.BIAS_4 = LFR_cdf.ELECTRICAL.data(i_first:i_last,:,1);
            %                             demuxer_input.BIAS_5 = LFR_cdf.ELECTRICAL.data(i_first:i_last,:,2);
            %                         elseif R0==1 && R1==1
            %                             errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR ,'Illegal combination of LFR CDF R0=1 and R1=1.')
            %                         end
            %
            %                         i_first = i_last + 1;
            %
            %                     end   % while
            
            %------------------------------------------------------------------------------------------------
            % Derives approximate DIFF_GAIN values for LFR SCI time stamps/records (one per snapshot), using
            % DIFF_GAIN values for BIAS HK time stamps.
            % NOTE: Not perfect, since one should ideally use time stamps for every LFR _sample_.
            %------------------------------------------------------------------------------------------------
            DIFF_GAIN = bicas.data_manager.nearest_interpolate_records(...
                BIAS_HK_cdf.ACQUISITION_TIME.data, ...
                BIAS_HK_cdf.HK_BIA_DIFF_GAIN.data, ...
                LFR_cdf    .ACQUISITION_TIME.data);
            
            %------------------------------------------------------------------------------------------------------
            % Choose where to get MUX_SET from
            % --------------------------------
            % NOTE: One can obtain MUX_SET from either (1) LFR, or (2) BIAS HK.
            % NOTE: Not perfect handling of time, since one should ideally use time stamps for every LFR _sample_.
            %------------------------------------------------------------------------------------------------------
            MUX_SET = bicas.data_manager.nearest_interpolate_records(...
                BIAS_HK_cdf.ACQUISITION_TIME.data, ...
                BIAS_HK_cdf.HK_BIA_MODE_MUX_SET.data, ...
                LFR_cdf    .ACQUISITION_TIME.data);      % Use BIAS HK.
            %MUX_SET = LFR_cdf.BIAS_MODE_MUX_SET.data;    % Use LFR SCI.
            
            
            
            %================================================================================
            % Run demuxer on sequences of records (1 snapshot/record) with constant settings
            %================================================================================
            i_first = 1;    % First record in sequence of records with constant settings.
            demuxer_output = struct(...
                'V1',     [], 'V2',     [], 'V3',     [], ...
                'V12',    [], 'V23',    [], 'V13',    [], ...
                'V12_AC', [], 'V23_AC', [], 'V13_AC', []);
            while i_first <= N_LFR_records;
                
                % Find sequence of records having identical settings.
                i_last = bicas.data_manager.find_last_same_sequence(...
                    i_first, ...
                    DIFF_GAIN, ...
                    MUX_SET, ...
                    LFR_cdf.FREQ.data);
                
                diff_gain        = DIFF_GAIN(i_first);
                mux_set          = MUX_SET(i_first);
                sample_frequency = bicas.data_manager.get_LFR_samples_in_snapshot_frequency(LFR_cdf.FREQ.data(i_first));
                
                demuxer_input_subseq = bicas.data_manager.select_subset_from_struct(demuxer_input, i_first, i_last);
                
                %=================================================
                % CALL DEMUXER - See method/function for comments
                %=================================================
                demuxer_output_subseq = bicas.data_manager.simple_demultiplex(demuxer_input_subseq, mux_set, diff_gain);
                
                demuxer_output = bicas.data_manager.add_components_to_struct(demuxer_output, demuxer_output_subseq);
                
                i_first = i_last + 1;
                
            end   % while
            
            
            
            %===============================================
            % Convert 1 snapshot/record --> 1 sample/record
            %===============================================
            ACQUISITION_TIME = bicas.data_manager.convert_snapshotsPR_to_samplesPR_ACQUISITION_TIME(...
                LFR_cdf.ACQUISITION_TIME.data, ...
                N_samples_per_snapshot, ...
                sample_frequency);
            Epoch = bicas.data_manager.convert_snapshotsPR_to_samplesPR_Epoch_TEMP(  LFR_cdf.Epoch.data, N_samples_per_snapshot  );
            fn_list = fieldnames(demuxer_output);
            for i = 1:length(fn_list)
                fn = fn_list{i};
                demuxer_output.(fn) = bicas.data_manager.convert_snapshotPR_to_samplePR_DATA(  demuxer_output.(fn), N_samples_per_snapshot  );
            end
            
            
            
            %==============================================================
            % Put together variables to be written to CDF file (elsewhere)
            %==============================================================
            output = [];
            output.Epoch = Epoch;
            output.ACQUISITION_TIME = ACQUISITION_TIME;
            output.V   = [demuxer_output.V1,     demuxer_output.V2,     demuxer_output.V3];
            output.E   = [demuxer_output.V12,    demuxer_output.V13,    demuxer_output.V23];
            output.EAC = [demuxer_output.V12_AC, demuxer_output.V13_AC, demuxer_output.V23_AC];
            output.QUALITY_FLAG     = cast(  zeros(size(Epoch)),  convert_CDF_type_to_MATLAB_class('CDF_UINT1', 'Only CDF data types')  );
            output.QUALITY_BITMASK  = cast(  zeros(size(Epoch)),  convert_CDF_type_to_MATLAB_class('CDF_UINT1', 'Only CDF data types')  );
            output.DELTA_PLUS_MINUS = cast(  zeros(size(Epoch)),  convert_CDF_type_to_MATLAB_class('CDF_INT8',  'Only CDF data types')  );
            
            output.IBIAS1 = NaN * zeros(size(demuxer_output.V1));
            output.IBIAS2 = NaN * zeros(size(demuxer_output.V2));
            output.IBIAS3 = NaN * zeros(size(demuxer_output.V3));
        end



        % Generic utility function.
        % Given a struct, select a subset of that struct defined by a range of column indicies for every field.
        function s = select_subset_from_struct(s, i_first, i_last)
            fn_list = fieldnames(s);
            for i=1:length(fn_list)
                fn = fn_list{i};
                
                s.(fn) = s.(fn)(i_first:i_last, :, :);
            end
        end
        
        

        % Generic utility function.
        % Add values to every struct field by adding components after their highest column index (let them grow in
        % the column index).
        function s = add_components_to_struct(s, s_amendment)
            fn_list = fieldnames(s_amendment);
            for i=1:length(fn_list)
                fn = fn_list{i};
                
                s.(fn) = [s.(fn) ; s_amendment.(fn)];
            end
        end



        function freq = get_LFR_samples_in_snapshot_frequency(LFR_FREQ)
            global CONSTANTS ERROR_CODES            
            switch(LFR_FREQ)
                case 0
                    freq = CONSTANTS.C.LFR.F0;
                case 1
                    freq = CONSTANTS.C.LFR.F1;
                case 2
                    freq = CONSTANTS.C.LFR.F2;
                case 3
                    freq = CONSTANTS.C.LFR.F3;
                otherwise
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal LFR_FREQ value.')
            end
        end



        %=====================================================================================================================
        % Finds the greatest i_last such that all varargin{k}(i) are equal for i_first <= i <= i_last separately for every k.
        % Useful for finding a continuous sequence of records with the same data.
        %
        % ASSUMES: varargin{i} are all column arrays of the same size.
        %=====================================================================================================================
        function i_last = find_last_same_sequence(i_first, varargin)
            % PROPOSAL: Better name?
            
            global ERROR_CODES
            
            % Check arguments
            for k = 1:length(varargin)
                if ~iscolumn(varargin{k})
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'varargins are not all column vectors.')
                end
            end
                
            N_records = size(varargin{1}, 1);
            i_vetted = i_first;
            while i_vetted+1 <= N_records;
                for k = 1:length(varargin)
                    if varargin{k}(i_first) ~= varargin{k}(i_vetted+1)
                        break
                    end
                end
                i_vetted = i_vetted + 1;
            end
            i_last = i_vetted;
        end
        
        
        
        function t = ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME)
            t = double(ACQUISITION_TIME(:, 1)) + double(ACQUISITION_TIME(:, 2)) / 65536;
        end
        
        
        
        % ASSUMES: ACQUISITION_TIME should be uint32.
        function ACQUISITION_TIME = linear_seconds_to_ACQUISITION_TIME(t)            
            t_floor = floor(t);            
            
            ACQUISITION_TIME = ones(size(t), 'uint32');
            ACQUISITION_TIME(:, 1) = uint32(t_floor);
            ACQUISITION_TIME(:, 2) = uint32((t - t_floor) * 65536);
        end
        
        
        
        %==================================================================================
        % ASSUMES: Argument ACQUISITION_TIME refers to the first sample in every snapshot.
        %
        % PR = Per Record
        % sample_frequency : Unit: Hz. Frequency of samples within a snapshot.
        %==================================================================================
        function ACQUISITION_TIME_2 = convert_snapshotsPR_to_samplesPR_ACQUISITION_TIME(  ACQUISITION_TIME_1, N_samples_per_snapshot, sample_frequency  )
            if size(ACQUISITION_TIME_1, 2) ~= 2
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Wrong dimensions on input argumet ACQUISITION_TIME_1.')
            end
            
            N_records = size(ACQUISITION_TIME_1, 1);
            
            % Derive the corresponding column and row vectors.
            t_1           = bicas.data_manager.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_1);  % Column vector
            tr_snapshot_1 = (0:(N_samples_per_snapshot-1)) / sample_frequency;    % Row vector. tr = time relative (does not refer to absolute point in time).
            
            % Derive the corresponding same-sized matrices (one row per snapshot).
            t_1_M           = repmat(t_1,           1,         N_samples_per_snapshot);
            tr_snapshot_1_M = repmat(tr_snapshot_1, N_records, 1                     );
            
            % Add matrices and convert to column vector.
            t_2 = reshape((t_1_M + tr_snapshot_1_M)', N_records*N_samples_per_snapshot, 1);
            
            ACQUISITION_TIME_2 = bicas.data_manager.linear_seconds_to_ACQUISITION_TIME(t_2);
        end


        
        % TEMPORARY FUNCTION - FOR TESTING
        % NOTE: No argument sample_frequency.
        function Epoch_2 = convert_snapshotsPR_to_samplesPR_Epoch_TEMP(  Epoch_1, N_samples_per_snapshot)
            Epoch_2 = repelem(Epoch_1, N_samples_per_snapshot);
        end
        
        
        
        %============================================================================================
        % Convert data from 1 snapshot/record to 1 sample/record (from a matrix to a column vector).
        % 
        % PR = Per Record
        %============================================================================================
        function data_2 = convert_snapshotPR_to_samplePR_DATA(data_1, N_samples_per_snapshot)
            % PROPOSAL: Abolish argument N_samples_per_snapshot? Is already implicit in size(data, 2).
            
            if ndims(data_1) ~= 2
                % NOTE: ndims always returns at least two, which is exactly what we want, also for empty and scalars, and row vectors.
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Wrong dimensions on input argumet ACQUISITION_TIME_1.')
            elseif size(data_1, 2) ~= N_samples_per_snapshot
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Dimensions on input argumet ACQUISITION_TIME_1 does not match N_samples_per_snapshot.')
            end
            
            data_2 = reshape(data_1, size(data_1, 1)*N_samples_per_snapshot, 1);
        end
        
        
        
        %====================================================================================
        % Demultiplex, with only constant factors for calibration.
        % NOTE: For development/testing until there is proper code for using transfer functions.
        %
        % This implements Table 3 and Table 4 in "RPW-SYS-MEB-BIA-SPC-00001-IRF", iss1rev16.
        % Variable names are chosen according to these tables.
        %
        % NOTE: "input"/"output" refers to input/output for the function, but (approximately) the opposite for the BIAS hardware.
        %====================================================================================
        % BIAS_physical_input = ...
        function output = simple_demultiplex(input, mux_set, diff_gain)
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
            
            global ERROR_CODES
            global CONSTANTS
            
            if numel(mux_set) ~= 1 || numel(diff_gain) ~= 1
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal argument value "mux_set" or "diff_gain". Must be scalars (not arrays).')
            end
            
            ALPHA    = CONSTANTS.C.approximate_demuxer.alpha;
            BETA     = CONSTANTS.C.approximate_demuxer.beta;
            switch(diff_gain)
                case 0
                    GAMMA = CONSTANTS.C.approximate_demuxer.gamma_lg;
                case 1
                    GAMMA = CONSTANTS.C.approximate_demuxer.gamma_hg;
                otherwise
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal argument value "diff_gain".')
            end
            
            % Set default values which will remain for
            % variables which are not set by the demuxer.
            NONE_VALUES = ones(size(input.BIAS_1)) * NaN;
            V1_LF     = NONE_VALUES;
            V2_LF     = NONE_VALUES;
            V3_LF     = NONE_VALUES;
            V12_LF    = NONE_VALUES;
            V13_LF    = NONE_VALUES;
            V23_LF    = NONE_VALUES;
            V12_LF_AC = NONE_VALUES;
            V13_LF_AC = NONE_VALUES;
            V23_LF_AC = NONE_VALUES;

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
                    V13_LF    = V12_LF    - V23_LF;
                    V2_LF     = V1_LF     - V12_LF;
                    V3_LF     = V2_LF     - V23_LF;                    
                    V13_LF_AC = V12_LF_AC - V23_LF_AC;
                    
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
                    V13_LF_AC = V12_LF_AC - V23_LF_AC;
                    
                case {5,6,7}
                    errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Not implemented yet for this value of mux_set.')
                    
                otherwise
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal argument value for mux_set.')
            end
            
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
            
        end  % simple_demultiplex
        
        
        
        % NESTED FUNCTION (FUNCTION WITHIN FUNCTION)
        %===========================================================================================        
        % Take CDF data (src) divided into records (points in time) and use that to produce data
        % divided into other records (other points in time).
        %
        % Will produce NaN for values of ACQUISITION_TIME_dest outside the range of
        % ACQUISITION_TIME_src.
        %
        % ASSUMES: data_src is a column vector (i.e. one scalar/record).
        %
        % NOTE: Returned data is double (i.e. not e.g. logical).
        %===========================================================================================        
        function data_dest = nearest_interpolate_records(ACQUISITION_TIME_src, data_src, ACQUISITION_TIME_dest)
            % PROPOSAL: Better name?
            % PROPOSAL: Type cast return variable?
            
            t_src  = bicas.data_manager.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_src);
            t_dest = bicas.data_manager.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_dest);
            
            % Vq = interp1(X,V,Xq,METHOD,EXTRAPVAL) replaces the values outside of the
            % interval spanned by X with EXTRAPVAL.  NaN and 0 are often used for
            % EXTRAPVAL.  The default extrapolation behavior with four input arguments
            % is 'extrap' for 'spline' and 'pchip' and EXTRAPVAL = NaN (NaN +NaNi for 
            % complex values) for the other methods.
            data_dest = interp1(t_src, double(data_src), t_dest);
        end
        
    end

end
