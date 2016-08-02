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
% uniquely referenced by a unique string, a "process data type". These process data variables are all related to each other in a
% web (an acyclic directed graph) describing their dependencies of each other. Initially these process data variables
% are all empty.
%
% There are three forms of process data:
%   (1) elementary input process data  : Supplied BY the user. In practise this should correspond to the content of a CDF.
%   (2) elementary output process data : Returned TO the user. In practise this should correspond to the content of a CDF.
%                                        Derived from a fixed set of other process data.
%   (3) intermediate process data      : Derived from a fixed set of other process data. The outside user never sees these.
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
%   unaware of S/W modes.
%      Ex: Updates to the RPW pipeline, datasets.
%      Ex: Can use a different set of modes when validating using the official RCS validation software at LESIA (not all
%      datasets are available for input there, and one gets error when the S/W descriptor describes modes which require
%      non-existing dataset IDs).
%      Ex: Implement inofficial modes for testing(?).
% -- Easy to keep support for legacy datasets (older versions of CDFs)?
% -- Indirect "caching".
% -- Can reuse processing code?
% -- External code can automatically derive which datasets are required to derive which. ==> Can automatically derive the
% S/W modes in the S/W descriptor.
%
%
%
% DEFINITIONS OF TERMS ---- NEEDS REWORKING
% -- process data = In practise, one (class) internal instance variable representing some form of data in some step of
% the "data processing process".
% -- (process) data type = Type of data that is represented by a string.
% Can be a BICAS input dataset (in memory), a BICAS output dataset (in memory), or some intermediate
% internal type of data. Input and output datasets are represented by their respective dataset IDs.
% -- Elementary input/output (process data) = Data that goes in and out of the data manager as a whole. Not an intermediate data type.
% The word "elementary" is added to distinguish input/output from input/output for a smaller operation.
%
% CONVENTIONS: 
% -- It is implicit that arrays representing CDF data, or "CDF-like" data, use the first MATLAB array index to represent
% CDF records.
%
% IMPLEMENTATION NOTE: This class is intended to:
% - Not be aware of S/W modes as defined in the RCS ICD.
% - Not deal with writing CDFs. It does however read CDF files as a way of setting elementary input process data since
% it is trivial.
%
classdef data_manager
%###################################################################################################
% PROPOSAL: Use other class name that implices processing. "processing_manager"? "proc_manager"?
%
% PROPOSAL: Better name for process data (all together or individual variable), process data type, elementary input/output.
%     PROPOSAL: "process data", "process data variables" = All process data in aggregate?
%     PROPOSAL: "process data variable" (PDV)
%     PROPOSAL: "process data ID" (PDID).
%
% PROPOSAL: Different name for internal data variable. data, datas, data_types, internal_data,
% data_network, data_chain, chain_data, process_data, processing_data. A one-letter name?
%
% QUESTION: Should this code have some kind of facility for handling multiple dataset versions in
% principle, eventhough not formally required? Could be useful.
%    NOTE: Could have multiple internal datatypes, remnants from earlier implementations.
%    PROPOSAL: Add versions to data type strings for data sets.
%
% PROPOSAL: Move out the reading of input CDF files?! Of reading master CDFs?!!
%
% QUESTION: How handle master CDFs?
%   NOTE: Code needs the fill/pad values for the final CDF format. Can not work entirely with NaN internally since some
%         MATLAB classes (non-floating point) type do not have NaN. Also, in principle, there is only one NaN which can not
%         both represent the pad value and the fill value.
%   NOTE: There is exactly one master CDF per elementary output.
%   PROPOSAL: No master CDFs (no data, no pad/fill values) inside data_manager. Convert NaN to pad or fill value outside
%             before writing to file.
%       CON: Some MATLAB classes (non-floating point) type do not have NaN.
%   PROPOSAL: Supply the master CDF with the call to get_process_data_recursively.
%       CON: Won't work if it requires other process data corresponding to another (output) CDF which requires its own
%            master CDF.
%   PROPOSAL: User gives the master CDF to the data_manager which stores it in instance variable (containers.Map for output process data).
%   PROPOSAL: data_manager has access to CONSTANTS/function which retrieves master CDFs.
%
%
% QUESTION: Should this class AT ALL handle reading and writing CDF _FILES_? Should that not be done in some separate layer?!
%    PROPOSAL: Only work with input/output data types that correspond tightly to CDF files.
%       CON: Would work against the idea of combining code with the S/W descriptor info.
%          CON: Not really. How?
%       CON: Might need access to fill/pad values from the CDF templates for the final CDF files.
%          CON: Could make conversions NaN-->pad/fill value outside of data_manager.
%
% PROPOSAL: Implement master CDFs as internal process data types?!!
%    PRO: Does actually need the fill/pad values.
%    CON: Elementary inputs should fit with S/W description. Master CDFs will mess that up.
%
% QUESTION: How obtain the SW root directory and the master CDFs?
%     Some sort of function?
%     Some sort of new constant?
%
% PROPOSAL: Do NOT use dataset IDs (or dataset IDs+version) for data types.
%    PRO: That scheme still assumes that the internal variables format (in memory) is unambiguous.
%         Might have to handle both dataobj and other.
%
%
% QUESTION: Is dataobj a handle class?! How (deep) copy? Needed?! If so, then matters for setting elementary input
% processing data directly.
%
% NOTE: Both BIAS HK and LFR SURV CWF contain MUX data (only LFR has one timestamp per snapshot). True also for other input datasets?
%###################################################################################################


    properties(Access=private)
        % Convention: properties is empty = property has not been set yet.
        
        
        % NOTE: Only set the map's keys here. This makes sure that no disallowed keys are ever used.
        % For data corresponding to CDFs: Use dataset IDs.
        process_data = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        ELEMENTARY_INPUT_DATA_TYPES  = {'ROC-SGSE_L2R_RPW-LFR-SURV-CWF_V01',   'ROC-SGSE_L2R_RPW-LFR-SURV-SWF_V01', 'ROC-SGSE_HK_RPW-BIA_V01', 'ROC-SGSE_L2R_TEST_V99'};
        ELEMENTARY_OUTPUT_DATA_TYPES = {'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01', 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E_V01',                          'ROC-SGSE_L2S_TEST_V99'};
    end
    
    %###############################################################################################

    methods(Access=public)
        
        %=============
        % Constructor
        %=============
        function obj = data_manager()
            ALL_DATA_TYPES = {obj.ELEMENTARY_INPUT_DATA_TYPES{:}, obj.ELEMENTARY_OUTPUT_DATA_TYPES{:}};
            
            validate_strings_unique(ALL_DATA_TYPES)
            
            for data_type = ALL_DATA_TYPES
                obj.process_data(data_type{1}) = [];
            end
        end

        
        
        % NOTE: Not preventing from setting data types/properties which are not CDF files.
        function set_input_CDF(obj, data_type, file_path)
            % PROPOSAL: "Validate" read CDF here.
            irf.log('n', sprintf('data_type=%s: file_path=%s', data_type, file_path))    % NOTE: irf.log adds the method name.
            global ERROR_CODES
            
            if ~ismember(data_type, obj.ELEMENTARY_INPUT_DATA_TYPES)
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Data type is not an elementary input data type.')
            end
            data = dataobj(file_path);
            obj.set_process_data(data_type, data);
        end



        %===========================================================================================
        % Like "get_process_data" but tries to fill in values recursively using other process data.
        %
        % This function derives process data from other process data, e.g. CDF data from other CDF data.
        %===========================================================================================
        function data = get_process_data_recursively(obj, data_type)
            % TODO: General functions for handling CDFs:
            %    -validating input/output CDFs: dataset IDs, (same variable sizes?, same nbr of
            %     records), obtain number of records(?!!).
            %    -setting Parents+Parent_version,
            %    -reading, finding master CDF
            
            % TODO: General functions for handling data.
            %    -Convert fill/pad values <---> NaN.
            
            global ERROR_CODES CONSTANTS
            
            %=======================================
            % Return process data if already exists.
            %=======================================
            data = get_process_data(obj, data_type);
            if ~isempty(data)
                return
            end
            


            %===============================
            % Obtain the input process data
            %===============================
            input_data_types = data_manager.get_input_process_data_types(data_type);
            inputs = [];
            fn_list = fieldnames(input_data_types);
            for i=1:length(fn_list)
                fn = fn_list{i};
                inputs.(fn) = obj.get_process_data_recursively(input_data_types.(fn));
            end
            
            
            
            switch(data_type)
                case 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01'

                    % Might be relevant. True/false depending on available BIAS_1/2/3? Should be in
                    % agreement with MUX_SET?
                    % input_BIAS_HK.BIAS_MODE_BIAS1_ENABLED; ...
                    % input_LFR.data.R0; ...
                    
                    BIAS_HK_cdf = inputs.HK.data;
                    LFR_cdf     = inputs.SCI.data;

                    N_LFR_records          = size(LFR_cdf.POTENTIAL.data, 1);
                    N_samples_per_snapshot = size(LFR_cdf.POTENTIAL.data, 2);
                    
                    % Change to standard names.
                    demuxer_input = [];
                    demuxer_input.BIAS_1 = LFR_cdf.POTENTIAL.data;
                    demuxer_input.BIAS_2 = LFR_cdf.ELECTRICAL.data(:,:,1);
                    demuxer_input.BIAS_3 = LFR_cdf.ELECTRICAL.data(:,:,2);
                    demuxer_input.BIAS_4 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;
                    demuxer_input.BIAS_5 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;

                    DIFF_GAIN = data_manager.nearest_interpolate_records(...
                        BIAS_HK_cdf.ACQUISITION_TIME.data, ...
                        BIAS_HK_cdf.HK_BIA_DIFF_GAIN.data, ...
                        LFR_cdf    .ACQUISITION_TIME.data);
                    
                    %--------------------------------------------------------------------------
                    % IMPLEMENTATION NOTE: One can obtain MUX_SET from (1) LFR or (2) BIAS HK.
                    % This code chooses which one is actually used.
                    %--------------------------------------------------------------------------
                    MUX_SET = data_manager.nearest_interpolate_records(...
                        BIAS_HK_cdf.ACQUISITION_TIME.data, ...
                        BIAS_HK_cdf.HK_BIA_MODE_MUX_SET.data, ...
                        LFR_cdf    .ACQUISITION_TIME.data);
                    %MUX_SET = LFR_cdf.BIAS_MODE_MUX_SET.data;

                    
                    
                    %============================================================
                    % Run demuxer on sequences or records with constant settings
                    %============================================================
                    i_first = 1;
                    demuxer_output = struct(...
                        'V1',  [],    'V2', [], 'V3', [], ...
                        'V12', [],    'V23', [], 'V13', [], ...
                        'V12_AC', [], 'V23_AC', [], 'V13_AC', []);
                    while i_first <= N_LFR_records;
                        
                        % Find sequence of records having identical settings.
                        i_last = data_manager.find_last_same_sequence(...
                            i_first, ...
                            DIFF_GAIN, ...
                            MUX_SET, ...
                            LFR_cdf.FREQ.data);
                        
                        diff_gain        = DIFF_GAIN(i_first);
                        mux_set          = MUX_SET(i_first);
                        sample_frequency = data_manager.get_LFR_samples_in_snapshot_frequency(LFR_cdf.FREQ.data(i_first));
                        
                        demuxer_input_subseq = data_manager.select_subset_from_struct(demuxer_input, i_first, i_last);
                        
                        
                        % QUESTION: How to structure the demuxing?
                        % ----------------------------------------
                        % QUESTION: How split by record? How put together again? How do in a way which
                        % works for transfer functions? How handle the many non-indexed outputs?
                        % QUESTION: How handle changing values of diff_gain, mux_set, bias-dependent calibration offsets?
                        % NOTE: LFR data can be both 1 sample/record or 1 snapshot/record.
                        % PROPOSAL: Work with some subset of in- and out-values of each type?
                        %   PROPOSAL: Work with exactly one value of each type?
                        %       CON: Slow.
                        %           CON: Only temporary implementation.
                        %       PRO: Quick to implement.
                        %   PROPOSAL: Work with only some arbitrary subset specified by array of indices.
                        %   PROPOSAL: Work with only one row?
                        %   PROPOSAL: Work with a continuous sequence of rows/records?
                        %   PROPOSAL: Submit all values, and return structure. Only read and set subset specified by indices.
                        demuxer_output_subseq = data_manager.simple_demultiplex(demuxer_input_subseq, mux_set, diff_gain);
                        
                        demuxer_output = data_manager.add_components_to_struct(demuxer_output, demuxer_output_subseq);
                        
                        i_first = i_last + 1;
                    end   % while



                    %===============================================
                    % Convert 1 snapshot/record --> 1 sample/record
                    %===============================================
                    ACQUISITION_TIME = data_manager.convert_snapshotsPR_to_samplesPR_ACQUISITION_TIME(...
                        LFR_cdf.ACQUISITION_TIME.data, ...
                        N_samples_per_snapshot, ...
                        sample_frequency);
                    Epoch = data_manager.convert_snapshotsPR_to_samplesPR_Epoch_TEMP(  LFR_cdf.Epoch.data, N_samples_per_snapshot  );
                    fn_list = fieldnames(demuxer_output);
                    for i = 1:length(fn_list)
                        fn = fn_list{i};
                        demuxer_output.(fn) = data_manager.convert_snapshotPR_to_samplePR_DATA(  demuxer_output.(fn), N_samples_per_snapshot  );
                    end

                    
                    
                    %============================================================
                    % Put together variables to written to CDF file (elsewhere).
                    %============================================================
                    data = [];
                    data.EPOCH = Epoch;
                    %data.DELTA_PLUS_MINUS = 
                    data.ACQUISITION_TIME = ACQUISITION_TIME;
                    data.V   = [demuxer_output.V1,     demuxer_output.V2,     demuxer_output.V3];
                    data.E   = [demuxer_output.V12,    demuxer_output.V13,    demuxer_output.V23];
                    data.EAC = [demuxer_output.V12_AC, demuxer_output.V13_AC, demuxer_output.V23_AC];                    
                    data.QUALITY_FLAG     = cast(  zeros(size(Epoch)),  convert_CDF_type_to_MATLAB_class('CDF_UINT1', 'Only CDF data types')  );
                    data.QUALITY_BITMASK  = cast(  zeros(size(Epoch)),  convert_CDF_type_to_MATLAB_class('CDF_UINT1', 'Only CDF data types')  );
                    data.DELTA_PLUS_MINUS = cast(  zeros(size(Epoch)),  convert_CDF_type_to_MATLAB_class('CDF_INT8', 'Only CDF data types')  );
                    
                otherwise
                    errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Can not produce this data type, "%s".', data_type)
                    
            end   % switch
        end
        
    end   % methods
    
    %###############################################################################################

    methods(Access=private)        
        
        % Return process data value. Does NOT try to fill it with (derived) data if empty.
        % Can be used to find out whether process data has been set.
        %
        % ASSERTS: Valid data_type.
        % DOES NOT ASSERT: Process data has been set.
        function data = get_process_data(obj, data_type)
            global ERROR_CODES
            
            if ~obj.process_data.isKey(data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is no such process data type, "%s".', data_type);
            end
            data = obj.process_data(data_type);
        end



        % Set a process data variable.
        %
        % ASSERTS: Process data has NOT been set yet.
        % ASSERTS: Valid data_type
        function set_process_data(obj, data_type, data)
            global ERROR_CODES
            
            if ~obj.process_data.isKey(data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is no such process data type, "%s".', data_type);
            elseif ~isempty(obj.process_data(data_type))
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is aldready process data for that data type, "%s".', data_type);
            end
            obj.process_data(data_type) = data;
        end
        
    end   % methods
    
    %###############################################################################################
    
    %methods(Static, Access=private)
    methods(Static, Access=public)

        %===========================================================================================================
        % For every process data type, return a structure with fields set to the data types required to derive that
        % process data type (if any; some data process data types are just set by the caller).
        %
        % This function effectively
        % 1) defines/describes the dependencies between different process data types, and
        % 2) implicitly distinguishes the inputs data types by suitable field names for them (used as field names in other code).
        %
        % Not recursive but designed so that other recursive code can recurse over the implicit "tree data structure" defined by it.
        % Could be used to derive the (elementary) input datasets from the (elementary) output datasets.
        %===========================================================================================================
        function inputs = get_input_process_data_types(data_type)
            % PROPOSAL: Encode data into some data structure (in constants class) instead?!
            global ERROR_CODES
            
            inputs = struct();
            
            switch(data_type)
                case 'ROC-SGSE_HK_RPW-BIA_V01'
                case 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF_V01'
                case 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF_V01'
                case 'ROC-SGSE_L2R_TEST_V99';   % TEST
                    
                case 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E_V01'
                    inputs.HK  = 'ROC-SGSE_HK_RPW-BIA_V01';
                    inputs.SCI = 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF_V01';
                    
                case 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E_V01'
                    inputs.HK  = 'ROC-SGSE_HK_RPW-BIA_V01';
                    inputs.SCI = 'ROC-SGSE_L2R_RPW-LFR-SURV-SWF_V01';
                    
                % TEST
                case 'ROC-SGSE_L2S_TEST_V99'
                    inputs.SCI2  = 'ROC-SGSE_L2R_TEST_V99';
                    
                otherwise
                    errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Can not produce this data type, "%s".', data_type)
            end
        end
        
        
        
        %==================================================================================================
        % Collect all the elementary input process data types needed to produce a given process data type.
        % Can be used to determine which process data types are needed for a S/W mode.
        % Recursive.
        %==================================================================================================
        function process_data_types = get_elementary_input_process_data_types(data_type)
            s = data_manager.get_input_process_data_types(data_type);
            
            process_data_types = {};
            
            if isempty(fieldnames(s))
                % CASE: Found input process data type.
                process_data_types = {data_type};
            else
                % CASE: Found NON-input process data type.
                fn_list = fieldnames(s);
                for i = 1:length(fn_list)
                    fn = fn_list{i};
                
                    new_types = data_manager.get_elementary_input_process_data_types(s.(fn));   % NOTE: RECURSIVE CALL
                    process_data_types = [process_data_types, new_types];
                end
            end

        end
        
        
        
        % Given a struct, select a subset of that struct defined by a range of column indicies for every field.
        function s = select_subset_from_struct(s, i_first, i_last)
            fn_list = fieldnames(s);
            for i=1:length(fn_list)
                fn = fn_list{i};
                
                s.(fn) = s.(fn)(i_first:i_last, :, :);
            end
        end
        
        

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
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Can not handle this LFR_FREQ value.')
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
            t_1           = data_manager.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_1);  % Column vector
            tr_snapshot_1 = (0:(N_samples_per_snapshot-1)) / sample_frequency;    % Row vector. tr = time relative (does not refer to absolute point in time).
            
            % Derive the corresponding same-sized matrices (one row per snapshot).
            t_1_M           = repmat(t_1,           1,         N_samples_per_snapshot);
            tr_snapshot_1_M = repmat(tr_snapshot_1, N_records, 1                     );
            
            % Add matrices and convert to column vector.
            t_2 = reshape((t_1_M + tr_snapshot_1_M)', N_records*N_samples_per_snapshot, 1);
            
            ACQUISITION_TIME_2 = data_manager.linear_seconds_to_ACQUISITION_TIME(t_2);
        end


        
        % TEMPORARY FUNCTION - FOR TESTING
        % NOTE: No argument sample_frequency.
        function Epoch_2 = convert_snapshotsPR_to_samplesPR_Epoch_TEMP(  Epoch_1, N_samples_per_snapshot)
            Epoch_2 = repelem(Epoch_1, N_samples_per_snapshot);
        end
        
        
        
        % Convert data from 1 snapshot/record to 1 sample/record (from a matrix to a column vector).
        % 
        % PR = Per Record
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
        % BIAS = ...
        %====================================================================================
        % BIAS_physical_input = ...
        function output = simple_demultiplex(input, mux_set, diff_gain)
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
            
            global ERROR_CODES
            global CONSTANTS
            
            if numel(mux_set) ~= 1 || numel(diff_gain) ~= 1
                errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal argument value mux_set or diff_gain. Must be scalars (not arrays).')
            end
            
            ALPHA    = CONSTANTS.C.approximate_demuxer.alpha;
            BETA     = CONSTANTS.C.approximate_demuxer.beta;
            switch(diff_gain)
                case 0
                    GAMMA = CONSTANTS.C.approximate_demuxer.gamma_lg;
                case 1
                    GAMMA = CONSTANTS.C.approximate_demuxer.gamma_hg;
                otherwise
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Illegal argument value diff_gain.')
            end
            
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
                    % QUESTION: Overdetermined if we always have BIAS1-3? If so, how select what to
                    % calculate?! What if results disagree?
                    
                    % Summarize what we have;
                    V2_DC  = input.BIAS_1;
                    V3_DC  = input.BIAS_2;
                    V23_DC = input.BIAS_3;
                    %V12_AC = input.BIAS_4;   % Not accessible(?)
                    V23_AC = input.BIAS_5;
                    % Convert that which is trivial.
                    V2_LF     = V2_DC / ALPHA;
                    V3_LF     = V3_DC / ALPHA;
                    V23_LF_AC = V23_DC / GAMMA
                    
                    errorp(ERROR_CODES.OPERATION_NOT_IMPLEMENTED, 'Not certain if correct implementation / incomplete implementation.')
                    
                case {2,3,4,5,6,7}
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
            
        end
        
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
            
            t_src  = data_manager.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_src);
            t_dest = data_manager.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_dest);
            
            % Vq = interp1(X,V,Xq,METHOD,EXTRAPVAL) replaces the values outside of the
            % interval spanned by X with EXTRAPVAL.  NaN and 0 are often used for
            % EXTRAPVAL.  The default extrapolation behavior with four input arguments
            % is 'extrap' for 'spline' and 'pchip' and EXTRAPVAL = NaN (NaN +NaNi for 
            % complex values) for the other methods.
            data_dest = interp1(t_src, double(data_src), t_dest);
        end
        
    end

end
