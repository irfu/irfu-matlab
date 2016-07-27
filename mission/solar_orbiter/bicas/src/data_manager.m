% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-06-10
%
% "Data manager". Does the actual processing of datasets.
%
% "data type" = Type of data that is represented by a string.
% Can be a BICAS input dataset (in memory), a BICAS output dataset (in memory), or some intermediate
% internal type of data. Input and output datasets are represented by their respective dataset IDs.
%
%
%
% NOTE: It is implicit that arrays representing CDF data, or "CDF-like" use the first MATLAB array index to represent
% CDF records.
%
classdef data_manager
%###################################################################################################
% PROPOSAL: Use other class name that implices processing. "processing_manager"?
%
% QUESTION: Should this class at all handle reading and writing cdf _FILES_?
%    Should that not be done in some separate layer?!
%    PROPOSAL: Only work with input/output data types that correspond tightly to cdf files.
%       CON: Would work against the idea of combining code with the S/W descriptor info.
%
% QUESTION: Should this code have some kind of facility for handling multiple dataset versions in
% principle, eventhough not formally required? Could be useful.
%    NOTE: Could have multiple internal datatypes, remnants from earlier implementations.
%    PROPOSAL: Add versions to data type strings for data sets.
%
% PROPOSAL: Different name for internal data variable. data, datas, data_types, internal_data,
% data_network, data_chain, chain_data, process_data, processing_data. A one-letter name?
%
% PROPOSAL: Implement master cdfs as internal process data types?!!
%    NOTE: Should fit with S/W description.
%
% NOTE: Both BIAS HK and LFR SURV CWF contain MUX data (only LFR has one timestamp per snapshot). True also for other input datasets?
%
% PROPOSAL: Move out the reading of input cdf files?! Of reading master cdfs?!!
%
% QUESTION: How obtain the SW root directory and the master cdfs?
%     Some sort of function?
%     Some sort of new constant?
%
% PROPOSAL: Do NOT use dataset IDs (or dataset IDs+version) for data types.
%    PRO: That scheme still assumes that the internal variables format (in memory) is unambiguous.
%         Might have to handle both dataobj and other.
%
%
%
% QUESTION: Is dataobj a handle class?! How (deep) copy? Needed?!
%###################################################################################################


    properties(Access=private)
        % Convention: properties is empty = property has not been set yet.
        
        
        % NOTE: Only set keys here. This makes sure that no disallowed keys are ever used.
        % For data corresponding to cdfs: Use dataset IDs.
        process_data = containers.Map('KeyType', 'char', 'ValueType', 'any');
        
        test_property = [123];
    end
    
    %###############################################################################################

    methods (Access=public)
        
        % Constructor
        function obj = data_manager()
            INPUT_CDFS  = {'ROC-SGSE_L2R_RPW-LFR-SURV-CWF',   'ROC-SGSE_L2R_RPW-LFR-SURV-SWF', 'ROC-SGSE_HK_RPW-BIA'};
            OUTPUT_CDFS = {'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E', 'ROC-SGSE_L2S_RPW-LFR-SURV-SWF-E'};
            ALL_DATA_TYPES = {INPUT_CDFS{:}, OUTPUT_CDFS{:}};
            validate_strings_unique(ALL_DATA_TYPES)
            
            for data_type = ALL_DATA_TYPES
                obj.process_data(data_type{1}) = [];
            end
        end

        
        
        % NOTE: Not preventing from setting data types/properties which are not cdf files.
        function set_input_cdf(obj, data_type, file_path)
            % PROPOSAL: "Validate" cdf here.
            global ERROR_CODES
            irf.log('n', sprintf('data_type=%s: file_path=%s', data_type, file_path))    % NOTE: irf.log adds the method name.
            
            data = dataobj(file_path);
            obj.set_process_data(data_type, data);
        end
        
        
        
        % Like get_process_data but tries to fill in values recursively using other process data.
        % PROPOSAL: Better name. get_data_recursively?
        function data = get_data(obj, data_type)
            % TODO: General functions for handling cdfs:
            %    -validating input/output cdfs: dataset IDs, (same variable sizes?, same nbr of
            %    records), obtain number of records(?!!).
            %    -setting Parents+Parent_version,
            %    -reading master cdf
            
            % TODO: General functions for handling data.
            %    -Convert fill/pad values <---> NaN.
            
            global ERROR_CODES CONSTANTS
            
            data = get_process_data(obj, data_type);   % Implicitly checks if the data type exists.
            if ~isempty(data)
                return
            end
            
            C = CONSTANTS.get_general();
            %C_output = CONSTANTS.get_cdf_outputs_constants({data_type});
            %C_output = output{1};
            
            
            
            switch(data_type)
                case 'ROC-SGSE_L2S_RPW-LFR-SURV-CWF-E'
                    % PROPOSAL: Change names RCS_input/output? Implies input/output files.

                    input_BIAS_HK = get_data(obj, 'ROC-SGSE_HK_RPW-BIA');
                    input_LFR     = get_data(obj, 'ROC-SGSE_L2R_RPW-LFR-SURV-CWF');
                    BIAS_HK_cdf = input_BIAS_HK.data;
                    LFR_cdf     = input_LFR.data;

                    % Might be relevant. True/false depending on available BIAS_1/2/3? Should be in
                    % agreement with MUX_SET?
                    % input_BIAS_HK.BIAS_MODE_BIAS1_ENABLED;
                    % input_BIAS_HK.BIAS_MODE_BIAS2_ENABLED;
                    % input_BIAS_HK.BIAS_MODE_BIAS3_ENABLED;
                    %LFR_R0                  = input_LFR.data.R0;
                    %LFR_R1                  = input_LFR.data.R1;
                    %LFR_R2                  = input_LFR.data.R2;
                    
                    N_LFR_records          = size(LFR_cdf.POTENTIAL.data, 1);
                    N_samples_per_snapshot = size(LFR_cdf.POTENTIAL.data, 2);
                    
                    % Change to standard names.
                    RCS_input = [];
                    RCS_input.BIAS_1 = LFR_cdf.POTENTIAL.data;
                    RCS_input.BIAS_2 = LFR_cdf.ELECTRICAL.data(:,:,1);
                    RCS_input.BIAS_3 = LFR_cdf.ELECTRICAL.data(:,:,2);
                    RCS_input.BIAS_4 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;
                    RCS_input.BIAS_5 = ones(size(LFR_cdf.POTENTIAL.data)) * NaN;

                    DIFF_GAIN = data_manager.nearest_interpolate_records(...
                        BIAS_HK_cdf.ACQUISITION_TIME.data, ...
                        BIAS_HK_cdf.HK_BIA_DIFF_GAIN.data, ...
                        LFR_cdf    .ACQUISITION_TIME.data);
                    
                    %------------------------------------------------------------------------
                    % IMPLEMENTATION NOTE: One can obtain MUX_SET from both LFR and BIAS HK.
                    % This code chooses which one is used.
                    %------------------------------------------------------------------------
                    MUX_SET = data_manager.nearest_interpolate_records(...
                        BIAS_HK_cdf.ACQUISITION_TIME.data, ...
                        BIAS_HK_cdf.HK_BIA_MODE_MUX_SET.data, ...
                        LFR_cdf    .ACQUISITION_TIME.data);
                    %MUX_SET = LFR_cdf.BIAS_MODE_MUX_SET.data;

                    
                    
                    %============================================================
                    % Run demuxer on sequences or records with constant settings
                    %============================================================
                    i_first = 1;
                    RCS_output = struct(...
                        'V1',  [],    'V2', [], 'V3', [], ...
                        'V12', [],    'V23', [], 'V13', [], ...
                        'V12_AC', [], 'V23_AC', [], 'V13_AC', []);
                    while i_first <= N_LFR_records;
                        
                        i_last = data_manager.find_last_same_sequence(...
                            i_first, ...
                            DIFF_GAIN, ...
                            MUX_SET, ...
                            LFR_cdf.FREQ.data);
                        
                        diff_gain        = DIFF_GAIN(i_first);
                        mux_set          = MUX_SET(i_first);
                        sample_frequency = data_manager.get_LFR_samples_in_snapshot_frequency(LFR_cdf.FREQ.data(i_first));
                        
                        RCS_input_subsequence = data_manager.select_subset_from_struct(RCS_input, i_first, i_last);
                        
                        
                        % How to structure the demuxing?
                        % ------------------------------
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
                        RCS_output_subsequence = data_manager.simple_demultiplex(RCS_input_subsequence, mux_set, diff_gain);
                        
                        RCS_output = data_manager.add_components_to_struct(RCS_output, RCS_output_subsequence);
                        
                        i_first = i_last + 1;
                    end   % while
                    
                    
                    %===============================================
                    % Convert 1 snapshot/record --> 1 sample/record
                    %===============================================
                    ACQUISITION_TIME = data_manager.convert_snapshotsPR_to_samplesPR_ACQUISITION_TIME(...
                        LFR_cdf.ACQUISITION_TIME.data, ...
                        N_samples_per_snapshot, ...
                        sample_frequency);                        
                    fn_list = fieldnames(RCS_output);
                    for i = 1:length(fn_list)
                        fn = fn_list{i};
                        RCS_output.(fn) = data_manager.convert_snapshotsPR_to_samplesPR_DATA(RCS_output.(fn), N_samples_per_snapshot);
                    end

                    % NOTE: Needs ~C.sw_modes{1}.outputs{1}.master_cdf_filename
                    master_cdf_path = [CONSTANTS.SW_root_dir(), filesep, C.master_cdfs_dir_rel, filesep, output.master_cdf_filename];
                    master_cdf = dataobj(master_cdf_path);
                    
                    %data = [];
                    
                    %master_cdf.data.ACQUISITION_TIME.data = LFR_ACQUISITION_TIME;  % Incorrect!
                    % TODO: Update nrec? Needed when writing to file?!
                    
                otherwise
                    errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'Can not produce this data type, "%s".', data_type)
                    
            end   % switch
        end
        
    end   % methods
    
    %###############################################################################################
    
    %methods(Static, Access=private)
    methods(Static, Access=public)

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
            C = CONSTANTS.get_general();
            
            switch(LFR_FREQ)
                case 0
                    freq = C.LFR.F0;
                case 1
                    freq = C.LFR.F1;
                case 2
                    freq = C.LFR.F2;
                case 3
                    freq = C.LFR.F3;
                otherwise
                    errorp(ERROR_CODES.ASSERTION_ERROR, 'Can not handle this LFR_FREQ value.')
            end
        end



        % Finds the greatest i_last such that all varargin{k}(i) are equal for i_first <= i <= i_last separately for every k.
        % Useful for finding a continuous sequence of records with the same data.
        %
        % ASSUMES: varargin{i} are all column arrays of the same size.
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
        
        
        
        % ASSUMES: Argument ACQUISITION_TIME refers to the first sample in every snapshot.
        % PR = Per Record
        % sample_frequency : Unit: Hz. Frequency of samples within a snapshot.
        function ACQUISITION_TIME_2 = convert_snapshotsPR_to_samplesPR_ACQUISITION_TIME(ACQUISITION_TIME_1, N_samples_per_snapshot, sample_frequency)
            N_records = size(ACQUISITION_TIME_1, 1);
            
            % Derive the corresponding column and row vectors.
            t_1           = data_manager.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_1);
            tr_snapshot_1 = (0:(N_samples_per_snapshot-1)) / sample_frequency;    % tr = time relative (does not refer to absolute point in time).
            
            % Derive the corresponding same-sized matrices (one row per snapshot).
            t_1_M           = repmat(t_1,           1,         N_samples_per_snapshot);
            tr_snapshot_1_M = repmat(tr_snapshot_1, N_records, 1                     );
            
            % Add matrices and convert to column vector.
            t_2 = reshape((t_1_M + tr_snapshot_1_M)', N_records*N_samples_per_snapshot, 1);
            
            ACQUISITION_TIME_2 = data_manager.linear_seconds_to_ACQUISITION_TIME(t_2);
        end
        
        
        
        % PR = Per Record
        function data_2 = convert_snapshotsPR_to_samplesPR_DATA(data_1, N_samples_per_snapshot)
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
            
            ALPHA    = CONSTANTS.approximate_demuxer.alpha;
            BETA     = CONSTANTS.approximate_demuxer.beta;
            switch(diff_gain)
                case 0
                    GAMMA = CONSTANTS.approximate_demuxer.gamma_lg;
                case 1
                    GAMMA = CONSTANTS.approximate_demuxer.gamma_hg;
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
        % Take cdf data (src) divided into records (points in time) and use that to produce data
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

    %###############################################################################################

    methods(Access=private)        
        
        % Return what is in the internal data structure, without trying to fill it with data.
        function data = get_process_data(obj, data_type)
            global ERROR_CODES
            
            if ~obj.process_data.isKey(data_type)
                errorp(ERROR_CODES.SW_MODE_PROCESSING_ERROR, 'There is no such process data type, "%s".', data_type);
            end
            data = obj.process_data(data_type)
        end



        % Set the value associated with a key in the process_data.
        % Will make sure that an unused key is not set.
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
    
end
