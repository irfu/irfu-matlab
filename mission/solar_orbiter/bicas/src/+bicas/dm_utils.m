classdef dm_utils
    % Collections of minor utility functions (in the form of static methods) used by data_manager ("dm").
    % The functions are put here to reduce the size of data_manager.
    %
    % Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
    % First created 2016-10-10

    % PROPOSAL: Move some functions to "utils".
    %   Ex: add_components_to_struct, select_subset_from_struct
    
    methods(Static, Access=public)

        function filtered_data = filter_rows(data, row_filter)
            % Function intended for filtering out data from a zVariable.
            % records       : Numeric array with N rows.             (Intended to represent zVariables with N records.)
            % record_filter : Numeric/logical 1D vector with N rows. (Intended to represent zVariables with N records.)
            % filtered_data : Array of the same size as "records", with 
            %                 filtered_data(i,:,:, ...) == NaN,                 for record_filter(i)==0.
            %                 filtered_data(i,:,:, ...) == records(i,:,:, ...), for record_filter(i)~=0.

            % Name? filter_rows? filter_records?
            
            % ASSERTIONS
            if length(row_filter) ~= numel(row_filter)
                % Not necessary to require row vector, only 1D vector.
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'row_filter is not a 1D vector.')  % Use "DatasetFormat"?
            elseif size(row_filter, 1) ~= size(data, 1)
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Numbers of records do not match.')    % Use "DatasetFormat"?
            end
            
            filtered_data = data;
            
            % IMPLEMENTATION NOTE: Command works empirically for filtered_data having any number of dimensions.
            % However, if row_filter and filtered_data have different numbers of rows, then
            % the final array may get the wrong dimensions (without triggering error!) since
            % new array components (indices) are assigned. ==> The corresponding ASSERTION is important!
            filtered_data(row_filter==0, :) = NaN;
        end



        function s = select_subset_from_struct(s, i_first, i_last)
        % Given a struct, select a subset of that struct defined by a range of column indicies for every field.
        % Generic utility function.
        
            fn_list = fieldnames(s);
            for i=1:length(fn_list)
                fn = fn_list{i};
                
                s.(fn) = s.(fn)(i_first:i_last, :, :);
            end
        end
        
        

        function s = add_components_to_struct(s, s_amendment)
        % Add values to every struct field by adding components after their highest column index (let them grow in
        % the column index).
            
            % Generic utility function.
            fn_list = fieldnames(s_amendment);
            for i=1:length(fn_list)
                fn = fn_list{i};
                
                s.(fn) = [s.(fn) ; s_amendment.(fn)];
            end
        end



        function freq = get_LFR_frequency(FREQ)
            % Convert LFR zVariable FREQ constant values to Hz.
            %
            % FREQ : The FREQ zVariable in LFR CDFs.
            % freq : Frequency in Hz.
            
            % PROPOSAL: FREQ assertion in separate function. General assertion-values-from-subset function?
            
            global CONSTANTS
            
            % ASSERTION
            unique_values = unique(FREQ);
            
            if ~all(ismember(unique_values, [0,1,2,3]))
                unique_values_str = sprintf('%d', unique_values);   % NOTE: Has to print without \n to keep all values on a single-line string.
                error('BICAS:dm_utils:Assertion:IllegalArgument:DatasetFormat', 'Found unexpected values in LFR_FREQ (unique values: %s).', unique_values_str)
            end
            
            % NOTE: Implementation that works for arrays of any size.
            freq = ones(size(FREQ)) * -1;
            freq(FREQ==0) = CONSTANTS.C.LFR.F0;
            freq(FREQ==1) = CONSTANTS.C.LFR.F1;
            freq(FREQ==2) = CONSTANTS.C.LFR.F2;
            freq(FREQ==3) = CONSTANTS.C.LFR.F3;
        end
        
        
        
        function Rx = get_LFR_Rx(R0, R1, R2, FREQ)
            % Return the relevant value of LFR CDF zVariables R0, R1, or R2, or a hypothetical "R3" which is always 1.
            %
            % R0, R1, R2, FREQ : LFR CDF zVariables. All must have identical sizes.
            % Rx : The relevant values taken from R0, R1, R2, or "R3" (=1).
            %
            % NOTE: Works for all array sizes.
            
            Rx = -ones(size(FREQ));
            
            I = (FREQ==0); Rx(I) = R0(I);
            I = (FREQ==1); Rx(I) = R1(I);
            I = (FREQ==2); Rx(I) = R2(I);
            I = (FREQ==3); Rx(I) = 1;
        end



        %=====================================================================================================================
        % Finds the greatest i_last such that all varargin{k}(i) are equal for i_first <= i <= i_last separately for every k.
        % Useful for finding a continuous sequence of records with the same data.
        %
        % ASSUMES: varargin{i} are all column arrays of the same size.
        %=====================================================================================================================
        function i_last = find_last_same_sequence(i_first, varargin)
            % PROPOSAL: Better name?
            
            % Check arguments
            for k = 1:length(varargin)
                if ~iscolumn(varargin{k})
                    error('BICAS:dm_utils:Assertion:IllegalArgument', 'varargins are not all column vectors.')
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



%         function t = ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME)
%             % Convert zVariable ACQUISITION_TIME (array size Nx2) to number of seconds since "some" epoch.
%             t = double(ACQUISITION_TIME(:, 1)) + double(ACQUISITION_TIME(:, 2)) / 65536;
%         end
        
        
        
%         function ACQUISITION_TIME = linear_seconds_to_ACQUISITION_TIME(t)            
%             % Convert number of seconds since epoch to zVariable ACQUISITION_TIME.
%             %
%             % ASSUMES: ACQUISITION_TIME should be uint32.
%         
%             t_floor = floor(t);
%             
%             ACQUISITION_TIME = ones(size(t), 'uint32');
%             ACQUISITION_TIME(:, 1) = uint32(t_floor);
%             ACQUISITION_TIME(:, 2) = uint32((t - t_floor) * 65536);
%         end
        
        
        
        function t_tt2000 = ACQUISITION_TIME_to_tt2000(ACQUISITION_TIME)
            % Convert time in from ACQUISITION_TIME to tt2000 which is used for Epoch in CDF files.
            % 
            % NOTE: Return value is in int64.
            
            global CONSTANTS
            
            % ASSERTIONS
            if ndims(ACQUISITION_TIME) ~= 2
                error('BICAS:dm_utils:Assertion:IllegalArgument:DatasetFormat', 'ACQUISITION_TIME has the wrong number of dimensions.')
            elseif size(ACQUISITION_TIME, 2) ~= 2
                error('BICAS:dm_utils:Assertion:IllegalArgument:DatasetFormat', 'ACQUISITION_TIME is does not have two columns.')
            end
            
            t = double(ACQUISITION_TIME(:, 1)) + double(ACQUISITION_TIME(:, 2)) / 65536;
            t_tt2000 = spdfcomputett2000(CONSTANTS.C.ACQUISITION_TIME_EPOCH_UTC) + int64(1e9 * t);   % NOTE: spdfcomputett2000 returns int64 (as it should).
        end
        
        
        
        function UTC_str = tt2000_to_UTC_str(t_tt2000)
            % Convert tt2000 value to UTC string with nanoseconds.
            %
            % Example: 2016-04-16T02:26:14.196334848
            % NOTE: This is the inverse to spdfparsett2000.
            
            v = spdfbreakdowntt2000(t_tt2000);
            UTC_str = sprintf('%04i-%02i-%02iT%02i:%02i:%2i.%03i%03i%03i', v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9));
        end
        
        
        
        %==================================================================================
        % ASSUMES: Argument ACQUISITION_TIME refers to the first sample in every snapshot.
        %
        % PR = Per Record
        % sample_frequency : Unit: Hz. Frequency of samples within a snapshot.
        %==================================================================================
        function ACQUISITION_TIME_2 = convert_snapshotsPR_to_samplesPR_ACQUISITION_TIME(  ACQUISITION_TIME_1, N_samples_per_snapshot, sample_frequency  )
            if size(ACQUISITION_TIME_1, 2) ~= 2
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Wrong dimensions on input argumet ACQUISITION_TIME_1.')
            end
            
            N_records = size(ACQUISITION_TIME_1, 1);
            
            % Derive the corresponding column and row vectors.
            t_1           = bicas.dm_utils.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_1);  % Column vector
            tr_snapshot_1 = (0:(N_samples_per_snapshot-1)) / sample_frequency;    % Row vector. tr = time relative (does not refer to absolute point in time).
            
            % Derive the corresponding same-sized matrices (one row per snapshot).
            t_1_M           = repmat(t_1,           1,         N_samples_per_snapshot);
            tr_snapshot_1_M = repmat(tr_snapshot_1, N_records, 1                     );
            
            % Add matrices and convert to column vector.
            t_2 = reshape((t_1_M + tr_snapshot_1_M)', N_records*N_samples_per_snapshot, 1);
            
            ACQUISITION_TIME_2 = bicas.dm_utils.linear_seconds_to_ACQUISITION_TIME(t_2);
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
            
            % NOTE: ndims always returns at least two, which is exactly what we want, also for empty and scalars, and row vectors.
            if ndims(data_1) ~= 2
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Input argumet ACQUISITION_TIME_1 has the wrong dimensions.')
            elseif size(data_1, 2) ~= N_samples_per_snapshot
                error('BICAS:dm_utils:Assertion:IllegalArgument', 'Dimensions on input argumet ACQUISITION_TIME_1 does not match N_samples_per_snapshot.')
            end
            
            data_2 = reshape(data_1, size(data_1, 1)*N_samples_per_snapshot, 1);
        end
        
        
%         
%         function data_dest = nearest_interpolate_records(ACQUISITION_TIME_src, data_src, ACQUISITION_TIME_dest)
%             % Take CDF data (src) divided into records (points in time) and use that to produce data
%             % divided into other records (other points in time).
%             %
%             % Will produce NaN for values of ACQUISITION_TIME_dest outside the range of
%             % ACQUISITION_TIME_src.
%             %
%             % ASSUMES: data_src is a column vector (i.e. one scalar/record).
%             %
%             % NOTE: Returned data is double (i.e. not e.g. logical).
%         
%             % PROPOSAL: Better name?
%             % PROPOSAL: Type cast return variable?
%             % PROPOSAL: ABOLISH? Should not use functions which are tied to a specific time format (ACQUSITION_TIME vs
%             % Epoch)
%             
%             t_src  = bicas.dm_utils.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_src);
%             t_dest = bicas.dm_utils.ACQUISITION_TIME_to_linear_seconds(ACQUISITION_TIME_dest);
% 
%             % "Vq = interp1(X,V,Xq,METHOD,EXTRAPVAL) replaces the values outside of the
%             % interval spanned by X with EXTRAPVAL.  NaN and 0 are often used for
%             % EXTRAPVAL."
%             % "'linear'   - (default) linear interpolation"
%             data_dest = interp1(t_src, double(data_src), t_dest, 'nearest', NaN);
%         end
%         
    end   % Static
    
end

