%
% Collection of minor utility functions (in the form of static methods) used for
% data processing.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-10-10
%
classdef utils
%
% PROPOSAL: POLICY: No functions which set "policy"/configure the output of
%           datasets.
%
% PROPOSAL: Split up in separate files?!
%   PROPOSAL: Move all functions that are used outside bicas.proc to
%             "bicas.utils" (class). -- DONE
%   PROPOSAL: Move functions that are entirely used within a group/package of processing modes
%           (L1L2, L2L2, L2L3) to their respective group.
%       Ex: L1L2
%           Ex: lfr.m, dts.m: derive_DELTA_PLUS_MINUS
%           Ex: lfr.m:        bicas.proc.utils.set_NaN_rows
%           Ex: tds.m:        set_NaN_end_of_rows
%           Ex: dc.m:
%                   select_row_range_from_cell_comps
%                   assert_cell_array_comps_have_same_N_rows (used by select_row_range_from_cell_comps)
%                   convert_matrix_to_cell_array_of_vectors
%           Ex: L1L2.m:       ACQUISITION_TIME_to_TT2000
%                             (not used in processing function directly, only indirectly)
%       CON: Useful to have all SMALL functions collected.
%       CON: Can forget that a function exists if it is needed somewhere else.
%       TODO-DEC: Move to subfile if only used there instead?
%           Ex: Move function to lfr.m instead of L1L2.m?
%
% PROPOSAL: More test code.
% PROPOSAL: Write test code for ACQUISITION_TIME_to_TT2000 and its inversion.
%
% N-->1 sample/record
%    NOTE: Time conversion may require moving the zero-point within the snapshot/record.
%    PROPOSAL: : All have nSamplesPerOldRecord as column vector.
%       PRO: LFR.
%    PROPOSAL: First convert column data to 2D data (with separate functions),
%              then reshape to 1D with one common function.
%       CON: Does not work for ACQUISITION_TIME since two columns.
%
% PROPOSAL: Replace functions
%               set_struct_field_rows()                     -- DONE
%               select_row_range_from_cell_comps()          -- ABOLISHED
%               assert_struct_num_fields_have_same_N_rows()
%       with new class that has a map from arbitrary value to arrays and/or
%       instances of same class (recursive).
%   PRO: Can enforce same number of rows.
%   PRO: Can simultaneously (1) iterate over fields, and (2) identify fields
%        using non-number, e.g. strings.
%   PROPOSAL: Permit cell arrays. Replace functions
%           select_row_range_from_cell_comps()
%           convert_matrix_to_cell_array_of_vectors()
%           convert_cell_array_of_vectors_to_matrix()
%           assert_cell_array_comps_have_same_N_rows() -- ABOLISHED
%       with methods.



    methods(Static, Access=public)



        % Wrapper around bicas.handle_struct_name_change() to be used
        % locally.
        %
        % ARGUMENTS
        % =========
        % inSciDsi : Input SCI DSI which contains the zVariable.
        % varargin : Passed on to bicas.handle_struct_name_change as its
        %            varargin.
        %
        function handle_ZV_name_change(fnChangeList, inSciDsi, Bso, L, varargin)
            anomalyDescrMsgFunc = @(oldFieldname, newFieldname) (sprintf(...
                ['Input dataset DSI=%s uses an alternative', ...
                ' but illegal(?) zVariable name "%s" instead of "%s".'], ...
                inSciDsi, oldFieldname, newFieldname));

            bicas.handle_struct_name_change(...
                fnChangeList, ...
                Bso, L, anomalyDescrMsgFunc, varargin{:})
        end
        
        
        
        function assert_increasing(v, isMonotonic, errorId, msg)
            assert(isnumeric(v) && isvector(v))
            assert(islogical(isMonotonic) && isscalar(isMonotonic))

            if isMonotonic
                option = 'strictascend';
            else
                option = 'ascend';
            end

            if ~issorted(v, option)
                % Indices at which Epoch does not increment
                iDecrArray = find(diff(v) < 0) + 1;
                decrStr = sprintf(...
                    '\nIndices with negative increment: [%s]', ...
                    strjoin(string(iDecrArray), ', ') ...
                );

                if ~isMonotonic
                    error(errorId, [msg, decrStr])
                else
                    iEqualArray = find(diff(v) == 0) + 1;
                    equalStr = sprintf(...
                        '\nIndices with zero increment: [%s]', ...
                        strjoin(string(iEqualArray), ', ') ...
                    );
                    error(errorId, [msg, equalStr, decrStr])
                end
            end

        end



        %################################
        % MODIFYING, DERIVING ZVARIABLES
        %################################



        % Function intended for filtering out data from a zVariable by setting
        % parts of it to NaN. Also useful for constructing aonymous functions.
        %
        %
        % ARGUMENTS
        % =========
        % zv
        %       Numeric array with N rows.
        % bRowFilter
        %       Numeric/logical column vector with N rows.
        %
        %
        % RETURN VALUE
        % ============
        % zv
        %       Array of the same size as argument "zv", such that
        %       zv(i,:,:) == NaN for bRowFilter(i)==true.
        function zv = set_NaN_rows(zv, bRowFilter)
        % PROPOSAL: Implement using set_NaN_end_of_rows().
        %   Ex: zv = set_NaN_after_snapshots_end(zv, size(zv, 2) .*  bRowFilter)
        %   CON: That function can not handle an arbitrary number of dimensions.

            % ASSERTIONS
            assert(isfloat(zv), ...
                'BICAS:Assertion:IllegalArgument', ...
                ['Argument "zv" is not a floating-point class (can', ...
                ' therefore not represent NaN).'])
            assert(islogical(bRowFilter))
            % Not really necessary to require row vector, only 1D vector.
            irf.assert.sizes(...
                zv,         [-1, NaN, NaN], ...
                bRowFilter, [-1])

            % Overwrite selected rows with NaN
            % --------------------------------
            % IMPLEMENTATION NOTE: Command works empirically for zv having any
            % number of dimensions. However, if rowFilter and zv have different
            % numbers of rows, then the final array may get the wrong dimensions
            % (without triggering error!) since new array components (indices)
            % are assigned. ==> Having a corresponding ASSERTION is important!
            zv(bRowFilter, :) = NaN;
        end



        % ARGUMENTS
        % =========
        % zv
        %       2-D floating-point array.
        % snapshotLengths
        %       1-D numeric array.
        %
        % RETURN VALUE
        % ============
        % zv
        %       2-D floating-point array copied from input argument, except that
        %       zv(iRecord, (snapshotLengths(iRecord)+1):end) = NaN;
        function zv = set_NaN_end_of_rows(zv, snapshotLengths)
            % ASSERTIONS
            assert(isfloat(zv))
            [nRecords, snapshotMaxLength] = irf.assert.sizes(...
                zv,              [-1, -2], ...
                snapshotLengths, [-1]);
            % IMPLEMENTATION NOTE: Add zero to vector (sent to max()) so that
            % max() gives a sensible value for empty snapshotLengths.
            assert(snapshotMaxLength >= max([snapshotLengths; 0]))

            % IMPLEMENTATION
            for iRecord = 1:nRecords
                zv(iRecord, (snapshotLengths(iRecord)+1):end) = NaN;
            end
        end
        
        
        
        % Interpret LRX flag as available BLTS's.
        %
        % IMPLEMENTATION NOTE: Function partly exists so that it can be extended
        % to do the authoritative interpretation of LRX for all of BICAS, e.g.
        % for converting input ZVs (in particular LFR) into BLTS-sorted ZVs.
        %
        % ARGUMENTS
        % =========
        % lrx
        %       Scalar float. May be NaN=Unknown.
        %
        function [iBltsAr] = interpret_LRX(lrx)
            % NEED: Handle LRX=unknown

            % PRO: One interpretation of values for all of BICAS.
            % PROBLEM: How fit with possible future implementation that stores
            %          all data in 3 channels with nominally always data
            %          (instead of 5 BLTS's, with some missing data)?
            %     PROPOSAL: Always return 3x1 array, with one value per such
            %               channel.
            % PROPOSAL: LRX (Nx1) --> iBltsArray (Nx3)
            %     PRO: Can possibly be used for copying CDF ZV elements
            %          into internal variables (split by BLTS).
            %     NOTE: Unnecessary if BICAS internally stores samples in
            %           3 data channels instead, since that is what both LFR and
            %           TDS L1R datasets effecitvely do.
            % PROPOSAL: Return variables describing conversion in "both
            %           directions" between 3 data channels and 5 BLTS's.
            %           LRX (1xN) --> 
            %               iBltsArray=[iBltsA, iBltsB, iBltsC]
            %                   iBltsX = 1,2,3,4,5,NaN/FP
            %               hasBltsData=[blts1HasData, ..., blts5HasData]
            %                   bltsXhasData = true,false,NaN/FP

            assert(isscalar(lrx) && isfloat(lrx) && ismember(lrx, [0, 1, NaN]))

            if lrx == 0
                iBltsAr = [1,4,5];
            elseif lrx == 1
                iBltsAr = [1,2,3];
            else
                iBltsAr = [1];
            end
        end



        function tt2000 = ACQUISITION_TIME_to_TT2000(...
                ACQUISITION_TIME, ACQUISITION_TIME_EPOCH_UTC)
        %
        % Convert time in from ACQUISITION_TIME to tt2000 which is used for
        % Epoch in CDF files.
        %
        %
        % ARGUMENTS
        % =========
        % ACQUSITION_TIME            : NOTE: Can be negative since it is uint32.
        % ACQUISITION_TIME_EPOCH_UTC :
        %               Numeric row vector. The time in UTC at
        %               which ACQUISITION_TIME == [0,0], expressed as
        %               [year, month, day, hour, minute, second,
        %                millisecond, microsecond(0-999), nanoseconds(0-999)].
        %
        %
        % RETURN VALUE
        % ============
        % tt2000 : NOTE: int64
        %

            bicas.utils.assert_ZV_ACQUISITION_TIME(ACQUISITION_TIME)

            % at = ACQUISITION_TIME
            ACQUISITION_TIME = double(ACQUISITION_TIME);
            atSeconds = ACQUISITION_TIME(:, 1) + ACQUISITION_TIME(:, 2) / 65536;
            % NOTE: spdfcomputett2000 returns int64 (as it should).
            tt2000 = spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC) ...
                + int64(atSeconds * 1e9);
        end



        function ACQUISITION_TIME = TT2000_to_ACQUISITION_TIME(...
                tt2000, ACQUISITION_TIME_EPOCH_UTC)
        %
        % Convert from tt2000 to ACQUISITION_TIME.
        %
        % ARGUMENTS
        % =========
        % t_tt2000
        %       Nx1 vector. Required to be int64 like the real zVar Epoch.
        %
        % RETURN VALUE
        % ============
        % ACQUISITION_TIME : Nx2 vector. uint32.
        %       NOTE: ACQUSITION_TIME can not be negative since it is uint32.
        %

            % ASSERTIONS
            bicas.utils.assert_ZV_Epoch(tt2000)

            % NOTE: Important to type cast to double because of multiplication
            % AT = ACQUISITION_TIME
            atSeconds = double(int64(tt2000) - ...
                spdfcomputett2000(ACQUISITION_TIME_EPOCH_UTC)) * 1e-9;

            % ASSERTION: ACQUISITION_TIME must not be negative.
            if any(atSeconds < 0)
                error(...
                    'BICAS:Assertion:IllegalArgument:DatasetFormat', ...
                    ['Can not produce ACQUISITION_TIME (uint32) with', ...
                    ' negative number of integer seconds.'])
            end

            atSeconds = round(atSeconds*65536) / 65536;
            atSecondsFloor = floor(atSeconds);

            ACQUISITION_TIME = uint32([]);
            ACQUISITION_TIME(:, 1) = uint32(atSecondsFloor);
            ACQUISITION_TIME(:, 2) = uint32((atSeconds - atSecondsFloor) * 65536);
            % NOTE: Should not be able to produce ACQUISITION_TIME(:, 2)==65536
            % (2^16) since atSeconds already rounded (to parts of 2^-16).
        end



        function zv_DELTA_PLUS_MINUS = derive_DELTA_PLUS_MINUS(freqHz, nSpr)
        %
        % Derive value for zVar DELTA_PLUS_MINUS.
        %
        % NOTE: All values on any given row (CDF record) of DELTA_PLUS_MINUS are
        % identical. Not sure why multiple values per row are needed but it is
        % probably intentional, as per YK's instruction.
        %
        %
        % ARGUMENTS
        % =========
        % freqHz
        %       Frequency column vector in s^-1.
        %       Can not handle freqHz=NaN since the output is an integer
        %       (assertion).
        % nSpr
        %       Number of samples/record.
        %
        %
        % RETURN VALUE
        % ============
        % DELTA_PLUS_MINUS
        %       Analogous to BIAS zVariable. CDF_INT8=int64.
        %       NOTE: Unit ns.

            ZV_DELTA_PLUS_MINUS_DATA_TYPE = 'CDF_INT8';

            % ASSERTIONS
            nRecords = irf.assert.sizes(freqHz, [-1]);
            assert(isfloat(freqHz) && all(isfinite(freqHz)), ...
                'BICAS:Assertion:IllegalArgument', ...
                'Argument "freqHz" does not consist of non-NaN floats.')
            assert(isscalar(nSpr), ...
                'BICAS:Assertion:IllegalArgument', ...
                'Argument "nSpr" is not a scalar.')



            zv_DELTA_PLUS_MINUS = zeros([nRecords, nSpr]);
            %DELTA_PLUS_MINUS = zeros([nRecords, 1]);    % Always 1 sample/record.
            for i = 1:nRecords
                % NOTE: Converts [s] (1/freqHz) --> [ns] (DELTA_PLUS_MINUS) so
                % that the unit is the same as for Epoch.
                % NOTE: Seems to work for more than 2D.
                % Unit: nanoseconds
                zv_DELTA_PLUS_MINUS(i, :) = 1./freqHz(i) * 1e9 * 0.5;
            end
            zv_DELTA_PLUS_MINUS = cast(zv_DELTA_PLUS_MINUS, ...
                irf.cdf.convert_CDF_type_to_MATLAB_class(...
                    ZV_DELTA_PLUS_MINUS_DATA_TYPE, 'Only CDF data types'));
        end



        %######################################################################
        % Convert between
        % (1) 2D non-cell array, and
        % (2) 1D column cell array of 1D arrays, one per (2D array) source row
        %######################################################################



        % Convert 2D array --> 1D cell array of 1D arrays, one per source row.
        %
        %
        % ARGUMENTS
        % =========
        % M
        %       2D matrix
        % nCopyColsPerRowArray
        %       1D column vector. Numeric.
        %       (i) = Number of elements to copy from M(i,:).
        %
        % RETURN VALUE
        % ============
        % ca
        %       Column cell array of 1D vectors.
        %       ca{i}(j). j = 1:nCopyColsPerRowArray(i)
        %
        function ca = convert_matrix_to_cell_array_of_vectors(M, nCopyColsPerRowArray)

            % ASSERTIONS
            assert(isnumeric(M))
            irf.assert.vector(nCopyColsPerRowArray)
            nRows = irf.assert.sizes(...
                M,                    [-1, NaN], ...
                nCopyColsPerRowArray, [-1, 1]);

            % Create "ca".
            ca = cell(nRows, 1);
            for iRow = 1:nRows
                ca{iRow} = M(iRow, 1:nCopyColsPerRowArray(iRow));
            end
        end



        % ARGUMENTS
        % =========
        % ca
        %       Column cell array of 1D vectors.
        % nMatrixColumns
        %       Scalar. Number of columns in M.
        %
        %
        % RETURN VALUES
        % =============
        % M
        %       Numeric 2D matrix.
        %       NOTE: Sets unset elements to NaN.
        % nCopyColsPerRowVec
        %       1D vector.
        %       {i}=Length of ca{i}=Number of elements copied to M{i,:}.
        function [M, nCopyColsPerRowVec] = ...
                convert_cell_array_of_vectors_to_matrix(ca, nMatrixColumns)
            
            assert(iscell(ca))
            assert(iscolumn(ca))
            assert(isscalar(nMatrixColumns))

            nRows              = numel(ca);
            nCopyColsPerRowVec = zeros(nRows, 1);   % Always column vector.
            M                  = nan(  nRows, nMatrixColumns);
            for iRow = 1:numel(nCopyColsPerRowVec)
                nCopyColsPerRowVec(iRow)            = numel(ca{iRow});
                M(iRow, 1:nCopyColsPerRowVec(iRow)) = ca{iRow};
            end

        end



        %######################################
        % Cell array of arrays (cell/non-cell)
        %######################################



        % Assert that all cell array components have the same number of rows.
        % This is useful when cell array components represent zVar-like data,
        % where rows represent CDF records.
        %
        % ARGUMENTS
        % =========
        % ca
        %       N-dim cell array.        
        %       size(ca{i}, 1) for all "i" must evaluate to the same value for
        %       assertion to pass.
        %
%         function assert_cell_array_comps_have_same_N_rows(ca)
%             nRowsArray = cellfun(@(v) (size(v, 1)), ca, 'UniformOutput', true);
%             irf.assert.all_equal( nRowsArray )
%         end



        % For every cell in a cell array, select a (non-cell array) index range
        % in the first dimension for every cell array component.
        %
        % RETURN VALUE
        % ============
        % ca2
        %       ca2{i} = ca1{i}(iFirst:iLast, :, :,:,:,:);
        %
%         function ca2 = select_row_range_from_cell_comps(ca1, iFirst, iLast)
% 
%             % ASSERTIONS
%             bicas.proc.utils.assert_cell_array_comps_have_same_N_rows(ca1)
% 
%             for i = 1:numel(ca1)
%                 ca2{i} = ca1{i}(iFirst:iLast, :, :,:,:,:);
%             end
%         end



        %####################
        % ~Struct of arrays
        %####################



        function nRows = assert_struct_num_fields_have_same_N_rows(S)
        % Assert that struct consisting of
        %   logical arrays
        %   numeric arrays
        %   bicas.utils.SameRowsMap
        % has the same number rows in all its arrays.
        %
        % Useful for structs where all fields represent CDF zVariables and/or
        % derivatives thereof, the size in the first index (number of CDF
        % record) should be equal.
        %
        %
        % ARGUMENTS
        % =========
        % S : Struct
        %       Fields may (ony) be of the following types. Number of rows must
        %       be identical for the specified data structure components
        %       (right-hand side).
        %       (1) numeric/logical fields  : The field itself.
        %       (2) cell fields             : Cell array components (not
        %                                     the cell array itself!).
        %       (3) struct field (in struct): The inner struct's fields (not
        %                                     recursive).
        %
        %
        % ACTUAL USAGE OF SPECIAL CASES FOR FIELDS (non-array fields)
        % ===========================================================
        % PostDc.Zv.AsrSamplesAVolt : bicas.utils.SameRowsMap.

        % NOTE: Function name somewhat bad.
        % PROPOSAL: Make recursive?!
        % PROPOSAL: Implement using new features in irf.assert.sizes().
        % TODO-NI: Function only used for cases where ALL fields should have
        %          same number of rows? (Due to previous refactoring.)
        % PROBLEM: Function is an obstacle to converting AsrSamplesAVoltSrm to
        %          a class.
        %   NOTE: An AsrSamplesAVoltSrm class should ensure INTERNALLY consistent array
        %         sizes, but not consistent with a parent struct.
        %   PROPOSAL: Special case for the specific class.
        %   PROPOSAL: Convert all "Zv" data structs into class that enforces what
        %             this function tests.
        %   PROPOSAL: Redefine as  ~assert_ZV_struct().
        %       PRO: Special case for future class AsrSamplesAVoltSrm is more
        %            natural.

            fieldNamesList1 = fieldnames(S);
            nRowsArray = [];
            for iFn1 = 1:length(fieldNamesList1)
                fieldValue = S.(fieldNamesList1{iFn1});

                if isnumeric(fieldValue) || islogical(fieldValue) || isa(fieldValue, 'bicas.utils.FPArray')
                    % CASE: Numeric & logical field.

                    nRowsArray(end+1) = size(fieldValue, 1);

%                 elseif iscell(fieldValue)
%                     % CASE: Cell array
% 
%                     for iCc = 1:numel(fieldValue)
%                         nRowsArray(end+1) = size(fieldValue{iCc}, 1);
%                     end

%                 elseif isstruct(fieldValue)
%                     % CASE: Struct
%                     % Check number of rows in every field (regardless of type).
% 
%                     fieldNamesList2 = fieldnames(fieldValue);
%                     for iFn2 = 1:length(fieldNamesList2)
%                         nRowsArray(end+1) = size(...
%                             fieldValue.(fieldNamesList2{iFn2}), ...
%                             1);
%                     end

                elseif isa(fieldValue, 'bicas.utils.SameRowsMap')
                    nRowsArray(end+1) = fieldValue.nRows();
                    
                else
                    % CASE: Other field value type.
                    error('BICAS:Assertion', ...
                        'Can not handle this type of struct field.')
                end
            end

            % NOTE: Empty vector if nRowsArray is empty.
            nRows = unique(nRowsArray);

            % NOTE: length==0 valid for struct containing zero numeric fields.
            if length(unique(nRowsArray)) > 1
                error('BICAS:Assertion', ...
                    ['Numeric fields and cell array components', ...
                    ' in struct do not have the same number', ...
                    ' of rows (likely corresponding to CDF zVar records).'])
            end
        end



    end   % Static



end
