function write_CDF_spdfcdfread(filePath, spdfcdfread_out, spdfcdfread_info, varargin)
% write_CDF_spdfcdfread(file_path, spdfcdfread_out, spdfcdfread_info, varargin)   Function which writes a CDF file.
%
% Attempt at a function which can easily write a CDF using variables on the same data format as
% returned by [out, info] = spdfcdfread(... , 'Structure', 1, 'KeepEpochAsIs', 1).
%
% IMPORTANT NOTE: As of 2016-10-20: This function is intended to be phased out and be replaced by write_CDF_dataobj.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-12
%
%
%
% ARGUMENTS
% =========
% spdfcdfread_out  : Data on the format returned from spdfcdfread.
% spdfcdfread_info : Data on the format returned from spdfcdfread.
% varargin         : Optionally 'fill_empty' : EXPERIMENTAL FEATURE.
%                    Empty CDF variable data is replaced by one pad
%                    value scalar (one record only) of the right type.
%                    NOTE: This is not meant for "serious use" (it is only a workaround) but
%                    for tests when writing CDFs, using data read from CDFs. Only works for numerical types.
%
% NOTE: As of 2016-07-28: In practice this function only uses "spdfcdfread_out.Data" in "spdfcdfread_out". Therefore the
% caller does not need to set the other "out" fields, like when reading a master file with empty fields for zero-record
% zVariables.
% 
%
%
% LIMITATIONS
% ===========
% NOTE PROBLEM: spdfcdfread and scpdfcdfinfo may crash MATLAB(!) when reading files written with spdfcdfwrite (which
% this function uses). It appears that this happens when spdfcdfwrite receives various forms of "incomplete" input data.
% spdfcdfwrite appears to often not give any warning/error message when and writes a file anyway. Before passing data to
% spdfcdfwrite, this function tries to give errors for, or correct such data, but can only do so as far as the problem
% is understood by the author. Submitting empty data for a CDF variable is one such case. Therefore, despite best
% efforts, this function might still produce nonsensical files instead of producing any warning or error message.
%
% NOTE PROBLEM(?): Can not select CDF encoding (Network/XDR, IBMPC etc) when writing files. The NASA SPDF MATLAB CDF
% Software distribution does not have this option (this has been confirmed with their email support 2016-07-14).
%
% BUG: SCALEMIN, SCALEMAX for tt2000/Epoch of the wrong type in the final CDF.
%
% BUG/NOTE: spdfcdfwrite has been observed to set the wrong pad value when writing 0 records. Observed when tried to
% copy
%     EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_HK_RPW-BIA_e1d0a9a_CNE_V01.cdf.
% where zVariable TIME_SYNCHRO_FLAG has zero records.
%
% NOTE: spdfcdfwrite always writes char as UCHAR (not CHAR) in the CDF.
%
%
%
% IMPLEMENTATION NOTES: COMPLICATIONS WHEN COMBINING spdfcdfread + cpdfcdfwrite
% =============================================================================
% NOTE: spdfcdfread and spdfcdfwrite interpret 'tt2000' (or at least "Epoch") data differently:
% - spdfcdfread by default converts tt2000 data to MATLAB's own time scalar (using datenum)
%   datestr converts value_default-->UTC. ("value_default" refers to the above).
%   Use spdfcdfread argument 'KeepEpochAsIs'=1 to read data that can be directly passed on to spdfcdfwrite.
%   spdfbreakdowntt2000 converts value_KeepEpochAsIs-->UTC.
% - spdfdatenumtott2000 converts value_default-->value_KeepEpochAsIs
%   However, in practice there is a slight difference in data (rounding?) between
%    (1) Read the CDF and get MATLAB's time scalar, then convert to tt2000, and
%    (2) Read the CDF with spdfcdfread+KeepEpochAsIs=1.
%
% NOTE: FILLVAL, VALIDMIN/-MAX, SCALEMIN/-MAX and pad values should presumably be thought of as having the
% same CDF data type (or MATLAB class) as the corresponding variable data. However, scpdfcdfread returns
% 1) FILLVAL, VALIDMIN/-MAX, SCALEMIN/-MAX as UTC strings for tt2000 data (or at least "Epoch").
% 2) pad value (info.Variables(i, 9)       as a tt2000 value (independent of spdfcdfread+KeepEpochAsIs flag).
%
% NOTE: String CDF variables can not be written back in the same format. spdfcdfread returns (0-dimensional and
% 1-dimensional) CDF char variables as matrices with
%   index 1=row   =component within a 1-dimensional record
%   index 2=column=#character in string (as usual)
%   index 3=#record.
% Empirically, strings should be specified without the (spdfcdfwrite) 'RecordBound' option. Multiple strings within a
% (1-dimensional) record should be specified as one string per row. Every such matrix is a component in a column cell array.
%
% NOTE: [out, info] = spdfcdfread(..., 'Structure', 1) returns out(i).VariableName/.Data/.Attributes==[] for zero-record
% zVariables. info.VariableAttributes.FIELDNAM{i,1} and info.Variables are still complete in those cases. Therefore, use
% info.Variables(:,1) instead. This is important when reading master CDFs which by their nature contain many zero-record
% zVariables.
%
% NOTE: The format of data returned by
%    [out, info] = spdfcdfread(... , 'Structure', 1, 'KeepEpochAsIs', 1)
% is peculiar. 
%  - It contains redundant data.
%  - The "out" array contains empty components (the struct fields are empty), apparently for CDF variables which contain no data.
%    However, the CDF variable names and variable attributes are still represented in info.Variables{#variable, #column} and
%    info.VariableAttributes.(var_attr_name) == [var_names_column, values_column].
%    Note: It appears that out(i) corresponds to info.Variables(i,:).
%



% PROPOSAL: Enable MD5 checksum.
% PROPOSAL: Option for filling empty variable with pad values. 1 record?
%    CON: Can not know the (non-record) dimensions of such a CDF variable.
%    CON: Using exactly one record automatically leads to the CDF labelling the CDF variable as record-invariant!
% PROPOSAL: Redo to work easily with dataobj instead? How handles time? Separat function? One wraps the other?
%    NOTE: Not necessarily difficult since dataobj seems to contain structures/arrays from spdfcdfread/-info.
%    QUESTION: Can dataobj be altered?!
%
% PROPOSAL: Call NASA SPDFs Java code instead?!! Should be possible to easily call Java code from inside MATLAB.
% PROPOSAL: Redefine argument "spdfcdfread_out" to "spdfcdfread_out_Data".
%   NOTE: Ambiguous format. Cell array (vector) of cell arrays?
%   CON: Not really analogue with anything that can be directly extracted from spdfcdfread.
%
% PROPOSAL: Some form of validation of input.
%    PRO: Would confirm that the input format is as it is understood.
%    PROPOSAL: Check that out(i) corresponds to info.Variables(i,:) (where "out" fields are not empty).
%

    VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE = {'VALIDMIN', 'VALIDMAX', 'SCALEMIN', 'SCALEMAX', 'FILLVAL'};
        
    % Rename variables for convenience.
    out  = spdfcdfread_out;
    info = spdfcdfread_info;



    % Shorten? Unnecessarily complicated if there is only one flag.
    fillEmptyVariables = 0;
    for i=1:length(varargin)
        if strcmp(varargin{i}, 'fill_empty')
            fillEmptyVariables = 1;
        else
            error('write_CDF:Assertion:IllegalArgument', 'Can not interpret argument option.')
        end
    end
    
    
    
    %============================================================
    % Construct variables that spdfcdfwrite accepts as arguments
    %============================================================
    VARIABLE_LIST = {};
    RECBNDVARS    = {};
    VARDATATYPE   = {};
    PADVALS       = {};
    
    for i=1:length(out)
        %---------------------------------------------------------------------------------------------------------------
        % IMPLEMENTATION NOTE: Not using (1) out(i).VariableName or (2) info.Variables(:,1) to obtain the variable name
        % since experience shows that components of (1) can be empty (contain empty struct fields) and (2) may not cover
        % all variables when obtained via spdfcdfread!!
        %---------------------------------------------------------------------------------------------------------------
        %zVarName             = info.VariableAttributes.FIELDNAM{i, 1};    % i.Variables{i, 1}
        zVarName             = info.Variables{i, 1};
        spdfcdfwriteDataType = info.Variables{i, 4};   % uint32, tt2000 etc. Change name to CDF_Data_type?! spdfcdfwrite_data_type?!
        padValue             = info.Variables{i, 9};
        
        %=========================================
        % Special case for zero record zVariables
        %=========================================
        if isempty(out(i).Data) % && ~isempty(out(i).VariableName) && ~isempty(out(i).Attributes)
            if ~fillEmptyVariables
                error('write_CDF:Assertion', 'Can not handle CDF variables with zero records (due to presumed bug in spdfcdfwrite).')
            else
                %----------------------------------------------------------------------------------------
                % EXPERIMENTAL SOLUTION: Store data with pad values instead.
                % NOTE: Incomplete since does not take CDF variable type, array dimensions into account.
                %----------------------------------------------------------------------------------------
                %nRecords = max([s{:}]);   % Get the greatest number of records that any variable has.
                nRecords = 1;
                matlabClass = bicas.utils.convert_CDF_type_to_MATLAB_class(spdfcdfwriteDataType, 'Permit MATLAB classes');
                try
                    zVarData = cast(ones(nRecords, 1), matlabClass) * padValue;       % NOTE: Hardcoded uint8
                catch exception
                    error('write_CDF:Assertion', 'Can not type cast variable to "%s" (CDF: "%s").', matlabClass, spdfcdfwriteDataType)
                end
            end
        else
            zVarData = out(i).Data;
        end
        
        
        
        %========================================================================================================
        % Check that the supplied zVariable data has a MATLAB class (type) which matches the specified CDF type
        % -----------------------------------------------------------------------------------------------------
        % IMPLEMENTATION NOTE:
        % (1) Empty data (empty arrays) from spdfcdfread are known to have the wrong data type (char).
        % Therefore, do this check after having dealt with empty data.
        % (2) Must do this after converting time strings (char) data to uint64/tt2000.
        %========================================================================================================
        matlabClass = class(zVarData);
        if ~strcmp(  bicas.utils.convert_CDF_type_to_MATLAB_class(spdfcdfwriteDataType, 'Permit MATLAB classes'),   matlabClass  )
            error('write_CDF:Assertion', 'MATLAB data type "%s" does not match the actual CDF data ("%s") for CDF variable "%s".', matlabClass, spdfcdfwriteDataType, zVarName)
        end
        
        
        
        %==========================================================================================
        % Special case for "char"
        % -----------------------
        % NOTE: Previous tt2000/UTC strings have already been converted to non-char at this stage.
        % QUESTION 2016-10-20: Is above true? Obsolete comment?
        %==========================================================================================
        if ischar(zVarData)
            %=======================================================================
            % Convert 3-D char matrices to column cell arrays of 2-D char matrices.
            %=======================================================================
            if ndims(zVarData) > 3
                error('write_CDF:Assertion:OperationNotImplemented', 'Can not handle more CDF char strings variables with more than 1 dimension (excluding the record dimension).')
            end
            zVarDataNew = {};
            for iRecord = 1:size(zVarData, 3)
                zVarDataNew{end+1,1} = zVarData(:, :, iRecord);
            end
            zVarData = zVarDataNew;
        else
            RECBNDVARS{end+1} = zVarName;            
        end
        

        
        %=================================================================
        % Convert specific VariableAttributes values.
        % 1) tt2000 values as UTC strings: converted to tt2000.
        % 2) All other: Convert to the (z)variable data type.
        %
        % IMPLEMENTATION NOTE: spdfcdfread (not spdfcdfwrite) can crash if not doing this!!! Think
        % it is tt2000 CDF variables which are the problem(?).
        %
        % BUG: Does not seem to work on SCALEMIN/-MAX specifically despite identical treatment, for unknown reason.
        %=================================================================
        for i_va = 1:length(VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE)
            varAttrName = VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE{i_va};
            if ~isfield(info.VariableAttributes, varAttrName)
                continue
            end
            
            % IMPLEMENTATION NOTE: Can NOT assume that every CDF variable is represented among
            % info.VariableAttributes.
            % Example: EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-CWF_e1d0a9a_CNE_V01.cdf
            iVAf = info.VariableAttributes.(varAttrName);
            i_va = find(strcmp(iVAf(:,1), zVarName));
            if length(i_va) == 0
                continue
            elseif length(i_va) > 1
                error('write_CDF:Assertion:OperationNotImplemented', 'Can not handle multiple variable name matches in info.VariableAttributes.%s.', varAttrName)
            end
            varAttrValue = iVAf{i_va, 2};
            if strcmp(spdfcdfwriteDataType, 'tt2000') && ischar(varAttrValue) && 1
                varAttrValue = spdfparsett2000(varAttrValue);
            else
                if ~strcmp(spdfcdfwriteDataType, class(varAttrValue))
                    error('write_CDF:Assertion', 'Found VariableAttribute %s for CDF variable %s whose data type did not match the declared one.', varAttrName, zVarName)
                end
                %fprintf('%s  %s  ', class(var_attr_value), spdfcdfwrite_data_type)   % DEBUG
                %var_attr_value = typecast(var_attr_value, spdfcdfwrite_data_type);                
                %fprintf('%s\n', class(var_attr_value))   % DEBUG                
            end
            iVAf{i_va, 2} = varAttrValue;
            info.VariableAttributes.(varAttrName) = iVAf;
        end

        VARIABLE_LIST(end+[1,2]) = {zVarName, zVarData              };        
        VARDATATYPE  (end+[1,2]) = {zVarName, spdfcdfwriteDataType};
        PADVALS      (end+[1,2]) = {zVarName, padValue              };
    end
    
    % info.VariableAttributes = rmfield(info.VariableAttributes, {'VALIDMIN', 'VALIDMAX', 'SCALEMIN', 'SCALEMAX', 'FILLVAL'});
    
    
%===================================================================================================
% RELEVANT spdfcdfwrite OPTIONS:
% (Relevant excerpts from spdfcdfwrite.m COPIED here for convenience.)
% --------------------------------------------------------------------
%   SPDFCDFWRITE(FILE, VARIABLELIST, ...) writes out a CDF file whose name
%   is specified by FILE.  VARIABLELIST is a cell array of ordered
%   pairs, which are comprised of a CDF variable name (a string) and
%   the corresponding CDF variable value.  To write out multiple records
%   for a variable, there are two ways of doing it. One way is putting the
%   variable values in a cell array, where each element in the cell array
%   represents a record. Another way, the better one, is to place the
%   values in an array (single or multi-dimensional) with the option
%   'RecordBound' being specified.  
%
%   SPDFCDFWRITE(..., 'RecordBound', RECBNDVARS) specifies data values in arrays
%   (1-D or multi-dimensional) are to be written into "records" for the given
%   variable. RECBNDVARS is a cell array of variable names. The M-by-N array
%   data will create M rows (records), while each row having N elements. For
%   examples, 5-by-1 array will create five (5) scalar records and 1-by-5 array
%   will write out just one (1) record with 5 elements. For 3-D array of
%   M-by-N-by-R, R records will be written, and each record with M-by-N
%   elements. Without this option, array of M-by-N will be written into a single
%   record of 2-dimensions. See sample codes for its usage.
%
%   SPDFCDFWRITE(..., 'GlobalAttributes', GATTRIB, ...) writes the structure
%   GATTRIB as global meta-data for the CDF.  Each field of the
%   struct is the name of a global attribute.  The value of each
%   field contains the value of the attribute.  To write out
%   multiple values for an attribute, the field value should be a
%   cell array. 
%
%   If there is a master CDF that has all the meta-data that the new CDF needs,
%   then SPDFCDFINFO module can be used to retrieve the infomation. The
%   'GlobalAttributes' field from the returned structure can be
%   passed in for the GATTRIB.
%
%   In order to specify a global attribute name that is illegal in
%   MATLAB, create a field called "CDFAttributeRename" in the 
%   attribute struct.  The "CDFAttribute Rename" field must have a value
%   which is a cell array of ordered pairs.  The ordered pair consists
%   of the name of the original attribute, as listed in the 
%   GlobalAttributes struct and the corresponding name of the attribute
%   to be written to the CDF.
%
%   SPDFCDFWRITE(..., 'VariableAttributes', VATTRIB, ...) writes the
%   structure VATTRIB as variable meta-data for the CDF.  Each
%   field of the struct is the name of a variable attribute.  The
%   value of each field should be an Mx2 cell array where M is the
%   number of variables with attributes.  The first element in the
%   cell array should be the name of the variable and the second
%   element should be the value of the attribute for that variable.
%
%   If there is a master CDF that has all the meta-data that the new CDF needs,
%   then SPDFCDFINFO module can be used to retrieve the infomation. The 
%   'VariableAttributes' field from the returned structure can
%   be passed in for the VATTRIB.
%
%   In order to specify a variable attribute name that is illegal in
%   MATLAB, create a field called "CDFAttributeRename" in the 
%   attribute struct.  The "CDFAttribute Rename" field must have a value
%   which is a cell array of ordered pairs.  The ordered pair consists
%   of the name of the original attribute, as listed in the 
%   VariableAttributes struct and the corresponding name of the attribute
%   to be written to the CDF.   If you are specifying a variable attribute
%   of a CDF variable that you are re-naming, the name of the variable in
%   the VariableAttributes struct must be the same as the re-named variable.
%
%   SPDFCDFWRITE(..., 'Vardatatypes', VARDATATYPE) specifies the variable's
%   data types. By default, this module uses each variable's passed data to
%   determine its corresponding CDF data type. While it is fine for the most
%   cases, this will not work for the CDF epoch types, i.e., CDF_EPOCH (a double),
%   CDF_EPOCH16 (an array of 2 doubles) and CDF_TIME_TT2000 (an int64). This
%   option can be used to address such issue. VARDATATYPE is a cell array of
%   variable names and their respective data types (in string).
%
%   The following table shows the valid type strings, either in CDF defined
%   forms, or alternatively in the forms presented at column 4 in the Variables
%   field of the structure returned from a SPDFCDFINFO module call to an
%   existing CDF or master CDF.  
%       type             CDF Types
%       -----            ---------
%       int8             CDF_INT1 or CDF_BYTE
%       int16            CDF_INT2
%       int32            CDF_INT4
%       int64            CDF_INT8
%       uint8            CDF_UINT1
%       uint16           CDF_UINT2
%       uint32           CDF_UINT4
%       single           CDF_FLOAT or CDF_REAL4
%       double           CDF_DOUBLE or CDF_REAL8
%       epoch            CDF_EPOCH
%       epoch16          CDF_EPOCH16
%       tt2000           CDF_TIME_TT2000
%       char             CDF_CHAR or CDF_UCHAR
%
%   Note: Make sure variable's data match to the defined type.
%
%   SPDFCDFWRITE(..., 'PadValues', PADVALS) writes out pad values for given
%   variable names.  PADVALS is a cell array of ordered pairs, which
%   are comprised of a variable name (a string) and a corresponding 
%   pad value.  Pad values are the default value associated with the
%   variable when an out-of-bounds record is accessed.  Variable names
%   that appear in PADVALS must appear in VARIABLELIST.
%===================================================================================================
irf.log('n', sprintf('Writing CDF file "%s".', filePath))
spdfcdfwrite(filePath, VARIABLE_LIST(:), 'RecordBound', RECBNDVARS, 'GlobalAttributes', info.GlobalAttributes, ...
    'VariableAttributes', info.VariableAttributes, 'Vardatatypes', VARDATATYPE, 'PadValues', PADVALS)

end
