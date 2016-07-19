% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-12
%
% Attempt at a function which can easily write a CDF using variables on the same data format as
% returned by spdfcdfread(... , 'Structure', 1, 'KeepEpochAsIs', 1).
%
%
% PROBLEM: spdfcdfread and scpdfcdfinfo may crash MATLAB(!) when reading files written with
% spdfcdfwrite (which this function uses). It appears that this happens when spdfcdfwrite receives
% various forms of "incomplete" input data. spdfcdfwrite appears to be very bad at giving direct
% errors for such data and writes the file anyway. Before passing data to spdfcdfwrite, this
% function tries to give errors for, or correct such data, but can only do so as far as the problem
% is understood. Submitting empty data for a CDF variable is one such case.
%
% PROBLEM(?): Can not select CDF encoding (Network/XDR, IBMPC etc) when writing files. NASA SPDF
% MATLAB CDF Software distribution does not have this option (confirmed with their email support
% 2016-07-14).
%
% BUG: SCALEMIN, SCALEMAX for tt2000/Epoch of the wrong type in the final cdf.
%
% NOTE: spdfcdfwrite has been observed to set the wrong pad value when writing 0 records. 
% Observed when tried to copy 
% EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_HK_RPW-BIA_e1d0a9a_CNE_V01.cdf.
% where zVariable TIME_SYNCHRO_FLAG has zero records.
%
% NOTE: spdfcdfwrite always writes char as UCHAR (not CHAR).
%
%
%
% COMPLICATIONS WHEN COMBINING spdfcdfread + cpdfcdfwrite
% -------------------------------------------------------
% NOTE: spdfcdfread and spdfcdfwrite interpret 'tt2000' (or at least "Epoch") data differently:
% - spdfcdfread by default converts tt2000 data to MATLAB's own time scalar (using datenum).
%   datestr converts value_default-->UTC. 
%   Use spdfcdfread argument 'KeepEpochAsIs'=1 to read data that can be directly passed on to spdfcdfwrite.
%   spdfbreakdowntt2000 converts value_KeepEpochAsIs-->UTC.
% - spdfdatenumtott2000 converts value_default-->value_KeepEpochAsIs
%   However, in practice there is a slight difference in data (rounding?) between
%    (1) Read the cdf and get MATLAB's time scalar, and convert to tt2000, and
%    (2) Read the cdf with spdfcdfread+KeepEpochAsIs=1.
%
% NOTE: FILLVAL, VALIDMIN/-MAX, SCALEMIN/-MAX and pad values should presumably be thought of as the
% same type as the corresponding variable data. However, scpdfcdfread returns
% 1) FILLVAL, VALIDMIN/-MAX, SCALEMIN/-MAX as UTC strings for tt2000 data (or at least "Epoch").
% 2) pad value (info.Variables(i, 9) as a tt2000 value (independent of spdfcdfread+KeepEpochAsIs flag).
%
% NOTE: [out, info] = spdfcdfread(..., 'Structure', 1) has been observed to return out(3).VariableName
% == [] despite that info.VariableAttributes.FIELDNAM{3,1} implied that there was a name which could
% be used in its place (and had to be used for the purpose of "copying" a master cdf).
% Use info.Variables(:,1) instead.
% Observed in file: EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_HK_RPW-BIA_e1d0a9a_CNE_V01.cdf.
%
% NOTE: String CDF variables can not be written back in the same format. spdfcdfread returns (0-D
% and 1-D) CDF char variables as matrices with index=#character, index 2=row=component within a
% 1-D record, index 3=record.
% Empirically, strings should be specified without the 'RecordBound' option. Multiple strings within
% a (1-D) record should be specified as one string per row. Every such matrix is a component in a
% column cell array.
%
% NOTE: 
%
%
%
% ARGUMENTS
% ---------
% varargin: Optionally 'fill_empty' : Empty CDF variable data is replaced by one pad value scalar
% (one record only) of the right type. This is not meant for "serious use" but for tests when
% writing CDFs, using data read from CDFs. Only works for numerical types.
%
function write_cdf(file_path, spdfcdfread_out, spdfcdfread_info, varargin)
%
% PROPOSAL: Check/assert that CDF variable DATA has the right type.
% PROPOSAL: Enable MD5 checksum.
% PROPOSAL: Option for filling empty variable with pad values. 1 record?
%    CON: Can not know the (non-record) dimensions of such a CDF variable.
%    CON: Using exactly one record automatically leads to the CDF labelling the CDF variable as record-invariant!
% PROPOSAL: Redo to work easily with dataobj instead? How handles time? Separat function? One wraps the other?
%
% PROPOSAL: Switch from type casting data to checking formats, with exceptions.
% PROPOSAL: Call NASA SPDFs Java code instead?!! Should be possible to easily call Java code from inside MATLAB.
%
    VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE = {'VALIDMIN', 'VALIDMAX', 'SCALEMIN', 'SCALEMAX', 'FILLVAL'};
        
    out  = spdfcdfread_out;
    info = spdfcdfread_info;
    
    fill_empty_variables = 0;
    for i=1:length(varargin)
        if strcmp(varargin{i}, 'fill_empty')
            fill_empty_variables = 1;
        else
            error('Can not interpret argument option.')
        end
    end
    
    
    
    %============================================================
    % Construct variables that spdfcdfwrite accepts as arguments
    %============================================================
    VARIABLE_LIST = {};
    RECBNDVARS = {};
    VARDATATYPE = {};
    PADVALS = {};
    
    for i=1:length(out)
        % IMPLEMENTATION NOTE: Not using (1) out(i).VariableName or (2) info.Variables(:,1) to
        % obtain the variable name since experience shows that (1) can be empty and (2) may not
        % cover all variables when obtained via spdfcdfread!!
        % TODO: Find examples!!! Assertions?
        %var_name  = info.VariableAttributes.FIELDNAM{i, 1};    % i.Variables{i, 1}
        var_name  = info.Variables{i, 1};
        spdfcdfwrite_data_type = info.Variables{i, 4};   % uint32, tt2000 etc. Change name to cdf_Data_type?! spdfcdfwrite_data_type?!
        pad_value = info.Variables{i, 9};
        
        if isempty(out(i).Data)            
            if ~fill_empty_variables
                error('Can not handle CDF variables with zero records (due to presumed bug in psdfcdfwrite).')
            else
                % EXPERIMENTAL SOLUTION: Store data with pad values instead.
                % NOTE: Incomplete since does not take CDF variable type, array dimensions into account.

                %s = info.Variables(:,3);   % Get number of records for all variables.
                %N_records = max([s{:}]);   % Get the greatest number of records that any variable has.
                N_records = 1;
                data = cast(ones(N_records, 1), spdfcdfwrite_data_type) * pad_value;       % NOTE: Hardcoded uint8
            end
        else
            data = out(i).Data;
        end
        
        % IMPLEMENTATION NOTE:
        % (1) Empty data (empty arrays) from spdfcdfread are known to have the wrong data type (char).
        % Therefore, do this check after having dealt with empty data.
        % (2) Must do this after converting tt2000 char data to uint64/tt2000.        
        data_data_type = class(data);
        if ~MATLAB_type_compatible_with_cdf_type(data_data_type, spdfcdfwrite_data_type)
            error(sprintf('MATLAB data type "%s" does not match the actual CDF data ("%s") for CDF variable "%s".', data_data_type, spdfcdfwrite_data_type, var_name))
        end
        
        % NOTE: Previous tt2000 strings have already been converted to non-char at this stage.
        if ischar(data)
            % Convert 3-D char matrices to column cell arrays of 2-D char matrices.
            if length(size(data)) > 3
                error('Can not handle more CDF char strings variables with more than 1 dimension (excluding the record dimension).')
            end
            data_new = {};
            for i_record = 1:size(data, 3)
                data_new{end+1,1} = data(:,:, i_record);
            end
            data = data_new;
        else
            RECBNDVARS{end+1}        = var_name;
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
            var_attr_name = VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE{i_va};
            if ~isfield(info.VariableAttributes, var_attr_name)
                continue
            end
            
            % IMPLEMENTATION NOTE: Can NOT assume that every CDF variable is represented among
            % info.VariableAttributes.
            % Example: EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-CWF_e1d0a9a_CNE_V01.cdf
            iVAf = info.VariableAttributes.(var_attr_name);
            i_va = find(strcmp(iVAf(:,1), var_name));
            if length(i_va) == 0
                continue
            elseif length(i_va) > 1
                error(sprintf('Can not handle multiple variable name matches in info.VariableAttributes.%s.', var_attr_name))
            end
            var_attr_value = iVAf{i_va, 2};
            if strcmp(spdfcdfwrite_data_type, 'tt2000') && ischar(var_attr_value) && 1
                var_attr_value = spdfparsett2000(var_attr_value);
            else
                if ~strcmp(spdfcdfwrite_data_type, class(var_attr_value))
                    error(sprintf('Found VariableAttribute %s for CDF variable %s whose data type did not match the declared one.', var_attr_name, var_name))
                end
                %fprintf('%s  %s  ', class(var_attr_value), spdfcdfwrite_data_type)   % DEBUG
                %var_attr_value = typecast(var_attr_value, spdfcdfwrite_data_type);                
                %fprintf('%s\n', class(var_attr_value))   % DEBUG                
            end
            iVAf{i_va, 2} = var_attr_value;
            info.VariableAttributes.(var_attr_name) = iVAf;
        end

        VARIABLE_LIST(end+[1,2]) = {var_name, data};        
        VARDATATYPE(end+[1,2])   = {var_name, spdfcdfwrite_data_type};
        PADVALS(end+[1,2])       = {var_name, pad_value};
    end
    
    % info.VariableAttributes = rmfield(info.VariableAttributes, {'VALIDMIN', 'VALIDMAX', 'SCALEMIN', 'SCALEMAX', 'FILLVAL'});
    
    
%===================================================================================================
% Relevant spdfcdfwrite options:
% ------------------------------
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
irf.log('n', sprintf('Writing CDF file "%s".', file_path))
spdfcdfwrite(file_path, VARIABLE_LIST(:), 'RecordBound', RECBNDVARS, 'GlobalAttributes', info.GlobalAttributes, ...
    'VariableAttributes', info.VariableAttributes, 'Vardatatypes', VARDATATYPE, 'PadValues', PADVALS)

end



% Checks whether a MATLAB type is compatible with a stated type that spdfcdfwrite accepts.
% NOTE: Not implemented for all types yet.
function result = MATLAB_type_compatible_with_cdf_type(MATLAB_type, cdf_type)

% Left column = Legal MATLAB types.
% Right column(s) = Legal types for scpdfcdfwrite which correspond to the ONE MATLAB type in the left column.
DATA = {...
    'int8',    {'int8',   'CDF_INT1', 'CDF_BYTE'};
    'int16',   {'int16',  'CDF_INT2'}; ...
    'int32',   {'int32',  'CDF_INT4'}; ...
    'int64',   {'int64',  'CDF_INT8', 'tt2000', 'CDF_TIME_TT2000'}; ...
    'uint8',   {'uint8',  'CDF_UINT1'}; ...
    'uint16',  {'uint16', 'CDF_UINT2'}; ...
    'uint32',  {'uint32', 'CDF_UINT4'}; ...
    'single',  {'single', 'CDF_FLOAT',  'CDF_REAL4'}; ...
    'double',  {'double', 'CDF_DOUBLE', 'CDF_REAL8'}; ...
    'char',    {'char',   'CDF_CHAR',   'CDF_UCHAR'}; ...
    };
%    'epoch',   {'CDF_EPOCH'}; ...
%    'epoch16', {'CDF_EPOCH16'}; ...

i = find(strcmp(DATA(:,1), MATLAB_type));
if numel(i) ~= 1
    error('Can not identify MATLAB type "%s".', MATLAB_type)
else
    result = any(strcmp(DATA{i,2}, cdf_type));
end

end