function write_CDF_dataobj(file_path, dataobj_GlobalAttributes, dataobj_data, dataobj_VariableAttributes, dataobj_Variables, varargin)
% write_CDF_dataobj(file_path, spdfcdfread_data, spdfcdfread_info, varargin)   Function which writes a CDF file.
%
% Attempt at a function which can easily write a CDF using variables on the same data format as returned by dataobj
% (irfu-matlab). Useful for reading a CDF file, modifying the contents somewhat, and then writing the modified contents
% to a CDF file. Originally based on write_cdf.m/write_cdf_spdfcdfread.m.
%
%
% IMPORTANT NOTE: As of 2016-10-20: This function is intended to replace write_CDF_spdfcdfread.
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2016-07-12 (write_cdf.m/write_cdf_spdfcdfread.m), 2016-10-20 (write_cdf_dataobj.m)
%
%
% ARGUMENTS
% =========
% file_path                   : Path to file to create.
% dataobj_GlobalAttributes,   
% dataobj_data,
% dataobj_VariableAttributes,
% dataobj_Variables           : The corresponding fields of an instantiated dataobj, i.e. do.data,
%                               do.GlobalAttributes and so on, where do=dataobj(...).
% varargin                    : EXPERIMENTAL FEATURE. Optionally string 'fill_empty':
%                               Empty CDF variable data is replaced by one pad
%                               value scalar (one record only) of the right type.
%                               NOTE: This is not meant for "serious use" (it is only a workaround) but
%                               for tests when writing CDFs, using data read from CDFs. Only works for numerical types.
%
% NOTE: 
% 
%
%
% LIMITATIONS
% ===========
% NOTE/PROBLEM: spdfcdfread and scpdfcdfinfo may crash MATLAB(!) when reading files written with spdfcdfwrite (which
% this function uses). It appears that this happens when spdfcdfwrite receives various forms of "incomplete" input data.
% spdfcdfwrite appears to often not give any warning/error message when receiving such data and writes a file anyway
% with neither error nor warning. Before passing data to spdfcdfwrite, this function tries to give errors for, or
% correct such data, but can only do so as far as the problem is understood by the author. Submitting empty data for a
% CDF variable is one such case. Therefore, despite best efforts, this function might still produce nonsensical files
% instead of producing any warning or error message.
%
% NOTE PROBLEM(?): Can not select CDF encoding (Network/XDR, IBMPC etc) when writing files. The NASA SPDF MATLAB CDF
% Software distribution does not have this option (this has been confirmed with their email support 2016-07-14).
%
% BUG: Variable attributes SCALEMIN, SCALEMAX for Epoch are stored as CDF_INT8 (not CDF_TT2000) in the final CDF file.
% The information stored seems correct though. Therefore, the same variable attributes are also represented as integers
% when reading the CDF file with dataobj.
%
% BUG/NOTE: spdfcdfwrite has been observed to set the wrong pad value when writing 0 records. Observed when tried to
% copy
%     EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_HK_RPW-BIA_e1d0a9a_CNE_V01.cdf.
% where zVariable TIME_SYNCHRO_FLAG has zero records.
%
% NOTE: spdfcdfwrite always writes char as UCHAR (not CHAR) in the CDF.
%
% NOTE: Might fail for matrices of four dimensions or higher (including records as one dimension). Untested.
%
%
%
% IMPLEMENTATION NOTE
% ===================
% The function does not accept a whole dataobj object since:
% (1) instances of dataobj are likely not meant to be modified after creation. It is possible to modify them though.
% Therefore the function only accepts the parts of a dataobj that it really needs, which still means it accepts a lot of
% redundant information in the arguments.
% (2) One might want to use the function for writing a CDF file without basing it on a dataobj (i.e. without basing it
% on an existing CDF file).

%=======================================================================================================================
% PROPOSAL: Implement using NASA SPDFs Java code instead?!! Should be possible to easily call Java code from inside MATLAB.
%   PRO?: Java interface might be more easy to work with and have fewer quirks & limitations.
%
% PROPOSAL: Enable MD5 checksum.
%
% PROPOSAL: Option for filling empty variable with pad values. 1 record?
%    CON: Can not know the (non-record) dimensions of such a CDF variable.
%    CON: Using exactly one record automatically leads to the CDF labelling the CDF variable as record-invariant!
%
% PROPOSAL: Rework into function that uses write_CDF_spdfcdfread (wrapper)?!! 
%    NOTE: Not necessarily difficult since dataobj seems to contain structures/arrays from spdfcdfread/-info.
%    PRO: Would reduce amount of code.
%    CON: write_CDF_spdfread needs knowledge of strange behaviour in spdfcdfread to mimic it (int the interface), which
%         this function then needs to have knowledge of to use write_CDF_spdfread.
%    PROPOSAL: Abolish write_CDF_spdfcdfread, or at least stop updating it.
%
% PROPOSAL: Some form of validation of input.
%    PRO: Would confirm that the input format is as it is understood.
%    PROPOSAL: Check that stated size (within record) fits data.
%       PROPOSAL: Flag
%    PROPOSAL: Check that stated nbr of records fits data.
%       PROPOSAL: Flag
%    PROPOSAL: Assertions for redundant data (within dataobj data/attributes).
%
% PROPOSAL: Accept a smaller subset of dataobj data, with less/no redundant data.



    VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE = {'VALIDMIN', 'VALIDMAX', 'SCALEMIN', 'SCALEMAX', 'FILLVAL'};



    % Shorten? Unnecessarily complicated if there is only one flag.
    fill_empty_variables = 0;
    for i=1:length(varargin)
        if strcmp(varargin{i}, 'fill_empty')
            fill_empty_variables = 1;
        else
            error('write_CDF_dataobj:Assertion:IllegalArgument', 'Can not interpret argument option.')
        end
    end
    
    
    
    %============================================================
    % Construct variables that spdfcdfwrite accepts as arguments
    %============================================================
    VARIABLE_LIST = {};
    RECBNDVARS = {};
    VARDATATYPE = {};
    PADVALS = {};
    
    %for i=1:length(data)
    for i=1:length(dataobj_Variables(:,1))
        %---------------------------------------------------------------------------------------------------------------
        % IMPLEMENTATION NOTE: Not using (1) data(i).VariableName or (2) info.Variables(:,1) to obtain the variable name
        % since experience shows that components of (1) can be empty (contain empty struct fields) and (2) may not cover
        % all variables when obtained via spdfcdfread!!
        %---------------------------------------------------------------------------------------------------------------
        zVar_name              = dataobj_Variables{i, 1};
        dataobj_stated_data_type = dataobj_Variables{i, 4};   % uint32, tt2000 etc. Change name to CDF_Data_type?! dataobj_stated_data_type?!
        pad_value              = dataobj_Variables{i, 9};   % This value can NOT be found in dataobj_data. Has to be read from dataobj_Variables.
        
        %=========================================
        % Special case for zero record zVariables
        %=========================================
        if isempty(dataobj_data.(zVar_name).data)
            if ~fill_empty_variables
                error('write_CDF_dataobj:Assertion', 'Can not handle CDF variables with zero records (due to presumed bug in spdfcdfwrite).')
            else
                %----------------------------------------------------------------------------------------
                % EXPERIMENTAL SOLUTION: Store data with pad values instead.
                % NOTE: Incomplete since does not take CDF variable type, array dimensions into account.
                %----------------------------------------------------------------------------------------
                %N_records = max([s{:}]);   % Get the greatest number of records that any variable has.
                N_records = 1;
                zVar_data_MATLAB_class = bicas.utils.convert_CDF_type_to_MATLAB_class(dataobj_stated_data_type, 'Permit MATLAB classes');
                try
                    zVar_data = cast(ones(N_records, 1), zVar_data_MATLAB_class) * pad_value;       % NOTE: Hardcoded uint8
                catch exception
                    error('write_CDF_dataobj:Assertion', 'Can not type cast variable to "%s" (CDF: "%s").', ...
                        zVar_data_MATLAB_class, dataobj_stated_data_type)
                end
            end
        else
            zVar_data = dataobj_data.(zVar_name).data;
        end
        
        
        
        %========================================================================================================
        % Check that the supplied zVariable data has a MATLAB class (type) which matches the specified CDF type
        % -----------------------------------------------------------------------------------------------------
        % IMPLEMENTATION NOTE:
        % (1) Empty data (empty arrays) from spdfcdfread are known to have the wrong data type (char).
        % Therefore, do this check after having dealt with empty data.
        % (2) Must do this after converting time strings (char) data to uint64/tt2000.
        %========================================================================================================
        zVar_data_MATLAB_class = class(zVar_data);
        if ~strcmp(  bicas.utils.convert_CDF_type_to_MATLAB_class(dataobj_stated_data_type, 'Permit MATLAB classes'),   zVar_data_MATLAB_class  )
            error('write_CDF_dataobj:Assertion', ...
                'MATLAB class (variable type) ("%s") does not match the actual CDF data type ("%s") for CDF variable "%s".', ...
                zVar_data_MATLAB_class, dataobj_stated_data_type, zVar_name)
        end
        
        
        
        if ischar(zVar_data)
            %===========================================================================================================
            % Special case for "char": Convert 3-D char matrices to column cell arrays of 2-D char matrices.
            % -----------------------
            % IMPLEMENTATION NOTE: It is not possible to permute indices for string as one can for non-char for ndim==3.
            %===========================================================================================================
            if ndims(zVar_data) > 3
                error('write_CDF_dataobj:Assertion:OperationNotImplemented', ...
                    'Can not handle more CDF char strings variables with more than 1 dimension (excluding the record dimension).')
            end
            zVar_data_new = {};
            for i_record = 1:size(zVar_data, 3)
                zVar_data_new{end+1,1} = zVar_data(:,:, i_record);
            end
            zVar_data = zVar_data_new;
        else
            RECBNDVARS{end+1} = zVar_name;
            
            %===========================================================================================================
            % Special behaviour for 3D matrices.
            % ----------------------------------
            % For 3D matrices, spdfcdfwrite interprets the last index (not the first index!) as the record number.
            % Must therefore permute the indices so that write_cdf2 is consistent for all numbers of dimensions.
            %     write_cdf2 data arguments        : index 1 = record.
            %     matrix passed on to spdfcdfwrite : index 3 = record.
            % NOTE: spdfcdfread (at least with argument "'Structure', 1, 'KeepEpochAsIs', 1") works like spdfcdfwrite in
            % this regard. 
            %
            % Excerpt from the comments in "spdfcdfwrite.m":
            % ----------------------------------------------
            %   """SPDFCDFWRITE(..., 'RecordBound', RECBNDVARS) specifies data values in arrays
            %   (1-D or multi-dimensional) are to be written into "records" for the given
            %   variable. RECBNDVARS is a cell array of variable names. The M-by-N array
            %   data will create M rows (records), while each row having N elements. For
            %   examples, 5-by-1 array will create five (5) scalar records and 1-by-5 array
            %   will write out just one (1) record with 5 elements. For 3-D array of
            %   M-by-N-by-R, R records will be written, and each record with M-by-N
            %   elements. Without this option, array of M-by-N will be written into a single
            %   record of 2-dimensions. See sample codes for its usage."""
            %===========================================================================================================
            if ndims(zVar_data) == 3
                zVar_data = permute(zVar_data, [2,3,1]);
            elseif ndims(zVar_data) > 3
                error('write_CDF_dataobj:Assertion:IllegalArgument:OperationNotImplemented', ...
                    'Found zVar data with more than three dimensions. Functionality has not been tested for this.')
            end
        end
        

        
        %===========================================================================================================
        % Convert specific VariableAttributes values.
        % Case 1: tt2000 values as UTC strings : Convert to tt2000.
        % Case 2: All other                    : Convert to the zVariable data type.
        % --------------------------------------------------------------------------
        % IMPLEMENTATION NOTE: spdfcdfread (not spdfcdfwrite) can crash if not doing this!!!
        % The tt2000 CDF variables are likely the problem(?).
        %
        % BUG: Does not seem to work on SCALEMIN/-MAX specifically despite identical treatment, for unknown reason.
        %===========================================================================================================
        %zVar_name
        for i_van = 1:length(VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE)   % van = variable attribute name
            var_attr_name = VARIABLE_ATTRIBUTES_OF_VARIABLE_TYPE{i_van};
            if ~isfield(dataobj_VariableAttributes, var_attr_name)
                continue
            end
            
            % IMPLEMENTATION NOTE: Can NOT assume that every CDF variable is represented among the cell arrays in
            % dataobj_VariableAttributes.(...).
            % Example: EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-CWF_e1d0a9a_CNE_V01.cdf
            VA_field = dataobj_VariableAttributes.(var_attr_name);
            i_vav = find(strcmp(VA_field(:,1), zVar_name));
            i_van = i_vav;
            if length(i_vav) == 0
                % CASE: The current zVariable does not have this attribute (var_attr_name).
                continue
            elseif length(i_vav) > 1
                error('write_CDF_dataobj:Assertion:OperationNotImplemented', ...
                    'Can not handle multiple variable name matches in dataobj_VariableAttributes.%s.', var_attr_name)
            end
            var_attr_value = VA_field{i_vav, 2};
            if strcmp(dataobj_stated_data_type, 'tt2000') && ischar(var_attr_value)
                %var_attr_name
                %var_attr_value
                var_attr_value = spdfparsett2000(var_attr_value);   % Convert char-->tt2000.
            elseif ~strcmp(dataobj_stated_data_type, class(var_attr_value))
                %zVar_name
                %dataobj_stated_data_type
                %var_attr_name
                %var_attr_value
                class(var_attr_value)
                error('write_CDF_dataobj:Assertion', ...
                    'Found VariableAttribute %s for CDF variable %s whose data type did not match the declared one.', ...
                    var_attr_name, zVar_name)
            end
            VA_field{i_van, 2} = var_attr_value;
            dataobj_VariableAttributes.(var_attr_name) = VA_field;
        end

        VARIABLE_LIST(end+[1,2]) = {zVar_name, zVar_data               };        
        VARDATATYPE  (end+[1,2]) = {zVar_name, dataobj_stated_data_type};
        PADVALS      (end+[1,2]) = {zVar_name, pad_value               };
    end
    
    % dataobj_VariableAttributes = rmfield(dataobj_VariableAttributes, {'VALIDMIN', 'VALIDMAX', 'SCALEMIN', 'SCALEMAX', 'FILLVAL'});
    
    
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
irf.log('n', sprintf('Writing CDF file "%s".', file_path))
spdfcdfwrite(file_path, VARIABLE_LIST(:), 'RecordBound', RECBNDVARS, 'GlobalAttributes', dataobj_GlobalAttributes, ...
    'VariableAttributes', dataobj_VariableAttributes, 'Vardatatypes', VARDATATYPE, 'PadValues', PADVALS)

end
