%
% Function which writes a CDF file.
%
% Attempt at a function which can easily write a CDF using variables on the same
% data format as returned by dataobj (irfu-matlab). Useful for reading a CDF
% file, modifying the contents somewhat, and then writing the modified contents
% to a CDF file. Originally based on write_cdf.m/write_cdf_spdfcdfread.m.
% Primarily used by BICAS (SolO; BIAS CAlibration Software).
%
% NOTE: This function is not fully generic, but quite (i.e. it can not handle
% all cases). See comments.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-07-12 (as write_cdf.m/write_cdf_spdfcdfread.m), 2016-10-20
% (as write_cdf_dataobj.m)
%
%
%
% ARGUMENTS
% =========
% filePath
%       Path to file to create.
% dataobj_GlobalAttributes
%       The corresponding field of an instantiated dataobj. Struct where
%       .<global attribute name>{i} = global attribute value as string (must be
%       strings?).
%       NOTE: Unclear if cell array must be 1D.
% dataobj_data
%       The corresponding field of an instantiated dataobj. Struct where
%       .<ZV name>.data = Array (numeric or char) with zVar content.
%           Indices for numeric arrays: (iRecord, i1, i2, ...)
%           Indices for char    arrays: Same as in dataobj.data.<zvName>.data,
%                                       i.e. inconsistent.
%       .<ZV name>.dim  : NOT USED. dataobj: Size of record. Row vector,
%                         at least length 2.
%       NOTE: May have other fields than ".data" which are then ignored.
% dataobj_VariableAttributes
%       The corresponding field of an instantiated dataobj. Struct where
%       .(<ZVA name>){iZv, 1} = ZV(!) name (ZVA name specified as field name)
%       .(<ZVA name>){iZv, 2} = ZVA value
% dataobj_Variables
%       The corresponding field of an instantiated dataobj. Nx12 cell array where
%       {iZv,  1} = ZV name
%       {iZv,  2} = Size of record. Row vector, at least length 2.
%       {iZv,  3} = NOT USED. dataobj: Uknown meaning. Scalar number.
%       {iZv,  4} = String representing data type (tt2000, single, char etc)
%       {iZv,  5} = NOT USED. dataobj: "Record variance", string representing
%                     on which dimensions the ZV changes. T=True,
%                     F=False.
%       {iZv,  6} = NOT USED. dataobj: Uknown meaning. String = 'Full' (always?)
%       {iZv,  7} = Compression algorithm, if any.
%       {iZv,  8} = NOT USED. dataobj: Uknown meaning. Scalar number
%       {iZv,  9} = Pad value
%       {iZv, 10} = NOT USED. dataobj: Uknown meaning. Scalar number or empty.
%       {iZv, 11} = NOT USED. dataobj: Uknown meaning. Scalar number or empty.
%       {iZv, 12} = NOT USED. dataobj: Uknown meaning. Scalar number or empty.
%       "NOT USED" = Not used by this function.
% varargin
%       Settings passed to irf.utils.interpret_settings_args(). See
%       implementation.
%
%
%
% LIMITATIONS
% ===========
% NOTE/PROBLEM: spdfcdfread() and spdfcdfinfo() may crash MATLAB(!) when reading
% files written with spdfcdfwrite() (which this function uses). It appears that
% this happens when spdfcdfwrite() receives various forms of "incomplete" input
% data. spdfcdfwrite() appears to often not give any warning/error message when
% receiving such data and writes a file anyway with neither error nor warning.
% Before passing data to spdfcdfwrite(), this function tries to give errors for,
% or correct such data, but can only do so as far as the problem is understood
% by the author. Submitting empty data for a CDF variable is one such case.
% Therefore, despite best efforts, this function might still produce nonsensical
% files instead of producing any warning or error message.
%
% NOTE PROBLEM(?): Can not select CDF encoding (Network/XDR, IBMPC etc) when
% writing files. The NASA SPDF MATLAB CDF Software distribution does not have
% this option (this has been confirmed with their email support 2016-07-14).
%
% BUG: Variable attributes SCALEMIN, SCALEMAX for Epoch are stored as CDF_INT8
% (not CDF_TT2000) in the final CDF file. The information stored seems correct
% though. Therefore, the same variable attributes are also represented as
% integers when reading the CDF file with dataobj.
%
% BUG/NOTE: spdfcdfwrite() has been observed to set the wrong pad value when
% writing 0 records.
%
% NOTE: spdfcdfwrite() always writes char as UCHAR (not CHAR) in the CDF.
%
% NOTE: The exact stated zVar dimensionality per record may be slightly wrong
% with regard to size=1 dimension:
% --The exact zvar size may be wrong. CDF file zVars of size 1 (scalar) per
% record may, in the CDF file, have specified record size "1:[1]" for 0--1
% records, and "0:[]" for >=2 records as displayed by cdfdump.
% --zVar dimensionality per record may have all its trailing ones removed
% (except for the case above).
%
% NOTE: dataobj may permute zVar dimensions for unknown reason (irfu-matlab
% commit 1a9a7c32, branch SOdevel).
% Ex: BIAS RCT zVar TRANSFER_FUNCTION_COEFFS (3 dimensions per record).
%
%
% IMPLEMENTATION NOTES
% ====================
% -- To keep the function as generic as possible, it does not contain any log
%    messages.
% -- The function does not accept a whole dataobj object since:
% (1) instances of dataobj are likely NOT meant to be modified after creation.
% Empirically, it is possible to modify them in practice though. Therefore the
% function only accepts the parts of a dataobj that it really needs, which still
% means it accepts a lot of redundant information in the arguments.
% (2) One might want to use the function for writing a CDF file without basing
% it on a dataobj (i.e. without basing it on an existing CDF file).
%
%
% IMPLEMENTATION NOTE: "spdfcdfwrite()" AND CHAR ZVARIABLES
% =========================================================
% The behaviour of spdfcdfwrite() when passing char arrays or cell arrays of
% strings for zVariables is very mysterious and hard to understand. Below is the
% empirical behaviour from passing such arrays to spdfcdfwrite() (RecordBound
% option disabled option, Singleton option enabled).
% --
% i = index within record. N,M,K>1
% Left  column = Size of CHAR ARRAY passed to scpdfcdfwrite.
% Right column = Result read from cdfdump (not dataobj).
% 0x0   : Error
% Mx1   : 1 record, M=i,     1 char/string!
% 1xN   : 1 record, N=strLen
% 1x1xK : 1 record, 1 char/string, all but first char lost!
% MxN   : 1 record, M=i, N=strLen
% MxNxK : 1 record, M=i, N=strLen, all but K index value=1 lost!
% NOTE: For the above, only 1 record is produced in all cases.
% NOTE: For the above, using RecordBound ONLY leads to more data being
% lost/ignored for some cases.
% --
% Left column  = Size of CELL ARRAY (of strings) passed to spdfcdfwrite().
% Right column = Result read from cdfdump (not dataobj).
% 0x0 : zVar is not written to file (still no error)!
% 1x1 : 1 record, 1 string/record
% 1x2 : 1 record, 2 strings/record
% 2x1 : 2 records, 1 string/record
% 3x2 : 3 records, 2 strings/record, BUT the strings are placed ILLOGICALLY/IN
% THE WRONG POSITIONS!
% 2x3 : 2 records, 3 strings/record. Strings are placed logically.
% NOTE: No alternative gives 1 record, with multiple strings.
% NOTE: For the above, using RecordBound makes no difference.
% --
% NOTE: dataobj always returns a char array but the meaning of indices seems to
% vary.
%
%
% DEFINITIONS
% ===========
% ZV  = CDF zVariable
% ZVA = CDF zVariable Attribute
% DO  = dataobj
%
function write_dataobj(filePath, ...
  dataobj_GlobalAttributes, ...
  dataobj_data, ...
  dataobj_VariableAttributes, ...
  dataobj_Variables, ...
  varargin)



%=======================================================================================================================
% PROPOSAL: Implement using NASA SPDFs Java code instead?!! Should be possible to easily call Java code from inside MATLAB.
%   PRO?: Java interface might be more easy to work with and have fewer quirks & limitations.
%
% PROPOSAL: Option for filling empty variable with pad values. 1 record?
%    CON: Can not know the (non-record) dimensions of such a CDF variable.
%    CON: Using exactly one record automatically leads to the CDF labelling the CDF variable as record-invariant!
%    CON: May not fit any zvar attribute DEPEND_x (zvar must have same length as other zvar).
%
% PROPOSAL: Reorganize into write_dataobj calling more genering function write_CDF which assumes more generic data
% structures.
%   PROPOSAL: Useful for combining with future generic function read_CDF which could replace dataobj.
%       NOTE/CON: spdfcdfread() returns some data structures similar to what dataobj contains so the gain might be small.
%
% PROPOSAL: Create analogous read_CDF+write_CDF (which use the same data structures). Combine with proper test code.
%   NOTE: This current code is based on writing a modified dataobj to disk, which is not necessarily desirable for a
%         general-purpose function write_CDF function.
%
%
% PROPOSAL: Some form of validation of input.
%    PROPOSAL: Assertions for redundant data (within dataobj data/attributes).
%       PROPOSAL: Check that both stated zvar sizes (within record) are consistent.
%       PROPOSAL: Check that stated nbr of records fits data.
%
% PROPOSAL: Shorten for-loop over zvars, by outsourcing tasks to functions.
%
% PROPOSAL: Write zVars using cell arrays of records (matrices) (spdfcdfwrite() permits it; see "RecordBound").
%
% PROPOSAL: Flag for different interpretations of indices in char arrays (dataobj or logical).
% PROPOSAL: Flag for assertion on preventing NaN.
%
% PROPOSAL: Check for illegal characters (or at least, characters which can not be handled) in global attributes:
%           åäöÅÄÖ, quotes(?).
%   NOTE: According to old notes, åäöÅÄÖ will be invisible global attributes, but the corresponding number of characters
%   will be cut out from the end.
%
% ~BUG: Zero-record numeric & char zVars are converted to [] (numeric, always 0x0) by dataobj. This code does not take
% this into account by internally converting the zVar variable back to the right class and size.
%
% PROPOSAL: Write automated tests for suitable helper functions:
%   construct_spdfcdfwrite_arguments()
%   prepare_ZV_value()
%   ensure_ZV_ZVA_data_types_consistent()
%   --
%   NOTE: Requires making inner functions, external functions.



% ZVAs that should reasonably have the same data type as the zVariable itself.
ZVA_ZV_SAME_DATA_TYPE_ZVA_NAMES_CA = {...
  'VALIDMIN', 'VALIDMAX', ...
  'SCALEMIN', 'SCALEMAX', ...
  'FILLVAL'};

DEFAULT_SETTINGS = struct();
% Whether zVariable value size per record must fit the submitted metadata
% specified in dataobj_Variables{i, 2}.
DEFAULT_SETTINGS.strictNumericZvSizePerRecord      = true;
% Default 1/true since dataobj is not strict about SIZE  of empty zVars.
DEFAULT_SETTINGS.strictEmptyNumericZvSizePerRecord = true;
% Default true since dataobj is not strict about CLASS of empty zVars.
DEFAULT_SETTINGS.strictEmptyZvClass                = true;
% Whether zVar attr should have same class as zVar.
% Exception: When zVar is TT2000 and zVar attr is char.
% Deactivation is useful for less stringent CDFs.
DEFAULT_SETTINGS.strictZvAttrClass                 = 'ERROR';   % Legal values: "ERROR", "WARNING", "IGNORE"
DEFAULT_SETTINGS.calculateMd5Checksum              = true;
%
Settings = irf.utils.interpret_settings_args(DEFAULT_SETTINGS, varargin);
irf.assert.struct(Settings, fieldnames(DEFAULT_SETTINGS), {})
assert(islogical(Settings.strictNumericZvSizePerRecord))
assert(islogical(Settings.strictEmptyNumericZvSizePerRecord))
assert(islogical(Settings.strictEmptyZvClass))
assert(islogical(Settings.calculateMd5Checksum))



% ASSERTION: ZV names are all unique and unambiguous
% ---------------------------------------------------------
% zvNameAllCa1 previously called for only non-char data. Why?
zvNameAllCa1 = dataobj_Variables(:, 1);
zvNameAllCa2 = fieldnames(dataobj_data);
%
irf.assert.castring_set(zvNameAllCa1)
irf.assert.castring_sets_equal(zvNameAllCa1, zvNameAllCa2)

zvNameAllCa = zvNameAllCa1;
clear zvNameAllCa1 zvNameAllCa2



[A, dataobj_VariableAttributes] = construct_spdfcdfwrite_arguments(...
  dataobj_Variables, ...
  dataobj_VariableAttributes, ...
  dataobj_data, ...
  Settings, ...
  ZVA_ZV_SAME_DATA_TYPE_ZVA_NAMES_CA);



%===================================================================================================
% RELEVANT spdfcdfwrite() OPTIONS:
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
%   attribute struct.  The "CDFAttribute Rename" field musdataobjStatedMatlabClasst have a value
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
%
%   SPDFCDFWRITE(..., 'Singleton', VARS, ...) indicates whether to keep the
%   singleton dimension(s) passed in from the multi-dimensional data. VARS is
%   a cell array of variable names, indicating each variable's singleton
%   dimension(s) is to be kept.
%   For example, variable with data dimensions like 10x1x100 will be written
%   as 2-dimensions (10x1) for 100 records if the record bound is specified.
%   For a row (1-by-M) or column (M-by-1) vector, the variable data will be
%   written as 2-dimension as is, unless the recordbound is specified.
%   The default setting is to have all singleton dimension(s) removed.
%   The above 10x1x100 variable will be written as 1-dimension
%   (with 10 elements).
%===================================================================================================
if Settings.calculateMd5Checksum ; checksumFlagArg = 'MD5';
else                             ; checksumFlagArg = 'None';
end

spdfcdfwrite(...
  filePath,             A.zvNameAndValueCa(:), ...
  'RecordBound',        A.zvNameRecordBoundCa, ...
  'GlobalAttributes',   dataobj_GlobalAttributes, ...
  'VariableAttributes', dataobj_VariableAttributes, ...
  'Vardatatypes',       A.zvNameAndDataTypeCa, ...
  'PadValues',          A.zvNameAndPadValueCa, ...
  'VarCompress',        A.zvNameAndCompressionCa, ...
  'Singleton',          zvNameAllCa, ...
  'Checksum',           checksumFlagArg)

end    % function







% Construct variables which spdfcdfwrite() accept as arguments. spdfcdfwrite()
% uses peculiarly formatted arguments.
%
% RETURN VALUE
% ============
% A
%       Struct with fields for different spdfcdfwrite() arguments.
%
function [A, dataobj_VariableAttributes] = construct_spdfcdfwrite_arguments( ...
  dataobj_Variables, ...
  dataobj_VariableAttributes, ...
  dataobj_data, ...
  Settings, ...
  ZVA_ZV_SAME_DATA_TYPE_ZVA_NAMES_CA)

zvNameRecordBoundCa = {};  % Refers to spdfcdfwrite() option "RecordBound".

% Lists where pairs of successive components contain (a) zVar name, and (b)
% corresponding zVar value/data type/pad value. These lists are needed as
% arguments to spdfcdfwrite(), which requires that very format.
zvNameAndValueCa       = {};
zvNameAndDataTypeCa    = {};
zvNameAndPadValueCa    = {};
zvNameAndCompressionCa = {};

for iZv = 1:length(dataobj_Variables(:,1))

  %===========================================================================
  % Extract ZV data from arguments
  % -------------------------------------
  % IMPLEMENTATION NOTE: Not using (1) data(i).VariableName or (2)
  % info.Variables(:,1) to obtain the variable name since experience shows
  % that components of (1) can be empty (contain empty struct fields) and (2)
  % may not cover all variables when obtained via spdfcdfread()!!
  %===========================================================================
  zvName                 = dataobj_Variables{iZv, 1};
  specifiedSizePerRecord = dataobj_Variables{iZv, 2};
  % "CdfDataType" refers to that the value should be interpreted as a CDF
  % standard string for representing data type (not a MATLAB class/type):
  % uint32, tt2000 etc.
  specifiedCdfDataType   = dataobj_Variables{iZv, 4};

  % Whether ZV is compressed or not, and by what method and degree.
  zvCompression          = dataobj_Variables{iZv, 7};

  % This value can NOT be found in dataobj_data. Has to be read from
  % dataobj_Variables.
  padValue               = dataobj_Variables{iZv, 9};

  zvValue                = dataobj_data.(zvName).data;
  specifiedMatlabClass   = irf.cdf.convert_CDF_type_to_MATLAB_class(...
    specifiedCdfDataType, 'Permit MATLAB classes');



  % ASSERTION: No zero-size dimensions (in size per record)
  %
  % IMPLEMENTATION NOTE: Code can not handle zero size dimensions (in size
  % per record).
  % In practice: #records > 0 with zero-size records ==> zero records
  % Not certain that the CDF files format is meant to handle this either.
  if prod(specifiedSizePerRecord) == 0
    error('write_dataobj:Assertion', ...
      ['Specified size per record contains zero-size dimension(s).', ...
      ' This function can not handle this case.'])
  end

  %zvValue = handle_zero_records(zvValue, padValue, dataobjStatedMatlabClass, turnZeroRecordsIntoOneRecord);



  %=========================================================================
  % ASSERTION:
  %   Check that the supplied ZV data variable has a MATLAB class
  %   (type) which matches the specified CDF type.
  % -------------------------------------------------------------
  % IMPLEMENTATION NOTE:
  % (1) Empty data (empty arrays) from spdfcdfread() are known to have the
  %     wrong data type (char). Therefore, do this check after having dealt
  %     dealt with empty data.
  % (2) Must do this after converting time strings (char) data to
  %     uint64/tt2000.
  %=========================================================================
  zvDataMatlabClass = class(zvValue);

  if ~strcmp( specifiedMatlabClass, zvDataMatlabClass ) && ...
      (Settings.strictEmptyZvClass || ~isempty(zvValue))

    error('write_dataobj:Assertion', ...
      ['The MATLAB class ("%s") of the variable containing zVariable', ...
      ' ("%s") data does not match specified CDF data type "%s".'], ...
      zvDataMatlabClass, zvName, specifiedCdfDataType)
  end



  [zvValue, isRecordBound] = prepare_ZV_value(...
    zvValue, specifiedSizePerRecord, Settings, zvName);
  if isRecordBound
    zvNameRecordBoundCa{end+1} = zvName;
  end



  dataobj_VariableAttributes = ensure_ZV_ZVA_data_types_consistent(...
    zvName, ...
    specifiedCdfDataType, ...
    specifiedMatlabClass, ...
    dataobj_VariableAttributes, ...
    ZVA_ZV_SAME_DATA_TYPE_ZVA_NAMES_CA);


  zvNameAndValueCa      (end+[1,2]) = {zvName, zvValue             };
  zvNameAndDataTypeCa   (end+[1,2]) = {zvName, specifiedCdfDataType};
  zvNameAndPadValueCa   (end+[1,2]) = {zvName, padValue            };
  zvNameAndCompressionCa(end+[1,2]) = {zvName, zvCompression       };
end    % for

% Construct function return value.
A = [];
A.zvNameAndValueCa       = zvNameAndValueCa;
A.zvNameAndDataTypeCa    = zvNameAndDataTypeCa;
A.zvNameAndPadValueCa    = zvNameAndPadValueCa;
A.zvNameRecordBoundCa    = zvNameRecordBoundCa;
A.zvNameAndCompressionCa = zvNameAndCompressionCa;

end







% Ensure that specific ZVAs have a data type that is consistent with the
% corresponding ZV data type.
%
% For specified ZVAs (any ZV):
% ----------------------------
% Case 1: ZV is tt2000 and ZVA is UTC string
%   ==> Convert ZVA to tt2000.
% Case 2: Anything else
%   ==> Assert that ZV and ZVA data types are identical.
%
% BUG: Does not seem to work on SCALEMIN/-MAX specifically despite
% identical treatment, for unknown reason.
%
% IMPLEMENTATION NOTE: spdfcdfread() (not spdfcdfwrite()) can crash if not
% doing this!!! The tt2000 CDF variables are likely the problem(?).
%
function dataobj_VariableAttributes = ensure_ZV_ZVA_data_types_consistent(...
  zvName, specifiedCdfDataType, specifiedMatlabClass, ...
  dataobj_VariableAttributes, ZVA_ZV_SAME_DATA_TYPE_ZVA_NAMES_CA)

% Ensure
% TT2000 ZV is TT2000 ZVA.
% ZV and ZVA data types consistent
% selected ZVA
% ensure_ZVA_ZV_data_types_consistent

for iZvaOfZvDataType = 1:length(ZVA_ZV_SAME_DATA_TYPE_ZVA_NAMES_CA)
  zvaName = ZVA_ZV_SAME_DATA_TYPE_ZVA_NAMES_CA{iZvaOfZvDataType};
  if ~isfield(dataobj_VariableAttributes, zvaName)
    % CASE: The current ZV zvName (iteration) does not have ZVA zvaName.
    continue
  end

  % IMPLEMENTATION NOTE: Can NOT assume that every CDF variable is
  % represented among the cell arrays in
  % dataobj_VariableAttributes.(...).
  % Example: EM2_CAL_BIAS_SWEEP_LFR_CONF1_1M_2016-04-15_Run1__e1d0a9a__CNES/ROC-SGSE_L2R_RPW-LFR-SURV-CWF_e1d0a9a_CNE_V01.cdf

  % Retrieve this ZVA (e.g. VALIDMIN) but for all zVariables (where
  % present) from argument standard struct dataobj_VariableAttributes.
  % Nx2 array.
  %
  % DOVAF = dataobj VariableAttributes Field
  dovafCa = dataobj_VariableAttributes.(zvaName);

  iDovafRow = find(strcmp(dovafCa(:,1), zvName));
  if isempty(iDovafRow)
    % CASE: The current zVariable does not have this attribute (zvaName).
    continue
  elseif length(iDovafRow) > 1
    error('write_dataobj:Assertion:OperationNotImplemented', ...
      ['Can not handle multiple zVariable name matches in', ...
      ' dataobj_VariableAttributes.%s.'], zvaName)
  end
  % CASE: iDovafRow is scalar.

  % =========================
  % Read and modify ZVA value
  % =========================
  zvaValue = dovafCa{iDovafRow, 2};
  if strcmp(specifiedCdfDataType, 'tt2000') && ischar(zvaValue)
    zvaValue = spdfparsett2000(zvaValue);   % Convert char-->tt2000.

  elseif ~strcmp(specifiedMatlabClass, class(zvaValue))
    msg = sprintf(...
      ['Found VariableAttribute %s for CDF zVariable "%s"', ...
      ' whose data type did not match the declared one.', ...
      ' specifiedCdfDataType="%s", specifiedMatlabClass="%s",', ...
      ' class(zvaValue)="%s"'], ...
      zvaName, zvName, specifiedCdfDataType, ...
      specifiedMatlabClass, class(zvaValue));

    switch(Settings.strictZvAttrClass)
      case 'ERROR'
        error('write_dataobj:Assertion', msg)
      case 'WARNING'
        warning('write_dataobj:Assertion', msg)
      case 'IGNORE'
        % Do nothing.
      otherwise
        error('Illegal setting Settings.strictZvAttrClass="%s".', ...
          Settings.strictZvAttrClass)
    end
  end

  % Modify dataobj_VariableAttributes correspondingly.
  dovafCa{iDovafRow, 2}                = zvaValue;
  dataobj_VariableAttributes.(zvaName) = dovafCa;
end    % for

end







% Convert a char array that dataobj returns into a char array that
% prepare_char_ZV_data() interprets the same way.
%
% ARGUMENTS
% =========
% charArray
%       Char array with indices (iRecord,iCharWithinString,)
%
function charArray = convert_dataobj_charZVValue_2_consistent_charZVValue(...
  charArray, nWrd1)

% ASSERTION
assert(isscalar(nWrd1), ...
  'write_dataobj:Assertion', 'Argument nWrd1 is not a scalar.')

if nWrd1 == 1
  charArray = permute(charArray, [2,1,3]);
else
  charArray = permute(charArray, [2,3,1]);
end
end







% Function for converting a char array representing a char zVariable using a
% consistent and logical indexing scheme, into the VERY HARD-TO-UNDERSTAND
% scheme that spdfcdfwrite() requires to produce the desired zVariable.
%
% NOTE: If one wants another order of indices for charArray, then one should use
% permute() rather than change the algorithm.
%
%
% ARGUMENTS
% =========
% charArray
%       Array of chars with indices (iCharWithinString, iRecord, iWrd).
%       WRD = Within-Record Dimension
%       Must not have more dimensions than 3.
%       Must not have 0 elements.
%       Must not have both multiple records AND multiple strings per
%       record(!).
%
%
% RETURN VALUES
% =============
% zvValue
%       The variable that should be passed to spdfcdfwrite(). Can be (1) char
%       array, or (2) cell array of strings.
% isRecordBound
%       True/false. Whether the zVariable should be passed to spdfcdfwrite()
%       with option "RecordBound" enabled.
%
function [zvValue, isRecordBound] = prepare_char_ZV_data(charArray)

% ASSERTIONS. Important to check that the code can actually handle the case.
assert(ischar(charArray), ...
  'write_dataobj:Assertion', ...
  'Argument charArray is not a char array.')
assert(ndims(charArray) <= 3, ...
  'write_dataobj:Assertion:OperationNotImplemented', ...
  ['Argument charArray has more than 3 dimension (2 per record).', ...
  ' Can not produce value for such zVariable.'])
assert(~isempty(charArray), ...
  'write_dataobj:Assertion:OperationNotImplemented', ...
  ['Argument charArray constains zero strings.', ...
  ' Can not produce value for empty zVariable.'])

% WRD1 = Within-Record Dimension 1.
nRecords = size(charArray, 2);   % CASE: >=1, because of assertion.
nWrd1    = size(charArray, 3);   % CASE: >=1, because of assertion.

if nRecords == 1
  if nWrd1 == 1
    zvValue = permute(charArray, [2, 1, 3]);
  else
    zvValue = permute(charArray, [3, 1, 2]);
  end
else
  if nWrd1 == 1
    zvValue = cell(nRecords, nWrd1);
    for iRecord = 1:nRecords
      for iWrd1 = 1:nWrd1
        zvValue{iRecord, nWrd1} = ...
          permute(charArray(:, iRecord, iWrd1), [2, 1, 3]);
      end
    end
  else
    error('write_dataobj:Assertion:OperationNotImplemented', ...
      ['Argument charArray represents multiple records containing', ...
      ' multiple strings per record. Can not produce zVariable', ...
      ' value for this case.']);
  end
end

isRecordBound = 0;    % Always!
end    % function







% Modify ZV value so that it can be passed to spdfcdfwrite() and be interpreted
% correctly.
%
% ARGUMENTS
% =========
% zValue
%       If numeric array: indices=(iRecord, i1, i2, ...)
%       If char array   : indices are the same as in
%                         dataobj.data.<zvName>.data, i.e. inconsistent.
% specifiedSizePerRecord
%       Size per record used for assertion.
%       For numeric: zValue size minus the first value, "size(zvValue)(2:end)".
%
function [zvValue, isRecordBound] = prepare_ZV_value(...
  zvValue, specifiedSizePerRecord, Settings, zvName)

if ischar(zvValue)
  %==========================================================================
  % CASE: char zVar: Convert 3-D char matrices to column cell arrays of
  % 2-D char matrices.
  % ------------------------------------------------------------------------
  % IMPLEMENTATION NOTE: It is not possible to permute indices for string
  % as one can for non-char for ndim==3.
  %==========================================================================

  zvValue = convert_dataobj_charZVValue_2_consistent_charZVValue(...
    zvValue, specifiedSizePerRecord(1));

  %=======================================
  % ASSERTION: Check zVar size per record
  %=======================================
  % NOTE: This check can not be perfect since zvValue with multiple
  % strings can be interpreted correctly for two different values of
  % specifiedSizePerRecord: 1 (multiple strings in one record) and non-1
  % (multiple records, with one string per record).
  temp          = size(zvValue);
  % NOTE: Throw away indices iRecord and iCharWithinString.
  sizePerRecord = temp(3:end);
  if ~isequal(...
      normalize_size_vec(specifiedSizePerRecord), ...
      normalize_size_vec(sizePerRecord))
    error('write_dataobj:Assertion', ...
      ['The zVariable data size (dataobj_data.(''%s'').data) does', ...
      ' not fit the stated size per record (dataobj_Variables).'], ...
      zvName)
  end

  [zvValue, isRecordBound] = prepare_char_ZV_data(zvValue);



elseif isnumeric(zvValue)

  nRecords = size(zvValue, 1);
  if Settings.strictNumericZvSizePerRecord || ...
      (Settings.strictEmptyNumericZvSizePerRecord && (nRecords == 0))
    % NOTE: dataobj zVar data is always (empirically) [] (i.e. numeric
    % 0x0) when nRecords=0, i.e. also for char-valued zVars, and also
    % for non-empty size per record. Therefore code often needs to be
    % tolerant of this. Note that the code can not (?) reconstruct an
    % original char zVar from dataobj for nRecords=0 since it does not
    % have the length of the strings.

    %=======================================
    % ASSERTION: Check zVar size per record
    %=======================================
    temp          = size(zvValue);
    sizePerRecord = temp(2:end);
    if ~isequal(...
        normalize_size_vec(specifiedSizePerRecord), ...
        normalize_size_vec(sizePerRecord))

      sizePerRecordStr       = ['[', ...
        strjoin(irf.str.sprintf_many('%i', sizePerRecord), ', '), ...
        ']'];
      specifiedSizePerRecordStr = ['[', ...
        strjoin(irf.str.sprintf_many('%i', specifiedSizePerRecord), ', '), ...
        ']'];

      error('write_dataobj:Assertion', ...
        ['The zVariable "%s" data size according to data variable', ...
        ' itself is not consistent with the stated size per record', ...
        ' in other argument.\n', ...
        '    Size per record according to data variable produced', ...
        ' by processing: %s\n', ...
        '    Size per record separately specified:              ', ...
        '                %s'], ...
        zvName, sizePerRecordStr, specifiedSizePerRecordStr)
    end
  end



  %===================================================================================
  % Special behaviour for numeric matrices with >=2D per record
  % -----------------------------------------------------------
  % For 3D matrices, spdfcdfwrite() interprets the last index (not the first
  % index!) as the record number. Must therefore permute the indices so that
  % write_cdf2 is consistent for all numbers of dimensions.
  %     write_dataobj data arguments : index 1 = record.
  %     matrix passed on to spdfcdfwrite() : index 3 = record.
  % NOTE: spdfcdfread() (at least with argument "'Structure', 1,
  % 'KeepEpochAsIs', 1") works like spdfcdfwrite() in this regard.
  %
  % Excerpt from the comments in "spdfcdfwrite.m":
  % ----------------------------------------------
  %   """"SPDFCDFWRITE(..., 'RecordBound', RECBNDVARS) specifies data values in arrays
  %   (1-D or multi-dimensional) are to be written into "records" for the given
  %   variable. RECBNDVARS is a cell array of variable names. The M-by-N array
  %   data will create M rows (records), while each row having N elements. For
  %   examples, 5-by-1 array will create five (5) scalar records and 1-by-5 array
  %   will write out just one (1) record with 5 elements. For 3-D array of
  %   M-by-N-by-R, R records will be written, and each record with M-by-N
  %   elements. Without this option, array of M-by-N will be written into a single
  %   record of 2-dimensions. See sample codes for its usage.""""
  %
  %   """"SPDFCDFWRITE(..., 'Singleton', VARS, ...) indicates whether to keep the
  %   singleton dimension(s) passed in from the multi-dimensional data. VARS is
  %   a cell array of variable names, indicating each variable's singleton
  %   dimension(s) is to be kept.
  %   For example, variable with data dimensions like 10x1x100 will be written
  %   as 2-dimensions (10x1) for 100 records if the record bound is specified.
  %   For a row (1-by-M) or column (M-by-1) vector, the variable data will be
  %   written as 2-dimension as is, unless the recordbound is specified.
  %   The default setting is to have all singleton dimension(s) removed.
  %   The above 10x1x100 variable will be written as 1-dimension
  %   (with 10 elements).""""
  %==================================================================================
  %if nRecords == 0
  %    zvValue = zeros(sizePerRecord);
  %else
  if nRecords == 1
    % Shift/permute indices "left" so that index 1 appears last (and
    % hence "disappears" since it is size=1 due to how MATLAB handles
    % indices).
    zvValue = shiftdim(zvValue, 1);
    isRecordBound = 0;
  else
    % CASE: First index size>=2.
    if ndims(zvValue) >= 3
      % Shift/permute indices "left" so that index 1 appears last
      % where it will be interpreted as number of records.
      zvValue = shiftdim(zvValue, 1);
    end
    isRecordBound = 1;
  end
else
  error('write_dataobj:Assertion', 'zvValue is neither char nor numeric.')
end

end    % function







% "Normalize" a size vector, i.e. 1D vector describing size (dimensions) of
% variable. Forces a row vector. Removes trailing ones (explicit size one for
% higher dimensions). Using this makes size vectors easily (and safely)
% comparable.
%
% NOTE: size() returns all dimensions up until the last non-size-one dimension.
% This also means that all zero-sized dimensions are included.
%       size() adds trailing ones up to 2D;
%
% []    ==> []  (size 1x0)
% [0]   ==> [0]
% [1 0] ==> [1 0]
% [0 1] ==> [0]
% [1 1] ==> [] (1x0)
%
function sizeVec = normalize_size_vec(sizeVec)
% IMPLEMENTATION NOTE: sizeVec = [] ==> find returns [] (not a number) ==> 1:[],
% but that gives the same result as 1:0 so the code works anyway.
sizeVec = sizeVec(1:find(sizeVec ~= 1, 1, 'last'));    % "Normalize" size vector.
end







% Handle special case for zero-record zVariables: (1) error, or (2) modify
% zvValue
%
% function zvValue = handle_zero_records(zvValue, padValue, specifiedMatlabClass, turnZeroRecordsIntoOneRecord)
%
%     if isempty(zvValue)
%         if ~turnZeroRecordsIntoOneRecord
%             error('write_dataobj:Assertion', ...
%                 ['Can not handle CDF zVariables with zero records', ...
%                 ' (due to presumed bug in spdfcdfwrite()).'])
%         else
%             %--------------------------------------------------------------------
%             % EXPERIMENTAL SOLUTION: Store 1 record of data with only pad values
%             % instead of zero records.
%             % NOTE: Incomplete since does not take CDF variable type, array
%             % dimensions into account.
%             %--------------------------------------------------------------------
%             nRecords = 1;
%             try
%                 zvValue = cast(ones(nRecords, 1), specifiedMatlabClass) * padValue;
%             catch exception
%                 error('write_dataobj:Assertion', ...
%                     ['Can not type cast zvar data variable to', ...
%                     ' MATLAB class "%s" (CDF: "%s").'], ...
%                     specifiedMatlabClass, dataobjStatedCdfDataType)
%             end
%         end
%     end
%
% end
