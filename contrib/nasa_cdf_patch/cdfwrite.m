function cdfwrite(filename, varcell, varargin)
%CDFWRITE Write data to a CDF file.
% 
%   CDFWRITE(FILE, VARIABLELIST, ...) writes out a CDF file whose name
%   is specified by FILE.  VARIABLELIST is a cell array of ordered
%   pairs, which are comprised of a CDF variable name (a string) and
%   the corresponding CDF variable value.  To write out multiple records
%   for a variable, put the variable values in a cell array, where each
%   element in the cell array represents a record.
%
%   CDFWRITE(..., 'PadValues', PADVALS) writes out pad values for given
%   variable names.  PADVALS is a cell array of ordered pairs, which
%   are comprised of a variable name (a string) and a corresponding 
%   pad value.  Pad values are the default value associated with the
%   variable when an out-of-bounds record is accessed.  Variable names
%   that appear in PADVALS must appear in VARIABLELIST.
%
%   CDFWRITE(..., 'GlobalAttributes', GATTRIB, ...) writes the structure
%   GATTRIB as global meta-data for the CDF.  Each field of the
%   struct is the name of a global attribute.  The value of each
%   field contains the value of the attribute.  To write out
%   multiple values for an attribute, the field value should be a
%   cell array.
%
%   In order to specify a global attribute name that is illegal in
%   MATLAB, create a field called "CDFAttributeRename" in the 
%   attribute struct.  The "CDFAttribute Rename" field must have a value
%   which is a cell array of ordered pairs.  The ordered pair consists
%   of the name of the original attribute, as listed in the 
%   GlobalAttributes struct and the corresponding name of the attribute
%   to be written to the CDF.
%
%   CDFWRITE(..., 'VariableAttributes', VATTRIB, ...) writes the
%   structure VATTRIB as variable meta-data for the CDF.  Each
%   field of the struct is the name of a variable attribute.  The
%   value of each field should be an mx2 cell array where m is the
%   number of variables with attributes.  The first element in the
%   cell array should be the name of the variable and the second
%   element should be the value of the attribute for that variable.
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
%   CDFWRITE(..., 'WriteMode', MODE, ...) where MODE is either 'overwrite'
%   or 'append' indicates whether or not the specified variables or attributes
%   should be appended to the CDF if the file already exists.  The 
%   default is 'overwrite', indicating that CDFWRITE will not append
%   variables and attributes.
%
%   CDFWRITE(..., 'Format', FORMAT, ...) where FORMAT is either 'multifile'
%   or 'singlefile' indicates whether or not the data is written out
%   as a multi-file CDF.  In a multi-file CDF, each variable is stored
%   in a *.vN file where N is the number of the variable that is
%   written out to the CDF.  The default is 'singlefile', which indicates
%   that CDFWRITE will write out a single file CDF.  When the 'WriteMode'
%   is set to 'Append', the 'Format' option is ignored, and the format
%   of the pre-existing CDF is used.
%
%   CDFWRITE(..., 'Version', VERSION, ...) where VERSION is a string which 
%   specifies the version of the CDF library to use in writing the file.
%   The default option is to use the latest version of the library 
%   (currently version 3.2+), and has to be specified '3.0'.  The 
%   other available version is version 2.7 ('2.7').  Note that 
%   versions of MATLAB before R2006b will not be able to read files 
%   which were written with CDF versions greater than 3.0.
%
%   CDFWRITE(..., 'RecordBound', RECBNDVARS) indicates whether
%   data values are written into "records" for the given variable
%   names. RECBNDVARS is a cell array of variable names.
%   This option applies to data values in 1-by-N or N-by-1 vector that
%   normally will be written into a simple record of N elements. With
%   this option, N records are written, each record being a scalar, instead. 
%   Previously, cell arrays are the only way to create multiple records. 
%   These two commands should write out the data values identically:
%      CDFWRITE('example', {'Epoch', num2cell(1:100)}, ...
%               'epochiscdfepoch', true);
%      CDFWRITE('example', {'Epoch', 1:100}, 'recordbound', {'Epoch'}, ...
%               'epochiscdfepoch', true);
%
%   CDFWRITE(..., 'ConvertDatenumToEpoch', TF, ...) converts MATLAB datenum
%   values to CDF epoch data if TF is true. This option works with the
%   designated variable 'Epoch' (case sensitive) in a CDF, which is epoch. 
%   There are two ways to write data for 'Epoch' variable of CDF_EPOCH into
%   a CDF. First, uses cdfepoch objects, each of which is from a datenum,
%   datestr, or even cdfepoch object, to pass data to cdfwrite: 
%      CDFWRITE(..., {'Epoch',cdfepoch([...], ...)}, ...) 
%   This option is time and space consuming if large datasets are involved. 
%   The second way uses 'ConvertDatenumToEpoch':
%      CDFWRITE(..., {'Epoch',[...], ...}, 'ConvertDatenumToEpoch', TF, ...) 
%   If TF is set to true, the passed data are assumed as MATLAB datenum
%   values and a data conversion to CDF epoch values is performed. 
%   Setting it to false (the default), 'Epoch' is treated just like a normal
%   variable, not an epoch type any more. 
%
%   CDFWRITE(..., 'ConvertDatenumToTT2000', TF, ...) converts MATLAB datenum
%   values to CDF TT2000 data if TF is true. This option works with the same
%   variable 'Epoch' as CDF Epoch for epoch values.
%   There are two ways to write data for 'Epoch' variable of CDF_TIME_TT2000
%   into a CDF. First, uses cdftt2000 objects, each of which is from a datenum,
%   datestr, or even cdftt2000 object, to pass data to cdfwrite:
%      CDFWRITE(..., {'Epoch',cdftt2000([...], ...)}, ...)
%   This option is also time and space consuming if large datasets are involved.
%   The second way uses 'ConvertDatenumToTT2000':
%      CDFWRITE(..., {'Epoch',[...], ...}, 'ConvertDatenumToTT2000', TF, ...)
%   If TF is set to true, the passed data are assumed as MATLAB datenum
%   values and a data conversion to CDF TT2000 values is performed.
%   Setting it to false (the default), 'Epoch' is treated just like a normal
%   variable, not an epoch type any more. This option is mutually exclusive
%   to ConvertDatenumToEpoch.
%
%   CDFWRITE(..., 'EpochIsCDFEpoch', TF, ...) controls whether the passed
%   data are actually CDF epoch values. This option applies only to the
%   designated variable 'Epoch'. If TF is set to true, the passed data are
%   actually CDF_EPOCH values so a CDF_EPOCH is to be created and all values
%   will be used as is - no MATLAB datenum to CDF epoch conversion is needed.
%   Without the option or the option is false, the 'Epoch' will be treated as
%   a regular variable, and written as of CDF_DOUBLE, data type, instead of 
%   CDF_EPOCH. This option is mutually exclusive to ConvertDatenumToEpoch 
%   and TT2000.
%
%   CDFWRITE(..., 'TT2000', TF, ...) controls whether the passed data are
%   actually CDF TT2000 values. This option applies only to the designated
%   variable 'Epoch'. If TF is set to true, the passed data are actually CDF
%   TT2000 values so they will be used as is- no conversion is needed. 
%   The variable 'Epoch' will be written to the file as type CDF_TIME_TT2000.
%
%
%
%   The CDF TT2000 values needs to be numeric of mxINT64_CLASS (int64).
%   This option is mutually exclusive to ConvertDatenumToEpoch and
%   EpochIsCDFEpoch.
%
%   Notes:
%
%     CDFWRITE creates temporary files when writing CDF files.  Both the
%     target directory for the file and the current working directory
%     must be writable.
%
%     To maximize performance, specify the 'ConvertDatenumToEpoch' 
%     parameter with true (nonzero) value while providing datenum values
%     for 'Epoch' variable. 
%
%     CDF epoch is the number of milliseconds since 1-Jan-0000 0h:0m:0s:000ms.
%     CDF TT2000 is the number of nanoseconds since 1-Jan-2000
%                12h:0m:0s:000000000ns with leap seconds included.
%
%
%   Examples:
%
%   % Write out a file 'example.cdf' containing a variable 'Longitude'
%   % with a single record (row) of a vector with 361 elements:
%
%   cdfwrite('example', {'Longitude', 0:360});
%
%   % Write out a file 'example.cdf' containing variables 'Longitude'
%   % and 'Latitude' with the variable 'Latitude' having a pad value
%   % of 10 for all out-of-bounds records that are accessed. The
%   % 'Longitude' variable has one record (row) of a vector with 361 elements
%   % 'Latitude' has one record of a vector with 11 elements.
%
%   cdfwrite('example', {'Longitude', 0:360, 'Latitude', 10:20}, ...
%            'PadValues', {'Latitude', 10});
%
%   % Write out a file 'example.cdf', containing a variable 'Longitude'
%   % with the value [0:360], and with a variable attribute of
%   % 'validmin' with the value 10:
%
%   varAttribStruct.validmin = {'Longitude' [10]};
%   cdfwrite('example', {'Longitude' 0:360}, ...
%            'VariableAttributes', varAttribStruct);
%
%   % Write out a file 'example.cdf' containing variables 'Longitude'
%   % and 'Latitude' with the variable 'Latitude' being a Record-bound:
%
%   cdfwrite('example', {'Longitude', 0:360, 'Latitude', 10:20}, ...
%            'RecordBound', {'Latitude'});
%
%   % Write out a file 'example.cdf', with multiple rows of time series data.
%   % Each row has a time and a sampled data, both being scalar. 
%   % The following sample shows 100 reords (rows): the time being of
%   % CDF_EPOCH data type and sampled value as double. 
%   % The epoch starts from 2010-01-01T00:00:00.000 with 1 second stride.
%   % Each sampled value starts from 0.0, with 5.0 stride. 
%
%   epoch=utc2cdfepoch(2010,1,1,0,0,0,0);
%   epochvalues = epoch+[0:99]*1000;
%   values = [0:99]*5.0;
%   cdfwrite('example', {'Epoch', epochvalues, 'Samples', values}, ...
%            'EpochisCDFEpoch', true, ...
%            'RecordBound', {'Epoch', 'Samples'});
%
%   % 'EpochisCDFEpoch' option is needed as 'Epoch' is to be created of
%   % CDF_EPOCH type.
%
%   % Alternatively, the same result can be accomplished by making the epoch
%   % vector to cell. No 'RECORDBOUND' option is needed.
%   epoch=utc2cdfepoch(2010,1,1,0,0,0,0);
%   epochvalues = epoch+[0:99]*1000;
%   epochscell = num2cell(epochvalues);
%   values = num2cell([0:99]*5.0);
%   cdfwrite('example', {'Epoch', epochvalues, 'Samples', values}, ...
%            'EpochisCDFEpoch', true);
%
%   % Write out a file 'example.cdf', with multiple rows of time sequenced,
%   % vectorized data. 
%
%   epoch=utc2cdfepoch(2010,1,1,0,0,0,0);
%   epochvalues = epoch+[0:2]*1000;
%   value = [0:3]*5.0;
%   values = [value; value+10; value+20];
%   cdfwrite('example', {'Epoch', epochvalues, 'Samples', values}, ...
%            'EpochisCDFEpoch', true, ...
%            'RecordBound', {'Epoch', 'Samples'});
%
%   % For writing out multi-dimensional epoch variable records, the following
%   % sample will write two records, each with a 2-dimension (2x3) array:
%
%   epoch=utc2cdfepoch(2010,1,1,0,0,0,0);
%   epochvalue1 = [[epoch epoch+1];[epoch+1000 epoch+1001];...
%                  [epoch+2000 epoch+2001]];
%   epochvalue2 = [[epoch+5000 epoch+5001];[epoch+6000 epoch+6001];...
%                  [epoch+7000 epoch+7001]];
%   epochs={epochvalue1;epochvalue2};
%   cdfwrite('example', {'Epoch', epochvalues}, ...
%            'EpochisCDFEpoch', true);
%
%   % Write out a file 'example.cdf', with 100 MATLAB datenum values,
%   % starting from 2010-01-01T00:00:00.000 with 1 second stride, 
%   % into variable: 'Epoch' of CDF_EPOCH type. The first record has a date:
%   %  01-Jan-2010 00:00:00.000, the second record:
%   %  01-Jan-2010 00:00:01.000, the third record:
%   %  01-Jan-2010 00:00:02.000, etc.
%   % 'Epoch' has a pad value of 01-Jan-0000 00:00:00.001.
%
%   datenum1=datenum(2010,1,1,0,0,0);
%   datenumvalues = datanum1+[0:99]/86400;
%   cdfwrite('example', {'Epoch', datenumvalues}, ...
%            'ConvertDatenumToEpoch', true, ...
%            'RecordBound', {'Epoch'}, ...
%            'PadValue', {'Epoch', 1});
%
%   % Write out a file 'example.cdf', with three records for the variable
%   % 'Epoch' of CDF_TIME_TT2000 type. The first record has a date:
%   % 2010-10-10T01:02:03.456, the second record: 2010-11-11T02:04:06.789
%   % and the third record: 2010-12-12T03:04:05.123.
%
%   dates = datenum([2010 10 10 1 2 3.456; 2010 11 11 2 4 6.789; ...
%                    2010 12 12 3 4 5.123]);
%   cdfwrite('example', {'Epoch', dates}, ...
%            'ConvertDatenumToTT2000', true, ...
%            'RecordBound', {'Epoch'});
%
%   % Write out a file 'example.cdf', with 6 records for the variable 
%   % 'Epoch', which will be of CDF_TIME_TT2000 data type. The data
%   $ crosses over a leap second. Convert the epoch in UTC string
%   $ to their TT2000 values before write out to the CDF file.
%
%   time = {'2008-12-31T23:59:58.123456789';
%           '2008-12-31T23:59:59.123456789';
%           '2008-12-31T23:59:60.123456789';
%           '2009-01-01T00:00:00.123456789';
%           '2009-01-01T00:00:01.123456789';
%           '2009-01-01T00:00:02.123456789'};
%   values = parsett2000(time);          
%   cdfwrite('example', {'Epoch', values}, ...
%            'TT2000', true, ...
%            'RecordBound', {'Epoch'});
%
%   See also CDFREAD, CDFUPDATE, CDFINFO, CDFEPOCH, CDFTT2000, ENCODETT2000,
%            COMPUTETT2000, PARSETT2000, BREAKDOWNTT2000, CDFLEAPSECONDSINFO.

%   Copyright 1984-2006 The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2013/04/15 15:21:41 $

% HISTORY:
%   February 12, 2009  Mike Liu     The following changes have been made to
%                                   cdfwritec.c:
%                                     - Added parameter 'RecordBound'.
%                                     - Added parameter 'ConvertDatenumToEpoch'.
%                                     - Added a logic to check CDF_EPOCH and 
%                                       CDF_EPOCH16 data.

%
% Process arguments.
%

if (nargin < 2)
    error('MATLAB:cdfwrite:inputArgumentCount', ...
          'CDFWRITE requires at least two input arguments.')
end

% parse_inputs sorts out all of the input args.  Its return values:
%
% * args - an array of structs.  args.VarNames contains the names
%          of the variables to be written to the CDF.  args.VarVals contains
%          the corresponding values.  args.PadVals contains corresponding pad
%          values. args.ConvertDatenum indicates whether to convert
%          MATLAB datenum to CDF epoch values. args.RecBnd contains variables
%          that are for Record-bound. args.EpochIsCDFEpoch is a flag
%          indicating whether datenum to cdf epoch conversion is needed.
% * isAppending - whether or not to delete this file or if we need to
%                 append to the file
% * isMultifile - whether or not to write out as a multi-file CDF
% * CDFversion - which version of the CDF library to use
% * varAttribStruct - a struct containing the variable attributes
% * globalAttribStruct - a struct containing the global CDF attributes
% * msg - an error message from parse_inputs that we pass on to the user.

[args, isAppending, isMultifile, CDFversion, varAttribStruct, ...
 globalAttribStruct, exception] = parse_inputs(varcell, varargin{:});

if (~isempty(exception.msg))
    error(exception.id, '%s', exception.msg)
end

%
% Create a proper filename for the CDF
%

% See if there is an extension
[pname, fname, ext] = fileparts(filename);

% If there is an extension, then remove it before passing it to CDFlib.
if (~isempty(ext))
    if (~isempty(strfind(ext, 'cdf')))
        filename = fullfile(pname, fname);
    end
end

%
% Call the underlying cdfwritec function which calls the CDFlib
%

cdfwritec(filename, args.VarNames, args.VarVals, args.PadVals, ...
          globalAttribStruct, varAttribStruct, isAppending, ...
          isMultifile, CDFversion, args.ConvertDatenum, args.RecBnd, ...
          args.EpochIsCDFEpoch, args.TT2000, args.ConvertDatenum2);

%%%
%%% Function parse_inputs
%%%

function [args, isAppending, isMultifile, CDFversion, varAttribStruct, ...
          globalAttribStruct, exception] = parse_inputs(varcell, varargin)

% Set default values
args.PadVals = {};
args.RecBnd = {};
args.EpochIsCDFEpoch = false;
args.TT2000 = false;
args.ConvertDatenum = false;
args.ConvertDatenum2 = false;
isAppending = 0;
isMultifile = 0;
varAttribStruct = struct([]);
globalAttribStruct = struct([]);
% The following value indicates no version preference.
CDFversion = -1.0;

exception.msg = '';
exception.id = '';

% First check that varcell meets all of our requirements
args.VarNames = {varcell{1:2:end}};
args.VarVals = {varcell{2:2:end}};
% Wrap the scalars non-empties in cell arrays.
for i = 1:length(args.VarVals)
    if ~isempty(args.VarVals{i}) && (ischar(args.VarVals{i}) || (numel(args.VarVals{i}) == 1))
        args.VarVals{i} = {args.VarVals{i}};    
    end
end

if length(args.VarNames) ~= length(args.VarVals)
    exception.msg = 'All variable names must have a corresponding variable value.';
    exception.id = 'MATLAB:cdfwrite:variableWithoutValue';
    return
end

% Check and make sure that all variable values are of the same
% datatype, but ignore empties
if ~isempty(args.VarVals)
    for i = 1:length(args.VarVals)
        a = args.VarVals{i};
        if iscell(a)
            nonEmpties = {a{~cellfun('isempty',a)}};
            if iscell(nonEmpties) && ~isempty(nonEmpties)
                dtype = class(nonEmpties{1});
                if ~all(cellfun('isclass',nonEmpties,dtype))
                    exception.msg = 'All record values for a given variable must be of the same type.';    
                    exception.id = 'MATLAB:cdfwrite:inconsistentRecordTypes';
                end
            end
        else
            % If it isn't a cell array, then it is an array and
            % all elements are of the same type.  This is a single
            % record value and must be placed in a cell array.
            args.VarVals{i} = {args.VarVals{i}};
        end
    end
end

args.PadVals = cell(1,length(args.VarNames));
args.RecBnd = cell(1,length(args.VarNames));

% Parse arguments based on their number.
if (nargin > 0)
    
    paramStrings = {'padvalues'
                    'globalattributes'
                    'variableattributes'
                    'writemode'
                    'format'
                    'recordbound'
                    'convertdatenumtoepoch'
                    'convertdatenumtott2000'
                    'epochiscdfepoch'
                    'version'
                    'tt2000'};
    
    % For each pair
    for k = 1:2:length(varargin)
        param = lower(varargin{k});
            
        if (~ischar(param))
            exception.msg = 'Parameter name must be a string.';
            exception.id = 'MATLAB:cdfwrite:paramNameNotString';
            return
        end
        
        idx = strmatch(param, paramStrings);
        
        if (isempty(idx))
            exception.msg = sprintf('Unrecognized parameter name "%s".', param);
            exception.id = 'MATLAB:cdfwrite:unrecognizedParam';
            return
        elseif (length(idx) > 1)
            exception.msg = sprintf('Ambiguous parameter name "%s".', param);
             exception.id = 'MATLAB:cdfwrite:ambiguousParam';
           return
        end
        
        switch (paramStrings{idx})
        case 'padvalues'
           padCell = varargin{k+1};
           % If we weren't passed an even pair, then a variable
           % name or value was left out.
           if rem(length(padCell), 2)
               exception.msg = ['Number of variables to write out with ' ...
                      'padding does not match number of pad values.'];
               exception.id = 'MATLAB:cdfwrite:paddingMismatch';
               return;
           end
           vars = {padCell{1:2:end}};
           padVals = {padCell{2:2:end}};
           % Check that vars are in the list above.
           if ~iscellstr(vars)
               exception.msg = 'All variable names must be strings.';
               exception.id = 'MATLAB:cdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(vars, args.VarNames))
               exception.msg = ['Variables listed in the PadValues ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:cdfwrite:notSavingVarForPadValue';
               return
           end
           for i = 1:length(padVals)
               padVal = padVals{i};
               if isnumeric(padVal) || ischar(padVal) || isa(padVal,'cdfepoch') || ...
                  isa(padVal,'cdfepoch16')
                   args.PadVals{strcmp(args.VarNames,vars{i})} = padVals{i};
               else
                   exception.msg = 'Pad values must be numbers, strings, cdfepochs or cdfepoch16s.';
                   exception.id = 'MATLAB:cdfwrite:badPadValue';
                   return
               end
           end
       case 'globalattributes'
           globalAttribStruct = varargin{k+1};
           if ~isstruct(globalAttribStruct)
               exception.msg = ['''GlobalAttributes''' ' must be a struct.'];
               exception.id = 'MATLAB:cdfwrite:globalAttributesNotStruct';
               return
           end
           attribs = fieldnames(globalAttribStruct);
           
           % If the global attribute isn't a cell, then stuff it in one.
           for i = 1:length(attribs)
               attribVal = globalAttribStruct.(attribs{i});
               if ~iscell(attribVal)
                   globalAttribStruct.(attribs{i}) = {attribVal};
               end
           end
        case 'variableattributes'
           varAttribStruct = varargin{k+1};
           if ~isstruct(varAttribStruct)
               exception.msg = ['''VariableAttributes''' ' must be a struct.'];
               exception.id = 'MATLAB:cdfwrite:variableAttributesNotStruct';
               return
           end
           attribs = fieldnames(varAttribStruct);
           
           % Check the VariableAttributes struct.
           for i = 1:length(attribs)
               % If the variable attribute isn't in a cell (because
               % it is scalar, then put it into a cell.
               attribVal = varAttribStruct.(attribs{i});
               s = size(attribVal);
               if ~iscell(attribVal)
                   varAttribStruct.(attribVal) = {attribVal};
               end
               % The variable attribute struct may have more than one
               % variable per attribute.  However, there must only be
               % one associated value of the attribute for each variable,
               % hence the 2.
               if (s(2) == 2)
                   % Transpose it because CDFlib reads the arrays column-wise.
                   varAttribStruct.(attribs{i}) = attribVal';   
               else
                   % We have ordered pairs.
                   varAttribStruct.(attribs{i}) = reshape(varAttribStruct.(attribs{i})(:),numel(varAttribStruct.(attribs{i})(:))/2, 2);
               end
               
%                % Don't forget to ignore the "CDFAttributeRename" attribute
%                completeSet = {args.VarNames{:} 'CDFAttributeRename'};
%                tmpVar = varAttribStruct.(attribs{i});
%                varsWithAttributes = {tmpVar{1,:}};
%                if ~all(ismember(varsWithAttributes, completeSet))
%                    exception.msg = ['Variables listed in the VariableAttributes ' ...
%                           'struct must be on the list of variables ' ...
%                           'to save.'];
%                    return
%                end               
           end
        case 'writemode'
            isAppending = varargin{k+1};
            if strcmpi(isAppending, 'overwrite')
                isAppending = 0;
            elseif strcmpi(isAppending, 'append')
                isAppending = 1;
            else
                exception.msg = ['''WriteMode''' ' must be either ' '''overwrite''' ... 
                       ' or ' '''append'''];
                exception.id = 'MATLAB:cdfwrite:badWriteModeValue';
                return
            end
        case 'format'
            isMultifile = varargin{k+1};
            if strcmpi(isMultifile, 'singlefile')
                isMultifile = 0;
            elseif strcmpi(isMultifile, 'multifile')
                isMultifile = 1;
            else
                exception.msg = ['''Format''' ' must be either ' '''singlefile''' ... 
                       ' or ' '''multifile'''];
                exception.id = 'MATLAB:cdfwrite:badFormatValue';
                return
            end
        case 'version'
            version = varargin{k+1};
            if ischar(version) && ...
                    (strcmp(version,'2.7') || strcmp(version,'3.0'))
                CDFversion = str2num(version);
            else
                exception.msg = '''Version'' must be either ''2.7'' or ''3.0''';
                exception.id = 'MATLAB:cdfwrite:badVersionValue';
                return
            end
        case 'recordbound'
           RecBndCell = varargin{k+1};
           % Check that vars are in the list above.
           if ~iscellstr(RecBndCell)
               exception.msg = 'All variable names must be strings.';
               exception.id = 'MATLAB:cdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(RecBndCell, args.VarNames))
               exception.msg = ['Variables listed in the RECORDBOUND ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:cdfwrite:notSavingVarForRecBnd';
               return
           end
           for i = 1:length(RecBndCell)
               RecBndVar = RecBndCell{i};
               args.RecBnd{strcmp(args.VarNames,RecBndCell{i})} = RecBndVar;
           end
        case 'convertdatenumtoepoch'
           if (k == length(varargin))
               msg = 'No datenum conversion value specified.';
               return
           else
               convert = varargin{k + 1};
               if (numel(convert) ~= 1)
                   msg = 'Datenum conversion value must be a scalar logical.';
               end

               if (islogical(convert))
                   args.ConvertDatenum = convert;
               elseif (isnumeric(convert))
                   args.ConvertDatenum = logical(convert);
               else
                   msg = 'Datenum conversion value must be a scalar logical.';
               end
           end
        case 'convertdatenumtott2000'
           if (k == length(varargin))
               msg = 'No datenum conversion value specified.';
               return
           else
               convert = varargin{k + 1};
               if (numel(convert) ~= 1)
                   msg = 'Datenum conversion value must be a scalar logical.';
               end

               if (islogical(convert))
                   args.ConvertDatenum2 = convert;
               elseif (isnumeric(convert))
                   args.ConvertDatenum2 = logical(convert);
               else
                   msg = 'Datenum conversion value must be a scalar logical.';
               end
           end
        case 'epochiscdfepoch'
           if (k == length(varargin))
               msg = 'No EpochIsCDFEpoch value specified.';
               return
           else
               noconvert = varargin{k + 1};
               if (numel(noconvert) ~= 1)
                   msg = 'EpochIsCDFEpoch value must be a scalar logical.';
               end
               if (islogical(noconvert))
                   args.EpochIsCDFEpoch = noconvert;
               elseif (isnumeric(noconvert))
                   args.EpochIsCDFEpoch = logical(noconvert);
               else
                   msg = 'Datenum conversion value must be a scalar logical.';
               end
           end
        case 'tt2000'
           if (k == length(varargin))
               msg = 'No TT2000 value specified.';
               return
           else
               noconvert = varargin{k + 1};
               if (numel(noconvert) ~= 1)
                   msg = 'TT2000 value must be a scalar logical.';
               end
               if (islogical(noconvert))
                   args.TT2000 = noconvert;
               elseif (isnumeric(noconvert))
                   args.TT2000 = logical(noconvert);
               else
                   msg = 'TT2000 conversion value must be a scalar logical.';
               end
           end
        end  % switch
    end  % for
    
    % Do a sanity check on the sizes of what we are passing back
    if ~isequal(length(args.VarNames), length(args.VarVals), ...
                length(args.PadVals))
        exception.msg = 'Number of variable names, values, and pad values do not match.';
        exception.id = 'MATLAB:cdfwrite:sanityCheckMismatch';
        return    
    end
    if ~isequal(length(args.VarNames), length(args.VarVals))
        exception.msg = 'Number of variable names and values do not match.';
        exception.id = 'MATLAB:cdfwrite:sanityCheckMismatch';
        return    
    end
    validate_inputs(args);

end  % if (nargin > 1)

function validate_inputs(args)
%VALIDATE_INPUTS   Ensure that the mutually exclusive options weren't provided.

if ((args.TT2000) && (args.EpochIsCDFEpoch))
    error('MATLAB:cdfwrite:TT2000', '%s\n%s', ...
          'You cannot currently specify these two options.', ...
          'Specify only one of ''TT2000'' and ''EpochIsCDFEpoch''.')
end

if ((args.TT2000) && (args.ConvertDatenum))
    error('MATLAB:cdfwrite:TT2000', '%s\n%s', ...
          'You cannot currently specify these two options.', ...
          'Specify only one of ''TT2000'' and ''ConvertDatenumToEpoch''.')
end

if ((args.ConvertDatenum) && (args.ConvertDatenum2))
    error('MATLAB:cdfwrite:TT2000', '%s\n%s', ...
          'You cannot currently specify these two options.', ...
          'Specify only one of ''ConvertDatenumToEpoch'' and ''ConvertDatenumToTT2000''.')
end

