function spdfcdfupdate(filename, varargin)
%SPDFCDFUPDATE Update data in an existing CDF file.
% 
%   SPDFCDFUPDATE(FILE, 'VariableData', VARIABLELIST, ...) updates variable 
%   data in a CDF file whose name is specified by FILE.  VARIABLELIST is a 
%   cell array of ordered pairs, which are comprised of a CDF variable name 
%   (a string) and a corresponding value cell array. The value cell array 
%   is comprised of a number of cells, each cell consisting of a record number,
%   a dimension indices, as well as the new data value. The record number and 
%   dimension indices are all zero(0)-based. The number of dimension indices 
%   depends on the dimensionality of the variable. The form for entering
%   multiple value changes for variables is as follows:
%   {'varname1' [ {rec_no_1, {index1, index2, ...}, value1}; ...
%                 {rec_no_2, {index2, index2, ...}, value2}; ...
%                 ...
%               ] ...
%    'varname2' [ {rec_no_1, {index1, index2, ...}, value1}; ...
%                 {rec_no_2, {index2, index2, ...}, value2}; ...
%                 ...
%               ] ...
%    ...
%   }                 
%
%
%   SPDFCDFUPDATE(..., 'GlobalAttributes', GATTRIB, ...) updates the structure
%   GATTRIB as global meta-data for the CDF.  Each field of the
%   struct is the name of a global attribute.  The value of each
%   field contains the value of the attribute.  To update
%   multiple values for an attribute, the field value should be a
%   cell array.
%
%   SPDFCDFUPDATE(..., 'VariableAttributes', VATTRIB, ...) updates the
%   structure VATTRIB as variable meta-data for the CDF.  Each
%   field of the struct is the name of a variable attribute.  The
%   value of each field should be an mx2 cell array where m is the
%   number of variables with attributes.  The first element in the
%   cell array is the name of the variable, a string, and the 
%   second element, also a string, is the CDF specific data type of 
%   the entry, and the third element is the new entry value of the 
%   attribute for that variable.
%
%   SPDFCDFUPDATE(..., 'CDFCompression', TF, ...) resets the CDF compression.
%   The CDF is set to use GZIP compression if TF is true. If TF is
%   is false, the CDF is set to be a uncompressed file.
%
%   SPDFCDFUPDATE(..., 'CDFChecksum', TF, ...) resets the CDF checksum.
%   The CDF is set to use MD5 checksum if TF is true. If TF is
%   is false, the CDF is set not to use checksum.
%
%   SPDFCDFUPDATE(..., 'CDFLeapSecondLastUpdated', value, ...) resets the CDF's
%   leap second last updated date.  For CDFs created prior to V 3.6.0, this 
%   field is not set. It is set to indicate what leap second table this CDF is
%   based upon. The value, in YYYYMMDD form, must be a valid entry in the
%   currently used leap second table, or zero (0) if the table is not used. 
%   CDF will automatically fill the value with the last entry date in the leap
%   second table if this option is not specified. 
%
%   Notes:
%
%     SPDFCDFUPDATE only updates the data values for the existing variables.
%     For Epoch and Epoch16 data type variables, the updated values should
%     be in the string form, e.g., 01-Jan-2008 07:06:05.004 and 
%     01-Jan-2008 07:06:05.004:100:200:300:400, respectively, instead of
%     the double form as they are stored.
%
%   Examples:
%
%   %  Change data values for Longitude (a 1-dim of CDF_INT2 data type): 
%   %  4th element in Record 0 to 100, 3rd element in Record 1 to 200; and 
%   %  change data value for variable Latitude (a 2-dim of CDF_REAL4 data type) 
%   %  at index [0,0] in Record 10 to 20.2; and change record 2 for variable
%   %  Epoch (a 0-dim of CDF_EPOCH data type) to 01-Jan-2008 07:06:05.004 in 
%   %  the file 'example.cdf'. 
%
%      spdfcdfupdate('example', 'VariableData', {'Longitude' [{0,{3},100};{1,{2},200}] ...
%                                            'Latitude' [{10,{0,0},20.2}] ...
%                                            'Epoch' [{1,{0},'01-Jan-2008 07:06:05.004'}] ...
%                                           } ...
%               );
%   Or,
%      var1 = {'Longitude' [{0,{3},100};{1,{2},200}]};
%      var2 = {'Latitude' [{10,{0,0},20.2}]};
%      var3 = {'Epoch' [{1,{0},'01-Jan-2008 07:06:05.004'}]};
%      spdfcdfupdate('example', 'VariableData', [var1 var2 var3]);
%
%   Make sure that the data values are of the data type for the variables.
%
%   % Update the meta-data for the entry of variable attribute 'validmin' for
%   % variable 'Longitude' while updating its data values in the file 'example.cdf': 
%
%   varAttribStruct.validmin = {'Longitude' [10]};
%   spdfcdfupdate('example', 'VariableData', {'Longitude' [{0,{3},100};{1,{2},200}]}, ...
%                        'VariableAttributes', varAttribStruct);
%
%   % Update the global attribute, 'Name': the first entry to 'MyName', a CDF_CHAR 
%   % data type, and the second entry to 10, of data type CDF_INT2, in the 
%   % file 'example.cdf':
%
%   globalAttribStruct.Name = [{0, 'CDF_CHAR', 'MyName'},{1, 'CDF_INT2', 10}];
%   spdfcdfupdate('example', 'GlobalAttributes', globalAttribStruct);
%
%   % Set the CDF file, 'example.cdf', to use the GZIP compression:
%
%   spdfcdfupdate('example', 'cdfcompression', true);
%
%   % Set the CDF file, 'example.cdf', to use the MD5 checksum:
%
%   spdfcdfupdate('example', 'cdfchecksum', true);
%
%   See also SPDFCDFREAD, SPDFCDFWRITE, SPDFCDFINFO, CDFEPOCH, CDFTT2000,
%            SPDFENCODEEPOCH, SPDFCOMPUTEEPOCH, SPDFPARSEEPOCH,
%            SPDFBREAKDOWNEPOCH, SPDFENCODEEPOCH16, SPDFCOMPUTEEPOCH16,
%            SPDFPARSEEPOCH16, SPDFBREAKDOWNEPOCH16, SPDFENCODETT2000,
%            SPDFCOMPUTETT2000, SPDFPARSETT2000, SPDFBREAKDOWNTT2000,
%            SPDFDATENUMTOEPOCH, SPDFDATENUMTOEPOCH16, SPDFDATENUMTOTT2000,
%            SPDFCDFLEAPSECONDSINFO

% HISTORY:
%   January 13, 2009  Mike Liu      The new function was created.
%   August  23, 2013  Mike Liu      Added setting checksum operation.

%
% Process arguments.
%

if (nargin < 2)
    error('MATLAB:spdfcdfupdate:inputArgumentCount', ...
          'CDFWRITE requires at least two input arguments.')
end

% parse_inputs sorts out all of the input args.  Its return values:
%
% * args - an array of structs.  args.VarNames contains the names
%          of the variables to be updated to the CDF.  args.VarVals contains
%          the cell array for updating the data values.  
% * msg - an error message from parse_inputs that we pass on to the user.

[args, varAttribStruct, globalAttribStruct, exception] = parse_inputs(varargin{:});

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
% Call the underlying spdfcdfupdatec function which calls the CDFlib
%

spdfcdfupdatec(filename, args.VarNames, args.VarRecs, args.VarIndices, ...
           args.VarDataVals, varAttribStruct, globalAttribStruct, ...
           args.isCDFCompressed, args.CDFchecksum, ...
           args.CDFleapsecondlastupdated);

%%%
%%% Function parse_inputs
%%%

function [args, varAttribStruct, globalAttribStruct, exception] = parse_inputs(varargin)

% Set default values
args.VarNames = {};
args.VarRecs = {};
args.VarIndices = {};
args.VarDataVals = {};
varcell = {};
varAttribStruct = struct([]);
globalAttribStruct = struct([]);
exception.msg = '';
exception.id = '';
args.isCDFCompressed = int32(0);
args.CDFchecksum = int32(0);
args.CDFleapsecondlastupdated = int32(-999);

% Parse arguments based on their number.
if (nargin > 0)
    
    paramStrings = {'variabledata'
                    'globalattributes'
                    'variableattributes'
                    'cdfcompression'
                    'cdfchecksum'
                    'cdfleapsecondlastupdated'};
    % For each pair
    for k = 1:2:length(varargin)
        param = lower(varargin{k});
            
        if (~ischar(param))
            exception.msg = 'Parameter name must be a string.';
            exception.id = 'MATLAB:spdfcdfupdate:paramNameNotString';
            return
        end
        
        idx = strmatch(param, paramStrings);
        
        if (isempty(idx))
            exception.msg = sprintf('Unrecognized parameter name "%s".', param);
            exception.id = 'MATLAB:spdfcdfupdate:unrecognizedParam';
            return
        elseif (length(idx) > 1)
            exception.msg = sprintf('Ambiguous parameter name "%s".', param);
             exception.id = 'MATLAB:spdfcdfupdate:ambiguousParam';
           return
        end
        
        switch (paramStrings{idx})
        case 'globalattributes'
           globalAttribStruct = varargin{k+1};
           if ~isstruct(globalAttribStruct)
               exception.msg = ['''GlobalAttributes''' ' must be a struct.'];
               exception.id = 'MATLAB:spdfcdfupdate:globalAttributesNotStruct';
               return
           end
           attribs = fieldnames(globalAttribStruct);
           
           % If the global attribute isn't a cell, then stuff it in one.
           for i = 1:length(attribs)
               attribVal = globalAttribStruct.(attribs{i});
               s = size(attribVal);
               if ~iscell(attribVal)
                   globalAttribStruct.(attribs{i}) = {attribVal};
               end
               % The global attribute struct may have more than one set of
               % entry id, data type and value.  However, there must only be
               % one associated value of the attribute for each entry,
               % hence the 3.
               if (s(2) == 3)
                   % Transpose it because CDFlib reads the arrays column-wise.
                   globalAttribStruct.(attribs{i}) = attribVal';   
               else
                   % We have ordered pairs.
                   globalAttribStruct.(attribs{i}) = reshape(globalAttribStruct.(attribs{i})(:),numel(globalAttribStruct.(attribs{i})(:))/3, 3);
               end
           end
        case 'variabledata'
           varcell = varargin{k+1};
           if ~iscell(varcell)
               exception.msg = ['''VariableData''' ' must be a cell.'];
               exception.id = 'MATLAB:spdfcdfupdate:variabledatanotcell';
               return
           end
           % First check that varcell meets all of our requirements
           args.VarNames = {varcell{1:2:end}};
           args.VarVals = {varcell{2:2:end}};

           if length(args.VarNames) ~= length(args.VarVals)
               exception.msg = 'All variable names must have a corresponding variable value cell array.';
               exception.id = 'MATLAB:spdfcdfupdate:variableWithoutValue';
               return
           end
           % args.VarVals
           % Check and make sure that all variable values are of the same
           % datatype, but ignore empties
           if ~isempty(args.VarVals)
               for i = 1:length(args.VarVals)
                   a = args.VarVals{i};
                   if (~iscell(a))
                     a = {a};
                   end
           	jj = numel(a)/length(a);
           	args.VarRecs{i} = {a{1:jj}};
           	args.VarIndices{i} = {a{jj+1:2*jj}};
           	args.VarDataVals{i} = {a{2*jj+1:3*jj}};
               end
           end
        case 'variableattributes'
           varAttribStruct = varargin{k+1};
           if ~isstruct(varAttribStruct)
               exception.msg = ['''VariableAttributes''' ' must be a struct.'];
               exception.id = 'MATLAB:spdfcdfupdate:variableAttributesNotStruct';
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
               % variable-value pair.  However, there must only be
               % one associated value of the attribute for each variable,
               % hence the 2.
               if (s(2) == 2)
                   % Transpose it because CDFlib reads the arrays column-wise.
                   varAttribStruct.(attribs{i}) = attribVal';   
               else
                   % We have ordered pairs.
                   varAttribStruct.(attribs{i}) = reshape(varAttribStruct.(attribs{i})(:),numel(varAttribStruct.(attribs{i})(:))/2, 2);
               end
           end
        case 'cdfcompression'
           args.isCDFCompressed = varargin{k+1};
           if (numel(args.isCDFCompressed) ~= 1)
               exception.msg = 'CDF compression value must be a scalar logical.';
               exception.id = 'MATLAB:spdfcdfupdate:cdfcompression';
               return
           end

           if (islogical(args.isCDFCompressed))
               if (args.isCDFCompressed)
                 args.isCDFCompressed = int32(2);
               else
                 args.isCDFCompressed = int32(1);
               end
           elseif (isnumeric(args.isCDFCompressed))
%              args.isCDFCompressed = logical(args.isCDFCompressed);
               if (logical(args.isCDFCompressed))
                 args.isCDFCompressed = int32(2);
               else
                 args.isCDFCompressed = int32(1);
               end
           else
               exception.msg = 'CDF compression value must be a scalar logical.';
               exception.id = 'MATLAB:spdfcdfupdate:cdfcompression';
               return
           end
        case 'cdfchecksum'
           args.CDFchecksum = varargin{k+1};
           if (numel(args.CDFchecksum) ~= 1)
               exception.msg = 'CDF checksum value must be a scalar logical.';
               exception.id = 'MATLAB:spdfcdfupdate:cdfchecksum';
               return
           end
           if (islogical(args.CDFchecksum))
               if (args.CDFchecksum)
                 args.CDFchecksum = int32(2);
               else
                 args.CDFchecksum = int32(1);
               end
           end
        case 'cdfleapsecondlastupdated'
           args.CDFleapsecondlastupdated = varargin{k+1};
           if (isnumeric(args.CDFleapsecondlastupdated))
               args.CDFleapsecondlastupdated = int32(args.CDFleapsecondlastupdated);
           else
               exception.msg = 'CDF leapsecondlastupdated value must be a numeric value (in YYYYMMDD form).';
               exception.id = 'MATLAB:spdfcdfupdate:cdfleapsecondlastupdated';
               return
           end
        end % switch
    end  % for
    
end

