function info = spdfcdfinfo(filename, varargin)
%SPDFCDFINFO Get details about a CDF file.
%   INFO = SPDFCDFINFO(FILE) gives information about a Common Data Format
%   (CDF) file.  INFO is a structure containing the following fields:
%
%     Filename             A string containing the name of the file
%
%     FileModDate          A string containing the modification date of
%                          the file
%
%     FileSize             An integer indicating the size of the file in
%                          bytes
%
%     Format               A string containing the file format (CDF)
%
%     FormatVersion        A string containing the version of the CDF
%                          library used to create the file
%
%     FileSettings         A structure containing the file settings
%                          describing the file
%
%     Subfiles             A cell array of filenames which contain the
%                          CDF file's data if it is a multifile CDF
%
%     Variables            A cell array containing details about the
%                          variables in the file (see below)
%
%     GlobalAttributes     A structure containing the global meta-data
%
%     VariableAttributes   A structure containing meta-data for the
%                          variables
%
%     LibVersion           A string containing the library version of 
%                          the CDF that is used to build the SPDFCDFinfo tool
%
%     PatchVersion         A string containing the patch version of 
%                          the MATLAB-CDF modules, e.g., spdfcdfinfo,
%                          spdfcdfread, spdfcdfwrite, and its release
%
%   Note: The CDF file name can be passed in as a null. In this case, INFO
%         contains only the CDF library version and MATLAB-CDF version.  
%
%   The "Variables" field contains a cell array of details about the
%   variables in the CDF file.  There are two presentations that can be
%   acquired. The original form is presented as the following:
%   each row represents a variable in the file.  The columns are:
%
%     (1) The variable's name as a string.
%
%     (2) The dimensions of the variable according to MATLAB's SIZE
%         function. 
%
%     (3) The number of records assigned for this variable.
%
%     (4) The variable's data type as it is stored in the CDF file.
%
%     (5) The record and dimension variance settings for the variable.
%         The value to the left of the slash designates whether values
%         vary by record; the values to the right designate whether
%         values vary at each dimension.
%
%     (6) The sparsity of the variables records.  Allowable values are
%         'Full', 'Sparse(padded)', and 'Sparse(previous)'.
%
%     (7) The compression of the variables.  Allowable values are
%         'None', 'RLE', 'HUFF', 'AHUFF', and 'GZIP.x' where x is the 
%         GZIP compression level.
%
%     (8) The blocking factors of the variables, if set.  A blocking factor is
%         the chunk of records to be pre-allocated for a standard variable when
%         a record being written does not already exist, or the number of
%         records to be compressed together for a compressed variable. For 
%         variables with large record size or number, providing a large 
%         blocking factor would have a significant impact on I/O performance
%         over the defaults.
%
%     (9) The pad values of the variables, if set. 
%    (10) The FILLVAL attribute entry values of the variables, if set. 
%    (11) The VALIDMIN attribute entry values of the variables, if set. 
%    (12) The VALIDMAX attribute entry values of the variables, if set. 
%
%   varinfo = spdfcdfinfo (FILE, 'VARIABLES', {'var1', 'var2', ...});
%   The optional 'VARIABLES' can be used to specify the variables that their
%   info is to be returned. The info includes only two varibale related fields:
%   Variables and VariableAttributes. No CDF info, e.g., Filename, Format, etc,
%   will be retrieved. A null is filled if a variable is not found in the CDF.
%
%   info = spdfcdfinfo (FILE, 'VARSTRUCT', TF);
%   The returned Variable field can also be presented in a structure form, if 
%   the option is provided with a true value.
%
%   The structure has the following fields, each one is a cell array.
%
%     Name                 A string containing the name of the variable
%
%     Dimensions           The dimensions of the variable
%
%     NumRecords           The number of written records for the variable
%
%     DataType             The data type of the variable
%
%     RecDimVariance       Record and dimensional variances of the variable
%
%     Sparseness           The sparseness of the variable
%
%     Compression          The Compression of the variablee
%
%     BlockingFactor       The blocking factor of the variable
%
%     PadValue             The pad value of the variable
%
%     FILLVAL              The FILLVAL attribute entry value for the variable
%
%     VALIDMIN             The VALIDMIN attribute entry value for the variable
%
%     VALIDMAX             The VALIDMAX attribute entry value for the variable
%
%   info = spdfcdfinfo (FILE, 'VALIDATE', TF);
%   This is to specify whether the CDF is to be validated when it's open. The
%   default is NOT to valdate the file so the processing can be faster. There
%   are two ways to set the data validation: setting the environment variable
%   CDF_VALIDATE to "yes" outside of the MATLAB environment, or using the
%   option 'VALIDATE' with true value when calling this module. If a CDF has 
%   been validated before, there is no need to validate it over and over again.
%
%   The "GlobalAttributes" and "VariableAttributes" structures contain a
%   field for each attribute.  Each field's name corresponds to the name
%   of the attribute, and the field contains a cell array containing the
%   entry values for the attribute.  For variable attributes, the first
%   column of the cell array contains the Variable names associated with
%   the entries, and the second contains the entry values.
%
%   NOTE: Attribute names which SPDFCDFINFO uses for field names in
%   "GlobalAttributes" and "VariableAttributes" may not match the names
%   of the attributes in the CDF file exactly.  Because attribute names
%   can contain characters which are illegal in MATLAB field names, they
%   may be translated into legal field names.  Illegal characters which
%   appear at the beginning of attributes are removed; other illegal
%   characters are replaced with underscores ('_').  If an attribute's
%   name is modified, the attribute's internal number is appended to the
%   end of the field name.  For example, '  Variable%Attribute ' might
%   become 'Variable_Attribute_013'.
%
%   To get the CDF library version that is used to build the current MATLAB
%   tool programs, e.g., spdfcdfread, spdfcdfwrite, spdfcdfinfo, you can
%   simply enter the command without providing a CDF file:
%   spdfcdfinfo()
%
%   Library version may differ from a CDF file's version. A CDF file is
%   assigned with the library version when it is created or modified.
%
%   Notes:
%
%     SPDFCDFINFO creates temporary files when accessing CDF files.  The
%     current working directory must be writable.
%
%
%   See also SPDFCDFREAD, SPDFCDFUPDATE, SPDFCDFWRITE, CDFEPOCH, CDFTT2000,
%            SPDFENCODEEPOCH, SPDFCOMPUTEEPOCH, SPDFPARSEEPOCH,
%            SPDFBREAKDOWNEPOCH, SPDFENCODEEPOCH16, SPDFCOMPUTEEPOCH16,
%            SPDFPARSEEPOCH16, SPDFBREAKDOWNEPOCH16, SPDFENCODETT2000,
%            SPDFCOMPUTETT2000, SPDFPARSETT2000, SPDFBREAKDOWNTT2000,
%            SPDFDATENUMTOEPOCH, SPDFDATENUMTOEPOCH16, SPDFDATENUMTOTT2000,
%            SPDFCDFLEAPSECONDSINFO

%
% HISTORY:
%   August 17, 2007   David Han      Modified to handle CDF_EPOCH16. Look for
%                                    'epoch16'.
%   July 17, 2009     Mike Liu       Added a new field 'LibVersion' for the
%                                    returned structure to show the library 
%                                    version. 'PatchVersion' will show the
%                                    patch version of the MATLAB release.
%   August 10, 2010   Mike Liu       Added INT8 and TT2000 data types.
%
%   November 14, 2019 Mike Liu       Added the latest leap second info from the
%                                    table to 'LibVersion'.

%
% Process arguments.
%

%if (argin < 1)
%  if ~(nargin == 0) || ~(length(strtrim(filename)) == 0)
%    error('MATLAB:spdfcdfinfo:inputArguments', 'SPDFCDFINFO requires at least one input argument.')
%  end
%end

if ~(nargout == 1)
  if ~(nargin == 0) && ~(length(strtrim(filename)) == 0)
    error('MATLAB:spdfcdfinfo:outputArguments', 'SPDFCDFINFO requires one output argument.')
  end
end

% CDFlib creates temporary files in the current directory.  Make sure PWD is
% writable.
%[attrib_success, attrib_mode] = fileattrib(pwd);
%
%if (~attrib_mode.UserWrite)
%    error('Cannot create temporary files.  The current directory must be writable.')
%end

%
% Returned structure
%

if (ispc)
% PC: set characterset to UTF8 before MATLAB 2021a versions 
  if (cstrcmp(version('-release'), '2021a') < 0)
    if (cstrcmp(feature('DefaultCharacterSet'), 'UTF-8') ~= 0)
      feature('DefaultCharacterSet', 'UTF8');
    end
  end
elseif (ismac)
% Mac: set characterset to UTF8 before MATLAB 2020a versions 
  if (cstrcmp(version('-release'), '2020a') < 0)
    if (cstrcmp(feature('DefaultCharacterSet'), 'UTF-8') ~= 0)
      feature('DefaultCharacterSet', 'UTF8');
    end
  end
end

info1.Filename = '';
info1.FileModDate = '';
info1.FileSize = '';
info1.Format = '';
info1.FormatVersion = '';
info1.FileSettings = [];
info1.FileSettings = {};
info1.Subfiles = {};
info1.Variables = {};
info1.GlobalAttributes = [];
info1.VariableAttributes = [];
info1.LibVersion = '';
info1.PatchVersion = '3.9.0.0';
info2.Variables = {};
info2.VariableAttributes = [];
args.VarStruct = false;
args.Variables = {};
args.Validate = false;
variables = false;

if (nargin == 0) || (length(strtrim(filename)) == 0)
    % Only for library info
    tmp = spdfcdfinfoc(' ');
    % Library version.
    theLibVersion = sprintf('%d.%d.%d', tmp.LibVersion.Version, ...
                                        tmp.LibVersion.Release, ...
                                        tmp.LibVersion.Increment);
    theLibLeapSecond = sprintf('%d', tmp.LibVersion.Latest_Leapsecond);
    disp(['LibVersion: ',theLibVersion]);
    disp(['Latest leap Second from Table: ',theLibLeapSecond]);
    disp(['PatchVersion: ', info1.PatchVersion]);
else    
  % Parse arguments based on their number.
  if (nargin > 0)
    paramStrings = {'variables'
                    'validate'
                    'varstruct'};

    % For each pair
    for k = 1:2:length(varargin)
       param = lower(varargin{k});

       idx = strmatch(param, paramStrings);

       if (isempty(idx))
           msg = sprintf('Unrecognized parameter name "%s".', param);
           return
       elseif (length(idx) > 1)
           msg = sprintf('Ambiguous parameter name "%s".', param);
           return
       end

       switch (paramStrings{idx})
         case 'varstruct'
               varstruct = varargin{k + 1};
               if (numel(varstruct) ~= 1)
                   msg = 'The "varstruct" value must be a scalar logical.';
               end

               if (islogical(varstruct))
                   args.VarStruct = varstruct;
               elseif (isnumeric(varstruct))
                   args.VarStruct = logical(varstruct);
               else
                   msg = 'The "varstruct" value must be a scalar logical.';
               end

         case 'validate'
               validate = varargin{k + 1};

               if (islogical(validate))
                   args.Validate = validate;
               elseif (isnumeric(validate))
                   args.Validate = logical(validate);
               else
                   msg = 'The "validate" value must be a scalar logical.';
               end

         case 'variables'
           if (k == length(varargin))
               msg = 'No variables specified.';
               return
           else
               args.Variables = varargin{k + 1};
               if (~iscell(args.Variables))
                   args.Variables = {args.Variables};
               end
               for p = 1:length(args.Variables)
                   if (~ischar(args.Variables{p}))
                       msg = 'All variable names must be strings.';
                       return
                   end
               end
	       if (length(args.Variables) > 0)
		 variables = true;
	       end
           end

       end  % switch
    end  % for
  end
  % Get full filename.
  fid = fopen(filename);

  %
  % Verify existence of filename.
  %
  if (fid == -1)
  
  % Look for filename with extensions.
      fid = fopen([filename '.cdf']);
    
      if (fid == -1)
          fid = fopen([filename '.CDF']);
      end
    
  end

  if (fid == -1)
      error('MATLAB:spdfcdfinfo:fileOpen', 'Couldn''t open file (%s).', filename)
  else
      filename = fopen(fid);
      fclose(fid);
  end
  if ~variables
    %
    % Record the file details.
    %
    d = dir(filename);
    % Set the positions of the fields.
    info1.Filename = d.name;
    info1.FileModDate = d.date;
    info1.FileSize = d.bytes;
    info1.Format = 'CDF';
  end
  % CDFlib's OPEN_ routine is flakey when the extension ".cdf" is used.
  % Strip the extension from the file before calling the MEX-file.

  if ((length(filename) > 4) && (isequal(lower(filename((end-3):end)), '.cdf')))
      filename((end-3):end) = '';
  end
  % Get the attribute, variable, and library details.

  tmp = spdfcdfinfoc(filename, args.VarStruct, args.Variables, args.Validate);

  if ~variables
    % Process file attributes.
    info1.FileSettings = parse_file_info(tmp.File);
    info1.FormatVersion = info1.FileSettings.Version;
    info1.FileSettings = rmfield(info1.FileSettings, 'Version');
    % Handle multifile CDF's.
    if isequal(info1.FileSettings.Format, 'Multifile')

      d = dir([filename '.v*']);
  
      for p = 1:length(d)
        info1.Subfiles{p} = d(p).name;
      end
    
    end
  end
  if ~args.VarStruct
     % Process variable table.
     vars = tmp.Variables;
     types = vars(:, 4);
     otypes = vars(:, 4);
     sp = vars(:, 6);
     pv = vars(:, 9);
     for p = 1:length(types)
       types{p} = find_datatype(types{p});
       sp{p} = find_sparsity(sp{p});
       pv{p} = find_padvalue(otypes{p}, pv{p});
     end
     vars(:, 4) = types;
     vars(:, 6) = sp;
     vars(:, 9) = pv;
     if ~variables
       info1.Variables = vars;
     else
       info2.Variables = vars;
     end
  else
     vars = tmp.Variables;
     types = vars.DataType;
     otypes = vars.DataType;
     sp = vars.Sparseness;
     pv = vars.PadValue;
     for p = 1:length(types)
       types{p} = find_datatype(types{p});
       sp{p} = find_sparsity(sp{p});
       pv{p} = find_padvalue(otypes{p}, pv{p});
     end
     vars.DataType = types;
     vars.Sparseness = sp;
     vars.PadValue = pv;
     if ~variables
       info1.Variables = vars;
     else
       info2.Variables = vars;
     end
  end
  % Assign rest.
  if ~variables
    info1.GlobalAttributes = tmp.GlobalAttributes;
  end
  if ~variables
    info1.VariableAttributes = tmp.VariableAttributes;
  else
    info2.VariableAttributes = tmp.VariableAttributes;
  end
  if ~variables
    info1.LibVersion = sprintf('%d.%d.%d', tmp.LibVersion.Version, ...
                                           tmp.LibVersion.Release, ...
                                           tmp.LibVersion.Increment);
  end
  if ~variables
    info = info1;
  else
    info = info2;
  end
end

function out = parse_file_info(in)

    out = in;

    % Format.
    if (in.Format == 2)
        out.Format = 'Multifile';
    else
        out.Format = 'Single-file';
    end

    % Encoding.
    out.Encoding = find_encoding(in.Encoding);

    % Majority.
    if (in.Majority == 1)
        out.Majority = 'Row';
    else
        out.Majority = 'Column';
    end

    % Version.
    out.Version = sprintf('%d.%d.%d', in.Version, in.Release, in.Increment);
    out = rmfield(out, {'Release', 'Increment'});

    % Compression.
    [comp_type, comp_param, comp_pct] = find_compression(in.Compression, ...
                                                  in.CompressionParam, ...
                                                  in.CompressionPercent);

    out.Compression = comp_type;
    out.CompressionParam = comp_param;
    out.CompressionPercent = comp_pct;

    % Checksum.
    if (in.Checksum == 1)
        out.Checksum = 'MD5';
    else
        out.Checksum = 'None';
    end

    % LeapSecondLastUpdated.
    if (in.LeapSecondLastUpdated >= 0)
      out.LeapSecondLastUpdated = in.LeapSecondLastUpdated;
    else
      out.LeapSecondLastUpdated = 'Not set';
    end

function str = find_datatype(num)
if (isempty(num))
  str = [];
else
  switch (num)
    case {1, 41}
      str = 'int8';
    case {2}
      str = 'int16';
    case {4}
      str = 'int32';
    case {8}
      str = 'int64';
    case {11}
      str = 'uint8';
    case {12}
      str = 'uint16';
    case {14}
      str = 'uint32';
    case {21, 44}
      str = 'single';
    case {22, 45}
      str = 'double';
    case {31}
      str = 'epoch';
    case {32}
      str = 'epoch16';
    case {33}
      str = 'tt2000';
    case {51, 52}
      str = 'char';
  end
end

function str = find_padvalue(num, value)
if (isempty(value)) str = [];
else
  switch (num)
    case {1, 41}
      str = int8(value);
    case {2}
      str = int16(value);
    case {11}
      str = uint8(value);
    case {12}
      str = uint16(value);
    case {21, 44}
      str = single(value);
    case {22, 45}
      str = double(value);
    case {51, 52}
      str = value;
    case {8, 33}
      str = int64(value);
    case {4}
      str = int32(value);
    case {14}
      str = uint32(value);
    case {31}
      str = double(value);
    case {32}
      str = value;
  end
end

function str = find_sparsity(num)
if (isempty(num))
  str = [];
else
  switch (num)
    case 0
      str = 'Full';
    case 1
      str = 'Sparse(padded)';
    case 2
      str = 'Sparse(previous)';
  end
end

function str = find_encoding(num)
if (isempty(num))
  str = [];
else
  switch (num)
    case 1
      str = 'Network';
    case 2
      str = 'Sun';
    case 3
      str = 'Vax';
    case 4
      str = 'DECStation';
    case 5
      str = 'SGI';
    case 6
      str = 'IBM-PC';
    case 7
      str = 'IBM-RS';
    case 8
      str = 'Host';
    case 9
      str = 'Macintosh';
    case 11
      str = 'HP';
    case 12
      str = 'NeXT';
    case 13
      str = 'Alpha OSF1';
    case 14
      str = 'Alpha VMS d';
    case 15
      str = 'Alpha VMS g';
    case 16
      str = 'Alpha VMS i';
    case 17
      str = 'ARM Little';
    case 18
      str = 'ARM Big';
    case 19
      str = 'IA64 VMS i';
    case 20
      str = 'IA64 VMS d';
    case 21
      str = 'IA64 VMS g';

  end
end



function [ctype, param, pct] = find_compression(ctype, param, pct)
if (isempty(ctype))
  ctype = [];
  param = [];
  pct = [];
else
  switch (ctype)
    case 0
      ctype = 'Uncompressed';
      param = '';
      pct = [];
    case 1
      ctype = 'Run-length encoding';
      if (param == 0)
          param = 'Encoding of zeros';
      else
          param = '';
      end
    case 2
      ctype = 'Huffman';
      if (param == 0)
          param = 'Optimal encoding trees';
      else
          param = '';
      end
    case 3
      ctype = 'Adaptive Huffman';
      if (param == 0)
          param = 'Optimal encoding trees';
      else
          param = '';
      end
    case 4
      ctype = 'Rice';
      param = '';
    case 5
      ctype = 'Gzip';
  end 
end


function c = cstrcmp(a,b)
% c style strcmp. No type checks, use only strings and char arrays.
if strcmp(a,b)
    c=0;
else
    [idx,idx]=sort({char(a),char(b)});
    c=idx(1)-idx(2);
end

