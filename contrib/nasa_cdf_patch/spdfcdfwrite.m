function spdfcdfwrite(filename, varcell, varargin)
%SPDFCDFWRITE Write data to a CDF file.
% 
%   SPDFCDFWRITE(FILE, VARIABLELIST, ...) writes out a CDF file whose name
%   is specified by FILE.  VARIABLELIST is a cell array of ordered
%   pairs, which are comprised of a CDF variable name (a string) and
%   the corresponding CDF variable value.  To write out multiple records
%   for a variable, there are two ways of doing it. One way is putting the
%   variable values in a cell array, where each element in the cell array
%   represents a record. Another way, the better one, is to place the
%   values in a vector (single or multi-dimensional) with the option
%   'RecordBound' being specified.
%   Note: For variable with string data, make sure all strings are of equal
%         length.
%
%   SPDFCDFWRITE(..., 'PadValues', PADVALS) writes out pad values for given
%   variable names.  PADVALS is a cell array of pairs of a variable name (a
%   string) and a corresponding pad value.  Pad values are the default value
%   associated with the variable when an out-of-bounds record is accessed.
%   Variable names that appear in PADVALS must be in VARIABLELIST.
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
%   'VariableAttributes' field from the returned structure can be passed
%   in for the VATTRIB.
%   Note: For string variable attributes, they can be either a single string
%         or a cell of strings. 
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
%   SPDFCDFWRITE(..., 'WriteMode', MODE, ...) where MODE is either 'overwrite'
%   or 'append' indicates whether or not the specified variables or attributes
%   should be appended to the CDF if the file already exists.  The 
%   default is 'overwrite', indicating that SPDFCDFWRITE will not append
%   variables and attributes.
%
%   SPDFCDFWRITE(..., 'Format', FORMAT, ...) where FORMAT is either 'multifile'
%   or 'singlefile' indicates whether or not the data is written out
%   as a multi-file CDF.  In a multi-file CDF, each variable is stored
%   in a *.vN file where N is the number of the variable that is
%   written out to the CDF.  The default is 'singlefile', which indicates
%   that SPDFCDFWRITE will write out a single file CDF.  When the 'WriteMode'
%   is set to 'Append', the 'Format' option is ignored, and the format
%   of the pre-existing CDF is used.
%
%   SPDFCDFWRITE(..., 'Version', VERSION, ...) where VERSION is a string which 
%   specifies the version of the CDF library to use in writing the file.
%   The default option is to use the latest version of the library 
%   (currently version 3.2+), and has to be specified '3.0'.  The 
%   other available version is version 2.7 ('2.7').  Note that 
%   versions of MATLAB before R2006b will not be able to read files 
%   which were written with CDF versions greater than 3.0.
%
%   SPDFCDFWRITE(..., 'RecordBound', RECBNDVARS) specifies data values in arrays
%   (1-D or multi-dimensional) are to be written into "records" for the given
%   variable. RECBNDVARS is a cell array of variable names. The M-by-N array
%   data will create M rows (records), while each row having N elements. For an
%   M-by-1 (column) or 1-by-M (row) vector, it will create M records, each 
%   being a scalar. For a 3-D array M-by-N-by-R, R records will be written,
%   and each record with M-by-N elements. Without this option, an array of 
%   M-by-N will be written into a single record of 2-dimensions. See sample
%   codes for its usage.
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
%   SPDFCDFWRITE(..., 'VarCompress', COMPRESSVARS) specifies which variables
%   will be compressed by what compression method (and level for GZIP).
%   COMPRESSVARS is a cell array of variable names and their respective
%   compression methods. For examples, to set the compression for 'var1' and
%   'var2' with AHUFF and GZIP.6 respectively, COMPRESSVARS should be given
%   as {'var1', 'ahuff', 'var2', 'gzip.6'}. GZIP.6 compression provides the
%   best compression rate and performance.
%
%   If there is a master CDF that has all the same variable info as the new CDF,
%   then SPDFCDFINFO module can be used to retrieve the infomation. The 
%   'Variables' field from the returned structure contain the compression info
%   (at element 7) for each variable,  Set up a cell to use such compression
%   info.
%
%   SPDFCDFWRITE(..., 'CDFCompress', COMPRESSVALUE) specifies the CDF will be
%   compressed at the file level by what compression method (and level for GZIP).
%   COMPRESSVALUE is a string of a valid compression method and level.
%
%   If there is a master CDF that can be used to define the new CDF's settings,
%   then SPDFCDFINFO module can be used to retrieve the infomation. The 
%   'Compression' and 'CompressionParam' fields from the returned structure's
%   'FileSettings' structure field contain the compression info. Combining these 
%   fields will make a valid compression. See the example.
%   Alternatively, COMPRESSVALUE can be provided as one of the following strings:
%   'none', 'gzip.x', 'rle' (Run-length encoding), 'huff' (Huffman), 'ahuff'
%   (Adaptive Huffman), where 'x' in 'gzip.x' is the level of 1-9 (6 being the
%   preferable). 
%
%   SPDFCDFWRITE(..., 'VarSparse', SPARSEVARS) specifies which variables have
%   sparse records. SPARSEVARS is a cell array of orderly pairs of a variable
%   name and its respective value. The valid value should be 'full',
%   'Sparse(padded)' or 'Sparse(previous)' (the 6th column from Variables field
%   from spdfcdfinfo).
%   Value 'full' is the default if a variable does not have the sparse records.
%   For examples, to set the sparse record for 'var1' and 'var2' with 'full'
%   and 'Sparse(previous)', SPARSEVARS should be specified as {'var1', 'full',
%   'var2', 'Sparse(padded)'}. (Actually, there is no need for 'var1'.) 
%   Variable's record data should be provided in cell form for the sparse record
%   variables. Data for any virtual records should be presented as []s. See the
%   sample code for its usage.  
%
%   SPDFCDFWRITE(..., 'BlockingFactor', BFVARS) specifies the blocking factor
%   for the variables. The blocking factor is the number of variable records
%   that will be pre-allocated when a new record is to be written. The default 
%   value could be too small for variables that have large record size or
%   number. A large blocking factor (up to a variable's maximum record number)
%   can make access to the full variable data more efficient as less blocks
%   are involved (thus less I/Os). It can also make the file less fragmented.
%   BFVARS is a cell array of pairs of variable name and its respective value.
%   The value should be a numeric.
%
%   If there is a master CDF that has all the same variable info as the new CDF,
%   then SPDFCDFINFO module can be used to retrieve the infomation. The 
%  'Variables' field from the returned structure contain the record sparseness
%   info (at element 6) for each variable,  Set up a cell to use such sparseness
%   info.
%
%   SPDFCDFWRITE(..., 'ConvertDatenumToEpoch', TF, ...) converts MATLAB datenum
%   values to CDF epoch data if TF is true. This option works with the
%   variable(s) that is of CDF_EPOCH type in a CDF. 
%   There are two ways to write data for epoch variable(s) of CDF_EPOCH into
%   a CDF. First, uses cdfepoch objects, each of which is from a datenum,
%   datestr, or even cdfepoch object, to pass data to spdfcdfwrite: 
%      SPDFCDFWRITE(..., {'Epoch',cdfepoch([...], ...)}, ...) 
%   This option is time and space consuming if large datasets are involved. 
%   The second way uses 'ConvertDatenumToEpoch':
%     SPDFCDFWRITE(..., {'Epoch',[...], ...}, 'ConvertDatenumToEpoch', TF, ...) 
%   If TF is set to true, the passed data are assumed as MATLAB datenum
%   values and a data conversion to CDF epoch values is performed. 
%   Setting it to false (the default), All data will be considered already in
%   CDF_EPOCH form and will be filled, as is, to
%   CDF_EPOCH data type variable(s). The CDF_EPOCH data need to be numeric of
%   mxDouble_CLASS (double).
%
%   SPDFCDFWRITE(..., 'ConvertDatenumToTT2000', TF, ...) converts MATLAB datenum
%   values to CDF TT2000 data if TF is true. This option works with the
%   variable(s) that is of CDF_TIME_TT2000 type in a CDF.
%   There are two ways to write data for epoch variable (s) of CDF_TIME_TT2000
%   into a CDF. First, uses cdftt2000 objects, each of which is from a datenum,
%   datestr, or even cdftt2000 object, to pass data to spdfcdfwrite:
%      SPDFCDFWRITE(..., {'Epoch',cdftt2000([...], ...)}, ...)
%   This option is also time and space consuming if large datasets are involved.
%   The second way uses 'ConvertDatenumToTT2000':
%     SPDFCDFWRITE(..., {'Epoch',[...], ...}, 'ConvertDatenumToTT2000', TF, ...)
%   If TF is set to true, the passed data are assumed as MATLAB datenum
%   values and a data conversion to CDF TT2000 values is performed.
%   Setting it to false (the default), All data will be considered already in
%   CDF_TIME_TT2000 form and will be filled, as is, to
%   CDF_TIME_TT2000 data type variable(s). The CDF TT2000 values needs to be
%   numeric of mxINT64_CLASS (int64).
%
%   SPDFCDFWRITE(..., 'Checksum', CHECKSUMVAL, ...) specifies whether the output
%   CDF should have its checksum computed. The valid values are 'MD5' or 'none'.
%   The default is 'none'.
%
%   SPDFCDFWRITE(..., 'CDFLeapSecondLastUpdated', value, ...) resets the CDF's
%   leap second last updated date.  For CDFs created prior to V 3.6.0, this 
%   field is not set. It is set to indicate what leap second table this CDF is
%   based upon. The value, in YYYYMMDD form, must be a valid entry in the
%   currently used leap second table, or zero (0) if the table is not used. 
%   CDF will automatically fill the value with the last entry date in the leap
%   second table if this option is not specified. 
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
%
%   Notes:
%
%     SPDFCDFWRITE creates temporary files when writing CDF files.  Both the
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
%   spdfcdfwrite('example', {'Longitude', 0:360});
%
%   % Write out a file 'example.cdf', containing a variable 'Longitude'
%   % with the value [0:360] (one record of 361 values), and with a variable
%   % attribute of 'validmin' with the value 10:
%
%   varAttribStruct.validmin = {'Longitude' [10]};
%   spdfcdfwrite('example', {'Longitude' 0:360}, ...
%                'VariableAttributes', varAttribStruct);
%
%   % Write out a file 'example.cdf' containing variables 'Longitude'
%   % and 'Latitude' with the variable 'Longitude' being a Record-bound
%   % (361 records to be written)::
%
%   spdfcdfwrite('example', {'Longitude', (0:360), 'Latitude', 10:20}, ...
%                'RecordBound', {'Longitude'});
%
%   % These two commands should write out the data values identically:
%      SPDFCDFWRITE('example', {'Epoch', num2cell(1:100)}, ...
%                   'epochiscdfepoch', true);
%      SPDFCDFWRITE('example', {'Epoch', (1:100)'}, ...
%                   'recordbound', {'Epoch'}, ...
%                   'epochiscdfepoch', true);
%
%   % Write out a file 'example.cdf' containing variables 'Longitude'
%   % and 'Latitude' with the variable 'Latitude' having a pad value
%   % of 10 for all out-of-bounds records that are accessed. The
%   % 'Longitude' variable has one record (row) of a vector with 361 elements
%   % 'Latitude' has one record of a vector with 11 elements.
%
%   spdfcdfwrite('example', {'Longitude', (0:360)', 'Latitude', (10:20)'}, ...
%                'RecordBound', {'Latitude', 'Longitude'}, ...
%                'PadValues', {'Latitude', 10});
%
%   % Write out a file 'example.cdf', with multiple rows of time series data.
%   % Each row has a time and a sampled data, both being scalar. 
%   % The following sample shows 100 records (rows): epoch is of CDF_EPOCH 
%   % data type and Samples of CDF_DOUBLE. 
%   % The epoch starts from 2010-01-01T00:00:00.000 with 1 second stride.
%   % Each sampled value starts from 0.0, with 5.0 stride. 
%
%   epoch=utc2cdfepoch(2010,1,1,0,0,0,0);
%   epochvalues = (epoch+[0:99]*1000)';
%   values = ([0:99]*5.0)';
%   spdfcdfwrite('example', {'Epoch', epochvalues, 'Samples', values}, ...
%                'EpochisCDFEpoch', true, ...
%                'RecordBound', {'Epoch', 'Samples'});
%
%   % 'EpochisCDFEpoch' option is needed as 'Epoch' is to be created of
%   % CDF_EPOCH type.
%
%   % Alternatively, the same result can be accomplished by using the 
%   % 'EpochType' option:
%
%   spdfcdfwrite('example', {'Epoch', epochvalues, 'Samples', values}, ...
%                'EpochType', {'Epoch'}, ...
%                'RecordBound', {'Epoch', 'Samples'});
%
%   % Alternatively, the same result can be accomplished by making the epoch
%   % vector to cell. No 'RECORDBOUND' option is needed.
%   epoch=utc2cdfepoch(2010,1,1,0,0,0,0);
%   epochvalues = epoch+[0:99]*1000;
%   epochscell = num2cell(epochvalues);
%   values = num2cell([0:99]*5.0);
%   spdfcdfwrite('example', {'Epoch', epochvalues, 'Samples', values}, ...
%                'EpochType', {'Epoch'});
%
%   % Write out a file 'example.cdf', with single or multiple rows of
%   % vectorized data. Variable 'one0' will have one record with 5 elements.
%   % Variable 'one1' has five records, each record having 1 value. Variable
%   % 'two0' has one 5-by-2 record, while Variable 'two2' has five (5) 
%   % records, each record having 2 elements. Variable 'three0' has a
%   % single 3-D (3-by-2-by-2) record, while Variable 'three3' has two (2)
%   % records, each record being a 3-by-2 matrix. 
%
%   data0=1:5;
%   data1=data0';
%   data2=[10 20;30 40;50 60;70 80;90 100];
%   data2a=data2+100;
%   data3=[1 2;3 4;5 6];
%   data3(:,:,2)=[11 22; 33 44; 55 66];
%   data3a=data3+100;
%   spdfcdfwrite('example',{'one0',data0,'one1',data1,'two0',data2, ...
%                'two2',data2a,'three0',data3,'three3',data3a}, ...          
%                'recordbound',{'one1','two2','three3'});
%
%   % For writing out two variables: 'Epoch' of CDF_EPOCH type, and
%   % 'Sample' of CDF_DOUBLE type. Four records are written for each. Epoch's
%   % record is a scalar, while Sample's is an 1-D with 4 elements.
%
%   epoch=utc2cdfepoch(2010,1,1,0,0,0,0);
%   epochvalues = (epoch+[0:3]*1000)';
%   value = [0:3]*5.0;
%   values = [value; value+10; value+20; value+30];
%   spdfcdfwrite('example', {'Epoch', epochvalues, 'Sample', values}, ...
%                'EpochType', {'Epoch'}, ...
%                'RecordBound', {'Epoch', 'Sample'});
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
%   datenumvalues = (datenum1+[0:99]/86400)';
%   spdfcdfwrite('example', {'Epoch', datenumvalues}, ...
%                'ConvertDatenumToEpoch', true, ...
%                'RecordBound', {'Epoch'}, ...
%                'PadValue', {'Epoch', 1});
%
%   % Write out a file 'example.cdf', with three records for the variable
%   % 'Epoch' of CDF_TIME_TT2000 type. The first record has a date:
%   % 2010-10-10T01:02:03.456, the second record: 2010-11-11T02:04:06.789
%   % and the third record: 2010-12-12T03:04:05.123.
%
%   dates = datenum([2010 10 10 1 2 3.456; 2010 11 11 2 4 6.789; ...
%                    2010 12 12 3 4 5.123]);
%   spdfcdfwrite('example', {'Epoch', dates'}, ...
%                'ConvertDatenumToTT2000', true, ...
%                'RecordBound', {'Epoch'});
%
%   % Write out a file 'example.cdf', with 6 records for the variable 
%   % 'Epoch', which will be of CDF_TIME_TT2000 data type. The data
%   $ crosses over a leap second. Convert the epoch in UTC string
%   $ to their TT2000 values before write out to the CDF file.
%
%   time = {'2008-12-31T23:59:58.123456789'; ...
%           '2008-12-31T23:59:59.123456789'; ...
%           '2008-12-31T23:59:60.123456789'; ...
%           '2009-01-01T00:00:00.123456789'; ...
%           '2009-01-01T00:00:01.123456789'; ...
%           '2009-01-01T00:00:02.123456789'};
%   values = spdfparsett2000(time);          
%   spdfcdfwrite('example', {'Epoch', values}, ...
%                'EpochType', {'Epoch'}, ...
%                'RecordBound', {'Epoch'});
%
%   % Write out a file 'sample.cdf', with two variables. One variable, 'var2'.
%   % is compressed with GZIP.6.
%
%   var1data=......;
%   var2data=......;
%   spdfcdfwrite('sample',{'Epoch',var1data,'var2',var2data}, ...
%                'recordbound',{'Epoch','var2'}, ...
%                'varcompress',{'var2','gzip.6'});
%
%   % Write out a file 'sparse.cdf', with a variable that has sparse records.
%   % The variable 'one' only has two (2) physical data: at record 0 and
%   % record 4. 
%
%   spdata={[123 321];[];[];[];[-321 -123]};
%   spdfcdfwrite('sparse',{'one',spdata}, ...
%                'varsparse', {'one','Sparse(previous)'});
%
%   % Write out a file 'real.cdf', based on the master cdf, which provides the
%   % file settings info, i.e., checksum, CDF file level compression, as well as
%   % meta-data for all of the global and variable attribute information.  It
%   % also provides the variable spec, e.g., data type, record variance, record
%   % sparseness, blocking factor, pad value and compression, for each variable.
%   % 'Recordbound' option is used for specifying record-variant (RV) variables.
%   % Among the variables, Variable number 15, a single floating point array, has
%   % sparse records [at record 1 and 5]. Its data has to be in the cell array.
%
%   info=spdfcdfinfo('master.cdf');
%   for p = 1:length(info.Variables(:,1))
%     compress{(2*p)-1} = info.Variables{p,1}; 	% Variable name
%     compress{2*p} = info.Variables{p,7};	% Variable compression
%     sparse{(2*p)-1} = info.Variables{p,1};	% Variable name
%     sparse{2*p} = info.Variables{p,6};	% Variable sparseness
%     bf{2*p-1} = info.Variables{p,1};		% Variable name
%     bf{2*p} = info.Variables{p,8};		% Variable blocking factor
%     pad{2*p-1} = info.Variables{p,1};		% Variable name
%     pad{2*p} = info.Variables{p,9};		% Variable pad value
%     datatypes{2*p-1} = info.Variables{p,1};	% Variable name
%     datatypes{2*p} = info.Variables{p,4};	% Variable data type
%   end
%   rbvars = {info.Variables{:,1}};		% Variable names for recordbound
%   for p = length(rbvars):-1:1
%     if (strncmpi(info.Variables{p,5},'f',1)==1)	% NRV variable
%       rbvars(:,p)=[]; 	  		% Remove it
%     end
%   end
%   if isnumeric(info.FileSettings.CompressionParam) % A number for Gzip parameter 
%     cdfcompress=strcat(info.FileSettings.Compression, '.', ... % Make it 'gzip.x'
%                        num2str(info.FileSettings.CompressionParam));
%   else
%     cdfcompress=strcat(info.FileSettings.Compression, '.', ... % None or non-gzip
%                        info.FileSettings.CompressionParam);
%   end
%
%   % fill data
%   for p = 1:length(info.Variables(:,1))
%     varsdata{2*p-1} = info.Variables{p,1};
%     if (p == 15)				% A sparse record variable 
%       var15data={single([123 321]);[];[];[];single([-321 -123])};
%       varsdata{(2*15)} = var15data;		% Sparse record data
%     else
%       varsdata{(2*p)} = [...];		% Normal data
%     end
%   end
%   spdfcdfwrite('real',varsdata, ...
%                'GlobalAttributes', info.GlobalAttributes, ...	% Global attributes
%                'VariableAttributes', info.VariableAttributes, ... %Variable attributes
%                'RecordBound', rbvars, ...			% Var record bound
%                'varcompress',compress, ...			% Var compression 
%                'varsparse', sparse, ...			% Var sparseness
%                'blockingfactor', bf, ...			% Var blocking factors 
%                'padvalues', pad, ...				% Var Pad values
%                'cdfcompress',cdfcompress, ...			% CDF compression
%                'checksum', info.FileSettings.Checksum, ... 	% Checksum
%                'VarDatatypes', datatypes);			% Var data types
%
%   Note: The compatible data types between MATLAB and CDF are as follows:
%         MATLAB                  CDF
%         ------                --------
%         int8                  CDF_BYTE or CDF_INT1
%         int16                 CDF_INT2
%         int32                 CDF_INT4
%         int64                 CDF_INT8 or CDF_TIME_TT2000
%         uint8                 CDF_UINT1
%         uint16                CDF_UINT2
%         uint32                CDF_UINT4
%         single                CDF_FLOAT
%         double                CDF_DOUBLE or CDF_EPOCH or CDF_EPOCH16
%         char/string           CDF_UCHAR or CDF_CHAR
%
%         DEFAULT_FILLED_EPOCH_VALUE and DEFAULT_FILLED_TT2000_VALUE
%         are defined in this module so they can be used for the attribute
%         "FILLVAL" for CDF_EPOCH/CDF_EPOCH16 and CDF_TIME_TT2000 epoch
%         variables.
%
%   See also SPDFCDFREAD, SPDFCDFUPDATE, SPDFCDFINFO, CDFEPOCH, CDFTT2000,
%            SPDFENCODEEPOCH, SPDFCOMPUTEEPOCH, SPDFPARSEEPOCH,
%            SPDFBREAKDOWNEPOCH, SPDFENCODEEPOCH16, SPDFCOMPUTEEPOCH16,
%            SPDFPARSEEPOCH16, SPDFBREAKDOWNEPOCH16, SPDFENCODETT2000,
%            SPDFCOMPUTETT2000, SPDFPARSETT2000, SPDFBREAKDOWNTT2000,
%            SPDFDATENUMTOEPOCH, SPDFDATENUMTOEPOCH16, SPDFDATENUMTOTT2000,
%            SPDFCDFLEAPSECONDSINFO

% HISTORY:
%   February 12, 2009  Mike Liu     The following changes have been made to
%                                   spdfcdfwritec.c:
%                                     - Added parameter 'RecordBound'.
%                                     - Added parameter 'ConvertDatenumToEpoch'.
%                                     - Added a logic to check CDF_EPOCH and 
%                                       CDF_EPOCH16 data.

%
% Process arguments.
%

if (nargin < 2)
    error('MATLAB:spdfcdfwrite:inputArgumentCount', ...
          'SPDFCDFWRITE requires at least two input arguments.')
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
% Call the underlying spdfcdfwritec function which calls the CDFlib
%

spdfcdfwritec(filename, args.VarNames, args.VarVals, args.PadVals, ...
          globalAttribStruct, varAttribStruct, isAppending, ...
          isMultifile, CDFversion, args.ConvertDatenum, args.RecBnd, ...
          args.EpochIsCDFEpoch, args.TT2000, args.ConvertDatenum2, ...
          args.EpochTp, args.VarCompVals, args.SparseVarVals, ...
          args.BFVarVals, args.DTVarVals, args.MD5, args.CDFComp, ...
          args.CDFleapsecondlastupdated, args.SingletonVars);

%%%
%%% Function parse_inputs
%%%

function [args, isAppending, isMultifile, CDFversion, varAttribStruct, ...
          globalAttribStruct, exception] = parse_inputs(varcell, varargin)

% Set default values
args.PadVals = {};
args.RecBnd = {};
args.VarCompVals = {};
args.SparseVarVals = {};
args.BFVarVals = {};
args.EpochTp = {};
args.EpochIsCDFEpoch = false;
args.TT2000 = false;
args.ConvertDatenum = false;
args.ConvertDatenum2 = false;
args.MD5 = false;
args.CDFComp = '';
args.CDFleapsecondlastupdated = int32(-999);
isAppending = 0;
isMultifile = 0;
varAttribStruct = struct([]);
globalAttribStruct = struct([]);
% The following value indicates no version preference.
CDFversion = -1.0;

exception.msg = '';
exception.id = '';

DEFAULT_FILLED_EPOCH_VALUE = -1.0E30;
DEFAULT_FILLED_TT2000_VALUE = int64(-9223372036854775808);

% First check that varcell meets all of our requirements
args.VarNames = {varcell{1:2:end}};
args.VarVals = {varcell{2:2:end}};
% Wrap the scalars non-empties in cell arrays.
for i = 1:length(args.VarVals)
%     if ~isempty(args.VarVals{i}) && (ischar(args.VarVals{i}) || (numel(args.VarVals{i}) == 1))
     if ~isempty(args.VarVals{i}) && ischar(args.VarVals{i})
         args.VarVals{i} = {args.VarVals{i}};    
     end
end

if length(args.VarNames) ~= length(args.VarVals)
    exception.msg = 'All variable names must have a corresponding variable value.';
    exception.id = 'MATLAB:spdfcdfwrite:variableWithoutValue';
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
                    exception.id = 'MATLAB:spdfcdfwrite:inconsistentRecordTypes';
                end
            end
%        else
            % If it isn't a cell array, then it is an array and
            % all elements are of the same type. 
%            args.VarVals{i} = (args.VarVals{i});
        end
    end
end

args.PadVals = cell(1,length(args.VarNames));
args.RecBnd = cell(1,length(args.VarNames));
args.EpochTp = cell(1,length(args.VarNames));
args.VarCompVals = cell(1,length(args.VarNames));
args.SparseVarVals = cell(1,length(args.VarNames));
args.BFVarVals = cell(1,length(args.VarNames));
args.DTVarVals = cell(1,length(args.VarNames));
args.SingletonVars = cell(1,length(args.VarNames));

% Parse arguments based on their number.
if (nargin > 0)
    
    paramStrings = {'padvalues'
                    'varcompress'
                    'cdfcompress'
                    'varsparse'
                    'blockingfactor'
                    'vardatatypes'
                    'globalattributes'
                    'variableattributes'
                    'writemode'
                    'format'
                    'recordbound'
                    'convertdatenumtoepoch'
                    'convertdatenumtott2000'
                    'epochiscdfepoch'
                    'version'
                    'tt2000'
                    'checksum'
                    'cdfleapsecondlastupdated'
                    'epochtype'
                    'singleton'};
    
    % For each pair
    for k = 1:2:length(varargin)
        param = lower(varargin{k});
            
        if (~ischar(param))
            exception.msg = 'Parameter name must be a string.';
            exception.id = 'MATLAB:spdfcdfwrite:paramNameNotString';
            return
        end
        
        idx = strmatch(param, paramStrings);
        
        if (isempty(idx))
            exception.msg = sprintf('Unrecognized parameter name "%s".', param);
            exception.id = 'MATLAB:spdfcdfwrite:unrecognizedParam';
            return
        elseif (length(idx) > 1)
            exception.msg = sprintf('Ambiguous parameter name "%s".', param);
             exception.id = 'MATLAB:spdfcdfwrite:ambiguousParam';
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
               exception.id = 'MATLAB:spdfcdfwrite:paddingMismatch';
               return
           end
           vars = {padCell{1:2:end}};
           padVals = {padCell{2:2:end}};
           % Check that vars are in the list above.
           if ~iscellstr(vars)
               exception.msg = 'All variable names 1 must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(vars, args.VarNames))
               exception.msg = ['Variables listed in the PadValues ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:spdfcdfwrite:notSavingVarForPadValue';
               return
           end
           for i = 1:length(padVals)
               padVal = padVals{i};
               if (isempty(padVal)) 
                   continue
               end
               if isnumeric(padVal) || ischar(padVal) || isa(padVal,'cdfepoch') || ...
                  isa(padVal,'cdftt2000')
                   args.PadVals{strcmp(args.VarNames,vars{i})} = padVal;
               else
                   exception.msg = 'Pad values must be numbers, strings, cdfepoch or cdftt2000.';
                   exception.id = 'MATLAB:spdfcdfwrite:badPadValue';
                   return
               end
           end
        case 'varcompress'
           varCompCell = varargin{k+1};
           % If we weren't passed an even pair, then a variable
           % name or value was left out.
           if rem(length(varCompCell), 2)
               exception.msg = ['Number of variables to compress ' ...
                      'does not match number of compression values.'];
               exception.id = 'MATLAB:spdfcdfwrite:compressingMismatch';
               return
           end
	   if iscell(varCompCell{1})
	       vars = varCompCell{1:2:end};
	   else
               vars = {varCompCell{1:2:end}};
	   end
           varCompVals = {varCompCell{2:2:end}};
           % Check that vars are in the list above.
           if ~iscellstr(vars)
               exception.msg = 'All variable names 2 must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(vars, args.VarNames))
               exception.msg = ['Variables listed in the VarCompValues ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:spdfcdfwrite:notSavingVarForVarCompValue';
               return
           end
           for i = 1:length(varCompVals)
               varCompVal = varCompVals{i};
               if (isempty(varCompVal))
                   continue
               end
               if ischar(varCompVal) ...
                   args.VarCompVals{strcmp(args.VarNames,vars{i})} = varCompVal;
               else
                   exception.msg = 'Compression 3 must be strings.';
                   exception.id = 'MATLAB:spdfcdfwrite:badVarCompValue';
                   return
               end
           end
       case 'globalattributes'
           globalAttribStruct = varargin{k+1};
           if ~isstruct(globalAttribStruct)
               exception.msg = ['''GlobalAttributes''' ' must be a struct.'];
               exception.id = 'MATLAB:spdfcdfwrite:globalAttributesNotStruct';
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
               exception.id = 'MATLAB:spdfcdfwrite:variableAttributesNotStruct';
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
                exception.id = 'MATLAB:spdfcdfwrite:badWriteModeValue';
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
                exception.id = 'MATLAB:spdfcdfwrite:badFormatValue';
                return
            end
        case 'cdfcompress'
            compression = varargin{k+1};
            args.CDFComp = find_compression(compression);
        case 'version'
            version = varargin{k+1};
            if ischar(version) && ...
                    (strcmp(version,'2.7') || strcmp(version,'3.0'))
                CDFversion = str2num(version);
            else
                exception.msg = '''Version'' must be either ''2.7'' or ''3.0''';
                exception.id = 'MATLAB:spdfcdfwrite:badVersionValue';
                return
            end
        case 'recordbound'
           RecBndCell = varargin{k+1};
           % Check that vars are in the list above.
           if ~iscellstr(RecBndCell)
               exception.msg = 'All variable names for recordbound must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(RecBndCell, args.VarNames))
               exception.msg = ['Variables listed in the RECORDBOUND ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:spdfcdfwrite:notSavingVarForRecBnd';
               return
           end
           for i = 1:length(RecBndCell)
               RecBndVar = RecBndCell{i};
               args.RecBnd{strcmp(args.VarNames,RecBndCell{i})} = RecBndVar;
           end
        case 'checksum'
           if (k == length(varargin))
               msg = 'Missing "checksum" value.';
               return
           else
               checksum = varargin{k + 1};
               if ~isstr(checksum)
                 exception.msg = 'Checksum value must be a string.';
                 exception.id = 'MATLAB:spdfcdfwrite:checksumvalue';
                 return
               end
               args.MD5 = logical(find_checksum(checksum));
           end
        case 'cdfleapsecondlastupdated'
           args.CDFleapsecondlastupdated = varargin{k+1};
           if (isnumeric(args.CDFleapsecondlastupdated))
               args.CDFleapsecondlastupdated = int32(args.CDFleapsecondlastupdated);
           else
               exception.msg = 'CDF leapsecondlastupdated value must be a numeric value (in YYYYMMDD form).';
               exception.id = 'MATLAB:spdfcdfwrite:cdfleapsecondlastupdated';
               return
           end
        case 'varsparse'
           SparseVarsCell = varargin{k+1};
           % If we weren't passed an even pair, then a variable
           % name or value was left out.
           if rem(length(SparseVarsCell), 2)
               exception.msg = ['Number of variables to sparse records ' ...
                      'does not match number of their values.'];
               exception.id = 'MATLAB:spdfcdfwrite:sparserecordMismatch';
               return
           end
           vars = {SparseVarsCell{1:2:end}};
           sparseVals = {SparseVarsCell{2:2:end}};
           % Check that vars are in the list above.
           if ~iscellstr(vars)
               exception.msg = 'All variable names must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(vars, args.VarNames))
               exception.msg = ['Variables listed in the SparseValues ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:spdfcdfwrite:notSavingVarForSparseValue';
               return
           end
           for i = 1:length(sparseVals)
               sparseVal = sparseVals{i};
               if (isempty(sparseVal))
                   continue
               end
               if ischar(sparseVal)
                   args.SparseVarVals{strcmp(args.VarNames,vars{i})} = sparseVal;
               else
                   exception.msg = 'Sparse Record must be strings.';
                   exception.id = 'MATLAB:spdfcdfwrite:badSparseValue';
                   return
               end
           end
        case 'blockingfactor'
           BFVarsCell = varargin{k+1};
           % If we weren't passed an even pair, then a variable
           % name or value was left out.
           if rem(length(BFVarsCell), 2)
               exception.msg = ['Number of variables to blocking factors ' ...
                      'does not match number of their values.'];
               exception.id = 'MATLAB:spdfcdfwrite:blockingfactorMismatch';
               return
           end
           vars = {BFVarsCell{1:2:end}};
           bfVals = {BFVarsCell{2:2:end}};
           % Check that vars are in the list above.
           if ~iscellstr(vars)
               exception.msg = 'All variable names must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(vars, args.VarNames))
               exception.msg = ['Variables listed in the BFValues ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:spdfcdfwrite:notSavingVarForBFValue';
               return
           end
           for i = 1:length(bfVals)
               bfVal = bfVals{i};
               if (isempty(bfVal))
		   continue
	       end
               if isnumeric(bfVal)
                   args.BFVarVals{strcmp(args.VarNames,vars{i})} = int32(bfVal);
               else
                   exception.msg = 'Blocking factor must be numeric.';
                   exception.id = 'MATLAB:spdfcdfwrite:badBlockingFactor';
                   return
               end
           end
        case 'vardatatypes'
           DTVarsCell = varargin{k+1};
           % If we weren't passed an even pair, then a variable
           % name or value was left out.
           if rem(length(DTVarsCell), 2)
               exception.msg = ['Number of variables to blocking factors ' ...
                      'does not match number of their values.'];
               exception.id = 'MATLAB:spdfcdfwrite:blockingfactorMismatch';
               return
           end
           vars = {DTVarsCell{1:2:end}};
           dtVals = {DTVarsCell{2:2:end}};
           % Check that vars are in the list above.
           if ~iscellstr(vars)
               exception.msg = 'All variable names must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
           end
           if ~iscellstr(dtVals)
               exception.msg = 'All variable data types must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:vardatatypeNotString';
               return
           end
           if ~all(ismember(vars, args.VarNames))
               exception.msg = ['Variables listed in the DTValues ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:spdfcdfwrite:notSavingVarForDTValue';
               return
           end
           for i = 1:length(dtVals)
               dtVal = dtVals{i};
               if (~isempty(dtVal))
                   args.DTVarVals{strcmp(args.VarNames,vars{i})} = ...
                                           int32(find_datatype(dtVal));
               end
           end
        case 'epochtype'
           EpochTpCell = varargin{k+1};
           % Check that vars are in the list above.
           if ~iscellstr(EpochTpCell)
               exception.msg = 'All variable names 7 must be strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
           end
           if ~all(ismember(EpochTpCell, args.VarNames))
               exception.msg = ['Variables listed in the EPOCHTYPE ' ...
                      'cell must be on the list of variables ' ...
                      'to save.'];
               exception.id = 'MATLAB:spdfcdfwrite:notSavingVarForEpochTp';
               return
           end
           for i = 1:length(EpochTpCell)
               EpochTpVar = EpochTpCell{i};
               args.EpochTp{strcmp(args.VarNames,EpochTpCell{i})} = EpochTpVar;
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
        case 'singleton'
           if (k == length(varargin))
               msg = 'No variables specified for Singleton.';
               return
           else
             singletonVars = varargin{k+1};
             if ~iscellstr(singletonVars)
               exception.msg = 'Singleton variable names must be a cell of strings.';
               exception.id = 'MATLAB:spdfcdfwrite:varNameNotString';
               return
             end 
             % Check that vars are in the list above.
             for i = 1:length(singletonVars)
               singletonVar = singletonVars{i};
               args.SingletonVars{strcmp(args.VarNames,singletonVar)} = singletonVar;
             end
           end
        end  % switch
    end  % for
    
    % Do a sanity check on the sizes of what we are passing back
    if ~isequal(length(args.VarNames), length(args.VarVals), ...
                length(args.PadVals))
        exception.msg = 'Number of variable names, values, and pad values do not match.';
        exception.id = 'MATLAB:spdfcdfwrite:sanityCheckMismatch';
        return    
    end
    if ~isequal(length(args.VarNames), length(args.VarVals))
        exception.msg = 'Number of variable names and values do not match.';
        exception.id = 'MATLAB:spdfcdfwrite:sanityCheckMismatch';
        return    
    end
%    validate_inputs(args);

end  % if (nargin > 1)

function validate_inputs(args)
%VALIDATE_INPUTS   Ensure that the mutually exclusive options weren't provided.
%
if ((args.TT2000) && (args.EpochIsCDFEpoch))
    error('MATLAB:spdfcdfwrite:TT2000', '%s\n%s', ...
          'You cannot currently specify these two options.', ...
          'Specify only one of ''TT2000'' and ''EpochIsCDFEpoch''.')
end

if ((args.TT2000) && (args.ConvertDatenum))
    error('MATLAB:spdfcdfwrite:TT2000', '%s\n%s', ...
          'You cannot currently specify these two options.', ...
          'Specify only one of ''TT2000'' and ''ConvertDatenumToEpoch''.')
end

if ((args.ConvertDatenum) && (args.ConvertDatenum2))
    error('MATLAB:spdfcdfwrite:convertdatenum', '%s\n%s', ...
          'You cannot currently specify these two options.', ...
          'Specify only one of ''ConvertDatenumToEpoch'' and ''ConvertDatenumToTT2000''.')
end
if (((args.epochtype) && (args.EpochIsCDFEpoch)) || ...
    ((args.epochtype) && (args.TT2000)))
    error('MATLAB:spdfcdfwrite:epochtype', '%s\n%s', ...
          'You cannot currently specify these two options.', ...
          'Specify only one of ''TT2000'', ''EpochIsCDFEpoch'' and ''EpochType''.')
end

function num = find_datatype(str)
if (strcmpi(str,'double') == 1 || strcmpi(str,'cdf_double') == 1 || ...
    strcmpi(str,'cdf_real8') == 1)
  num = 45;
elseif (strcmpi(str,'single') == 1 || strcmpi(str,'cdf_float') == 1 || ...
        strcmpi(str,'cdf_real4') == 1)
  num = 44;
elseif (strcmpi(str,'int8') == 1 || strcmpi(str,'cdf_int1') == 1 || ...
        strcmpi(str,'cdf_byte') == 1)
  num = 1;
elseif (strcmpi(str,'int16') == 1 || strcmpi(str,'cdf_int2') == 1)
  num = 2;
elseif (strcmpi(str,'int32') == 1 || strcmpi(str,'cdf_int4') == 1)
  num = 4;
elseif (strcmpi(str,'int64') == 1 || strcmpi(str,'cdf_int8') == 1)
  num = 8;
elseif (strcmpi(str,'uint8') == 1 || strcmpi(str,'cdf_uint1') == 1)
  num = 11;
elseif (strcmpi(str,'uint16') == 1 || strcmpi(str,'cdf_uint2') == 1)
  num = 12;
elseif (strcmpi(str,'uint32') == 1 || strcmpi(str,'cdf_uint4') == 1)
  num = 14;
elseif (strcmpi(str,'epoch') == 1 || strcmpi(str,'cdf_epoch') == 1)
  num = 31;
elseif (strcmpi(str,'epoch16') == 1 || strcmpi(str,'cdf_epoch16') == 1)
  num = 32;
elseif (strcmpi(str,'tt2000') == 1 || strcmpi(str,'cdf_time_tt2000') == 1)
  num = 33;
elseif (strcmpi(str,'char') == 1 || strcmpi(str,'cdf_char') == 1 || ...
        strcmpi(str,'cdf_uchar') == 1)
  num = 52;
else
  error('MATLAB:spdfcdfwrite:enteredvardatatype', '%s:%s\n', ...
          'One of the entered variable data types is not valid. ', str);
end

function num = find_checksum(str)
if (strcmpi(str,'none') == 1)
  num = 0;
elseif (strcmpi(str,'md5') == 1)
  num = 1;
else
  error('MATLAB:spdfcdfwrite:checksumvalue', '%s:%s\n', ...
          'The checksum value is not valid. ', str);
end

function compress = find_compression(str)
compress = '';
lstr=lower(str);
item1=findstr(lstr, 'uncompressed');
item2=findstr(lstr, 'none');
if ~isempty(item1) || ~isempty(item2)
   compress = 'none';
else
  item1=findstr(lstr, 'gzip');
  if ~isempty(item1)
     len = length(lstr);
     if (isstrprop(lstr(len),'digit') && len < 7)
       level = str2num(lstr(len));
       if (level > 0 && level < 10) 
         if (len == 6)
           compress = lstr;
         elseif (len == 5)
           compress = strcat(lstr(1,4), ',', lstr(5));
         end
       end
     end
  else
    item1=findstr(lstr, 'run-length');
    item2=findstr(lstr, 'rle');
    if ~isempty(item1) || ~isempty(item2)
       compress = 'rle';
    else
      item1=findstr(lstr, 'adaptive');
      item2=findstr(lstr, 'ahuff');
      if ~isempty(item1) || ~isempty(item2)
        compress = 'ahuff';
      else
        item1=findstr(lstr, 'huffman');
        item2=findstr(lstr, 'huff');
        if ~isempty(item1) || ~isempty(item2)
          compress = 'huff';
        end
      end
    end
  end
end
if (length(compress) == 0)
  error('MATLAB:spdfcdfwrite:compression', '%s:%s\n', ...
          'The compression value is not valid. ', str);
end




