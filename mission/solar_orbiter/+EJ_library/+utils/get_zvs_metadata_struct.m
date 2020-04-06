%
% Convert "Variables" and "VariableAttributes" structs for CDF file into a single more easy-to-use struct suitable for
% looking up metadata for a specific zVariable.
%
%
% NOTES
% =====
% NOTE: Some Epoch zVariable attributes are represented in two different ways:
% * "VariableAttributes": VALIDMIN/-MAX, SCALEMIN/-MAX, FILLVAL on UTC format, i.e. not on the same format as
%                         the zVariable itself.
% * "Variables"         : VALIDMIN/-MAX, FILLVAL using the same data type as the zVariable iteself.
%   NOTE: The above might possibly(?) be a Solar Orbiter RPW dataset bug. This function returns both to be on the safe
%   side.
% --
% NOTE: spdfcdfinfo can return Variables in a struct format using
%   S = spdfcdfinfo(filePath, 'VARSTRUCT', true)
% It consists of a 1x1 struct where each field contains the content of one column. That is thus not a replacement for
% the struct returned by this function.
% --
% NOTE: The CDF file format permits zVariable names and zVariable attribute names that are not legal MATLAB struct
% fieldnames and may contain e.g. blanks. See "CDF USer's Guide".
% From the spdfcdfinfo help page:
% """"NOTE: Attribute names which spdfcdfinfo uses for field names in
%     "GlobalAttributes" and "VariableAttributes" may not match the names
%     of the attributes in the CDF file exactly.  Because attribute names
%     can contain characters which are illegal in MATLAB field names, they
%     may be translated into legal field names.  Illegal characters which
%     appear at the beginning of attributes are removed; other illegal
%     characters are replaced with underscores ('_').  If an attribute's
%     name is modified, the attribute's internal number is appended to the
%     end of the field name.  For example, '  Variable%Attribute ' might
%     become 'Variable_Attribute_013'.""""
% * If a zVariable           has a name that can not be used as a struct fieldname, then this function will not work (assertion).
% * If a zVariable attribute has a name that can not be used as a struct fieldname, then a modified fieldname will be
%   used, since zVariable attribute names in "VariableAttributes" are always modified that way by the "spdfcdfinfio".
%
%
% ARGUMENTS
% =========
% Variables          : Data structure on the format used by
%                      D.Variable, where D = spdfcdfinfo(....) (with 'VARSTRUCT' false i.e. the default), and
%                      D.Variable, where D = dataobj(...) .
% VariableAttributes : Data structure on the format used by
%                      D.VariableAttributes, where D = spdfcdfinfo(....), and
%                      D.VariableAttributes, where D = dataobj(...) .
%
%
% RETURN VALUE
% ============
% S : Recursive struct in three levels. Has fields
%       .(zvName)
%           .Attributes              : Struct. The reformatted content of "VariableAttributes".
%               .(zvAttributeName)
%           .Other                   : Struct. The reformatted content of "Variables".
%               .(metadataField)
%
%
% Author: Erik P G Johansson
% First created 2020-03-25
%
function S = get_zvs_metadata_struct(Variables, VariableAttributes)
    %
    % NOTE: Returned struct can be naturally generalized to incorporate content if entire CDF file (add global
    % attributes, zVariable data, and some other metadata from spdfcdfinfo(?)).
    %
    % NOTE: S = spdfcdfinfo also returns some other data, but that does not seem to contain any data from inside the CDF
    % (except "Subfiles"?).
    % >> S
    %
    %   struct with fields:
    %
    %               Filename: 'solo_L1_rpw-bia-current_20200116T164936-20200116T191610_V01_sample2020-01-29.cdf'
    %            FileModDate: '26-mar-2020 17:09:45'
    %               FileSize: 83375
    %                 Format: 'CDF'
    %          FormatVersion: '3.7.1'
    %           FileSettings: [1×1 struct]
    %               Subfiles: {}
    %              Variables: {5×12 cell}
    %       GlobalAttributes: [1×1 struct]
    %     VariableAttributes: [1×1 struct]
    %             LibVersion: '3.7.1'
    %           PatchVersion: '3.7.1.0'
    % 
    % TODO: Make execute_sw_mode (fill & pad values), validate_BIAS_master_CDFs use function.
    %
    % PROPOSAL: Make capable of handling all zVariable names.
    % PROPOSAL: Include "CDF" in filename, unless part of package "cdf".
    % PROPOSAL: Create function for reverse conversion.
    %   PRO: Useful for writing to file.
    % PROPOSAL: Assertions for agreement between "Variables" and "VariablesAttributes".
    %   Attributes & "Other": VALIDMIN, VALIDMAX, FILLVAL
    % PROPOSAL: Other name for Struct field "Other".
    %   PROPOSAL: "Metadata"    

    % ASSERTIONS: Guard against confusing the two arguments.
    assert(iscell(Variables), 'Argument "Variables" is not a cell array.')
    assert(isstruct(VariableAttributes))



    % Correct default value for the special case of having no zVars.
    S = struct();        
    
    
    
    %==================
    % Read "Variables"
    %==================
    vTable = Variables;    % Change variable name.
    for iZv = 1:size(vTable, 1)
        % NOTE: The meaning of "Variables" columns can be found on the help page for "spdfcdfinfo".
        zvName         = vTable{iZv,  1};
        sizeOfRecord   = vTable{iZv,  2};
        nRecords       = vTable{iZv,  3};
        dataType       = vTable{iZv,  4};   % NOTE: Not same as MATLAB class.
        recordVariance = vTable{iZv,  5};
        sparsity       = vTable{iZv,  6};
        compression    = vTable{iZv,  7};
        blockingFactor = vTable{iZv,  8};
        padValue       = vTable{iZv,  9};
        FILLVAL        = vTable{iZv, 10};
        VALIDMIN       = vTable{iZv, 11};
        VALIDMAX       = vTable{iZv, 12};

        % ASSERTIONS: Check values against assumptions on what they mean.
        EJ_library.assert.castring(zvName)
        assert(isnumeric(sizeOfRecord))
        EJ_library.assert.vector(sizeOfRecord)
        assert(all(sizeOfRecord >= 1))
        assert(all(round(sizeOfRecord) == sizeOfRecord))
        EJ_library.assert.castring(dataType)
        EJ_library.assert.castring(recordVariance)
        EJ_library.assert.castring(sparsity)
        EJ_library.assert.castring(compression)
        assert(isscalar(blockingFactor) && isnumeric(blockingFactor))
        assert(isscalar(padValue) || ischar(padValue))
        assert(isempty(FILLVAL)   || isscalar(FILLVAL)  || ischar(FILLVAL))
        assert(isempty(VALIDMIN)  || isscalar(VALIDMIN) || ischar(VALIDMIN))
        assert(isempty(VALIDMAX)  || isscalar(VALIDMAX) || ischar(VALIDMAX))
        
        % ~ASSERTION:
        % IMPLEMENTATION NOTE: Using try-catch to:
        %   (1) Give proper error message, instead of hard-to-understand error.
        %   (2) Make sure error happens here, rather than risking any kind of zVar mismatch when reading VariableAttributes.
        try
            % IMPLEMENTATION NOTE: Always create empty sub-struct "Attributes", just in case it is not filled with any
            % zVariable attributes later.
            S.(zvName) = struct('Other', struct(), 'Attributes', struct());
        catch Exception
            error('Argument "Variables" contains a zVariable name "%s" that can not be used as a struct fieldname.', zvName)
        end
        
        % IMPLEMENTATION NOTE: Field "zvName" added because it could be useful if sub-struct is passed around/copied
        % outside its parent struct.
        % NOTE: Not including all source variables yet. Add as the are needed.
        Other = struct();
        Other.name         = zvName;
        Other.nRecords     = nRecords;
        Other.sizeOfRecord = sizeOfRecord;
        Other.dataType     = dataType;
        Other.padValue     = padValue;
        Other.VALIDMIN     = VALIDMIN;
        Other.VALIDMAX     = VALIDMAX;
        Other.FILLVAL      = FILLVAL;
        S.(zvName).Other = Other;
    end

    
    
    % CASE: Sub-structs for all zVars have already been created.
    
    
    
    %===========================
    % Read "VariableAttributes"
    %===========================
    zvAttrNameList = fieldnames(VariableAttributes);
    for i = 1:numel(zvAttrNameList)
        zvAttrName = zvAttrNameList{i};
        
        zvAttrTable = VariableAttributes.(zvAttrName);
        
        for jZv = 1:size(zvAttrTable, 1)
            zvName      = zvAttrTable{jZv, 1};
            zvAttrValue = zvAttrTable{jZv, 2};
            
            % ASSERTION: Sub-struct for zVar already exists.
            % "VariableAttributes" shoud contain a subset of the zVars in "Variables". Do not want to just assume this,
            % but assert it so as to not mistakenly create more sub-structs.
            if ~isfield(S, zvName)
                error(...
                    'Argument "VariableAttributes" contains reference to a zVariable "%s" not present in argument "Variables".', ...
                    zvName);
            end
            S.(zvName).Attributes.(zvAttrName) = zvAttrValue;
        end
    end
    
    
    
    % Overkill?
    % ASSERTIONS
    for iZv = 1:size(vTable, 1)
        zvName = vTable{iZv,  1};
        EJ_library.assert.struct(S.(zvName), {'Attributes', 'Other'}, {})
    end
    
end
