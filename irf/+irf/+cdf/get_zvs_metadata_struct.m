%
% Convert "Variables" and "VariableAttributes" structs for CDF file into a
% single more easy-to-use struct suitable for looking up metadata for a specific
% zVariable.
%
%
% NOTES
% =====
% NOTE: Some zVariable attributes are represented in two different ways:
% VALIDMIN, VALIDMAX, FILLVAL. Empirically, the following is true.
%
% * "VariableAttributes":
%   * Epoch VALIDMIN/-MAX, FILLVAL are on UTC format (string), i.e. not on the
%     same format as the zVariable itself.
%   * Regular numeric VALIDMIN/-MAX, FILLVAL use the MATLAB class corresponding
%     to the zVar itself.
%
% * "Variables"
%   * VALIDMIN/-MAX, FILLVAL using scalars (double) the same data type as
%     the zVariable iteself.
%   * Regular numeric VALIDMIN/-MAX, FILLVAL always(?) use double, regardless of
%     the MATLAB class corresponding to the zVar itself.
%   * NOTE: padValue (only in "Variables" does seem to have the correct MATLAB
%     class.
%
% --
% NOTE: spdfcdfinfo can return Variables in a struct format using
%   S = spdfcdfinfo(filePath, 'VARSTRUCT', true)
% It consists of a 1x1 struct where each field contains the content of one
% column. This is thus not a replacement for the struct returned by this
% function.
% --
% NOTE: The CDF file format permits zVariable names and zVariable attribute
% names that are not legal MATLAB struct fieldnames and may contain e.g. blanks.
% See "CDF USer's Guide".
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
% * If a ZVARIABLE           has a name that can not be used as a struct
%   fieldname, then this function will not work (assertion).
% * If a ZVARIABLE ATTRIBUTE has a name that can not be used as a struct
%   fieldname, then a modified fieldname will be used, since zVariable attribute
%   names in "VariableAttributes" are always modified that way by the
%   "spdfcdfinfio".
%
% DESIGN INTENT
% =============
% Not to be dependent on the implementation of irfu-matlab's dataobj so that it
% can be used with dataobj and without it, e.g. if designing an alternative to
% dataobj.
%
%
% ARGUMENTS
% =========
% Variables          : Data structure on the format used by
%                      D.Variable, where D = spdfcdfinfo(....) (with
%                      'VARSTRUCT' false i.e. the default), and
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
%           .Attributes       : Struct. The reformatted content of "VariableAttributes".
%               .(zvAttributeName)
%           .Other            : Struct. The reformatted content of "Variables".
%               .(metadataField)
%               NOTE: .FILLVAL, .VALIDMIN, .VALIDMAX (under .Other) are empty if
%               CDF does not contain values.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-03-25
%
function S = get_zvs_metadata_struct(Variables, VariableAttributes)
    %
    % NOTE: Returned struct can be naturally generalized to incorporate content
    % if entire CDF file (add global attributes, zVariable data, and some other
    % metadata from spdfcdfinfo(?)).
    %
    % NOTE: S = spdfcdfinfo also returns some other data, but that does not seem
    % to contain any data from inside the CDF (except "Subfiles"?).
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
    %
    % PROPOSAL: Better function name.
    %   PROPOSAL: get_zvs_metadata (no "struct").
    %       CON: Does not imply better data structure that default.
    %       CON: Does not imply conversion from one format to another.
    %   PROPOSAL: reformat_zvs_metadata
    %   PROPOSAL: convert_zvs_metadata
    %
    % PROPOSAL: Convert into function for reading entire CDF.
    %   PROPOSAL: Add field S.(zVarName).value    or .data .
    %   PROPOSAL: Support (option) for converting fill value and pad value to NaN.
    %       NOTE: Fill value give in two different places which might differ (for tt2000).
    %       PROBLEM/TODO-NI: How does fill value work for strings/chars?
    %       PROBLEM: Can not be done for non-floats. ==> Must convert data type.
    %           Ex: Epoch
    %           Ex: IBIAS_1
    %           PROBLEM: How specify when to convert?
    %               PROPOSAL: Convert all (numeric) zVars to double.
    %                   ~CON: Epoch?
    %               PROPOSAL: Specify list of zVar names.
    %               PROPOSAL: Use class (for entire CDF content). Use method(s) for reading. Can supply options
    %                   (1) raw (with correct type/class)
    %                   (2) double with fill value, pad value replaced with NaN.
    %                   PRO: Can implement more support over time (assert when supports).
    %
    % PROPOSAL: Class for the content of a CDF file:
    %   (1) content of a loaded CDF file
    %   (2) content of a CDF file to be created
    %   (3) partly correct content of CDF file to be created (not necessarily
    %       corrupt/inconsistent data, just different from any pre-existing or future CDF file)
    %   PRO: Easy to load, modify, and write CDF file.
    %   PRO: zVar attributes on easy-to use format.
    %   PRO: Can implement assertions on content directly when setting, loading.
    %       Ex: DEPENDS_0=Epoch ==> Must have same nRecords.
    %       Ex: zVar attributes that should point to other zVars actually do.
    %   PRO: Can enforce that empty zVar values have the correct size (except char?)
    %   --
    %   PROPOSAL: Methods for adding/removing zVars
    %   PROPOSAL: Methods for reading/writing zVar value
    %       (1) raw content: exact CDF type/MATLAB class disregarding/ignoring pad & fill value
    %       (2) converted content: replace NaN for fill value (pad value?), convert double-->CDF type
    %           Check feasability: min-max, Inf, NaN.
    %   TODO-DEC: How reference zvars? (since not struct of structs)
    %       PROPOSAL: zVarName as method argument
    %           CON: Less convenient than structs.
    %               PROPOSAL: Method for returning simples struct of fields with zVar values.
    %   --
    %   TODO-NI: How easy is it to read a CDF file NOT using dataobj?!!
    %
    % PROPOSAL: Only select and return one version of duplicated variable
    % attributes (from either "Variables" or "VariableAttributes").
    %   PROPOSAL: Assertions on values not used to check the assumptions the
    %   code makes, e.g. that TT2000 always expects VALIDMIN/-MAX & FILLVAL
    %   begin UTC and should use "Variables".
    %   PROPOSAL: Convert TT2000 zVar char string to TT2000 (int64).



    % ASSERTIONS: Guard against confusing the two arguments.
    assert(iscell(Variables), 'Argument "Variables" is not a cell array.')
    assert(isstruct(VariableAttributes) && isscalar(VariableAttributes))
    nZvars = irf.assert.sizes(Variables, [-1, 12]);



    % Correct default value for the special case of having no zVars.
    S = struct();        
    
    
    
    %==================
    % Read "Variables"
    %==================
    % Change variable name, since original name is deceiving and comes from
    % spdfcdfinfo's a naming convention.
    vTable = Variables;
    clear Variables
    for iZv = 1:nZvars
        % NOTE: The meaning of "Variables" columns can be found on the help page
        % for "spdfcdfinfo".
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
        irf.assert.castring(zvName)
        assert(isnumeric(sizeOfRecord))
        irf.assert.vector(sizeOfRecord)
        assert(all(sizeOfRecord >= 1))
        assert(all(round(sizeOfRecord) == sizeOfRecord))
        irf.assert.castring(dataType)
        irf.assert.castring(recordVariance)
        irf.assert.castring(sparsity)
        irf.assert.castring(compression)
        assert(isscalar(blockingFactor) && isnumeric(blockingFactor))
        assert(isscalar(padValue) || ischar(padValue))
        assert(isempty(FILLVAL)   || isscalar(FILLVAL)  || ischar(FILLVAL))
        assert(isempty(VALIDMIN)  || isscalar(VALIDMIN) || ischar(VALIDMIN))
        assert(isempty(VALIDMAX)  || isscalar(VALIDMAX) || ischar(VALIDMAX))
        
        % ~ASSERTION:
        % IMPLEMENTATION NOTE: Using try-catch to:
        %   (1) Give proper error message, instead of hard-to-understand error.
        %   (2) Make sure error happens here, rather than risking any kind of
        %       zVar mismatch when reading VariableAttributes.
        try
            % IMPLEMENTATION NOTE: Always create empty sub-struct "Attributes",
            % just in case it is not filled with any zVariable attributes later.
            S.(zvName) = struct('Other', struct(), 'Attributes', struct());
        catch Exception
            error(...
                ['Argument "Variables" contains a zVariable name "%s" that', ...
                ' can not be used as a struct fieldname.'], zvName)
        end
        
        % IMPLEMENTATION NOTE: Field "zvName" added because it could be useful
        % if sub-struct is passed around/copied outside its parent struct.
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
            % "VariableAttributes" should contain a subset of the zVars in
            % "Variables". Do not want to just assume this, but assert it so as
            % to not mistakenly create more sub-structs.
            if ~isfield(S, zvName)
                error(...
                    ['Argument "VariableAttributes" contains reference to a', ...
                    ' zVariable "%s" not present in argument "Variables".'], ...
                    zvName);
            end
            S.(zvName).Attributes.(zvAttrName) = zvAttrValue;
        end
        
    end
    
    
    
    % Overkill?
    % ASSERTIONS
    for iZv = 1:nZvars
        zvName = vTable{iZv,  1};
        irf.assert.struct(S.(zvName), {'Attributes', 'Other'}, {})
        
        Zva = S.(zvName).Attributes;
        
        ZVAR_ATTR_CA = {'FILLVAL', 'VALIDMIN', 'VALIDMAX'};
        for i = 1:numel(ZVAR_ATTR_CA)
            zvAttrName = ZVAR_ATTR_CA{i};
            if isfield(Zva, zvAttrName)
                zvAttrValue = Zva.(zvAttrName);
                
                % ASSERTION
                % Could be a bad assertion. Copied from reading "Variables".
                %assert(isempty(zvAttrValue) || isscalar(zvAttrValue)  || ischar(zvAttrValue))
                assert(isscalar(zvAttrValue)  || ischar(zvAttrValue))
            end
        end
        
    end
    
end
