%
% Standalone "validation" tool for BIAS master CDFs, i.e. master (template) CDF
% files for any of the datasets that BICAS produces.
%
% The code tries to verify that multiple BICAS output dataset master CDFs
% satisfy selected criteria as an alternative to manually checking them. Checks
% are thus meant to be added as they are needed. The code is primarily intended
% as a standalone tool, separate from the pipeline. It is neither intended to be
% "complete", nor to ever be constant and "finished".
%
% IMPLEMENTATION NOTE: Uses whitelists for global attributes, zVariable names,
% and zVariable attributes. This is useful for making sure there are no
% obsoleted names, misspellings, or inconsistent spelling (between datasets or
% ZVs).
%
%
% NOTE 2020-11-10: Code has not been used in a long time (years). Idea still
% relevant but code needs to be updated. Has only been used for L2, but never
% L3.
%
%
% RATIONALE
% =========
% This is useful when
% (1) the formal specification for master CDFs is modified/reinterpreted.
% (2) manually creating master CDFs based on other CDFs (e.g. from other teams,
%     other levels), or
% (3) one is particularly uncertain about something (one can then just add
% another test).
%
%
% ARGUMENTS
% =========
% dirPath        : Directory path
% filenameRegexp : Regular expression which match the entire filenames of the
%                  files to validate in the directory. The
%                  function adds ^ (beginning of string) and $ (end of string)
%                  are added automatically to the regular
%                  expression.
%
%
% RETURN VALUE
% ============
% varargout : Optional single return value. Cell array of "dataobj" objects for
%             the validated CDF files. This is useful for manually inspecting
%             the CDF files to investigate what is wrong.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2016-07-06
%
function [varargout] = validate_BIAS_master_CDFs(dirPath, filenameRegexp)
    % PROPOSAL: Validate values of VDC/EDC/EAC_LABEL in case they have been deleted by ROC bug.
    %
    % PROPOSAL: Check fill values & pad values vs. zVar data type
    %   Ex: CDF_UINT2: fill value=65535 /ISTP
    %   Ex: CDF_UINT2: pad  value=65534 /CDF User's Guide
    % 
    % PROPOSAL: Check for ISTP compliance.
    %   PROPOSAL: Check length of variable attributes.
    %       FIELD_NAM: Max 30 chars (ISTP)
    %       CATDESC:   Approx 80 chars (ISTP)
    %       LABLAXIS:  See ISTP
    % PROPOSAL: Check for SOL-SGS-TN-0009-MetadataStandard-2.4.pdf compliance
    %   QUALITY_FLAG: CDF_UINT1
    %   QUALITY_BITMASK: CDF_UINT2
    % PROPOSAL: Check that {VALID,SCALE}{MIN,MAX},FILLVAL have same data type as zVar itself.
    % PROPOSAL: QUALITY_FLAG in range 0-4 (not 0-3).
    %
    % PROPOSAL: Check zVar attributes against list(s) of mandatory (exact?) attributes?
    %   NOTE: Not sure if all have FILLVAL.
    %
    % PROPOSAL: Somehow select which checks should be made for which datasets, on which parts.
    %   (1) Checks on all ZVs in all datasets.
    %       ISTP, SolO metadata standard, ROC/RPW standards(?)
    %   (2) Checks on some ZVs (glob.attr.) in some datasets.
    %       	(2a) Checks on a specific level
    %
    % PROPOSAL: Check that no duplicate global attributes, if ignoring case.
    %   Ex: CALIBRATION_VERSION (still used 2021-09-09) + Calibration_version (obsoleted/removed)
    %
    %
    %
    % TODO-NI: Useful also for datasets (master CDFs filled with data)?
    %   NOTE: It is ultimately datasets that should be compliant and used, not
    %   master CDFs.
    %
    % PROPOSAL: Start over with a version 2 of code.
    %   PRO: Can implement checks immediately.
    %
    % PROPOSAL: Wrapper that accepts multiple paths instead of regexp. One
    %           function for validating one CDF.
    %   CON: Regexp still useful for execution from inside MATLAB (instead of
    %   bash wrapper).
    % PROPOSAL: Class with static functions.
    %   PRO: Publicly accessible functions better for testing.
    %
    % PROPOSAL: Class with constants.
    %   PRO: Might be shared among functions.
    %
    % NOTE: It appears that dataobj does not contain the CDF datatype for
    %       zVar attributes. ==> Can not (strictly) check their data types, e.g.
    %       that they match the zVar itself.

    
    
    % ASSERTION
    % IMPLEMENTATION NOTE: Natural mistake to assume that there is only one
    % argument for a single path.
    assert(nargin == 2, 'Wrong number of arguments.')

    % Needed for irf.
    % Needed for sdfpcdfread?
    % NOTE: Might run automatically run as part of user-configured MATLAB
    % initialization.
    irf('check_path');
    
    

    %=====================================
    % Iterate over all files in directory
    %=====================================
    % NOTE: irf.fs.glob_files_dirs does not add ^ and $.
    filenameRegexp = ['^', filenameRegexp, '$'];
    oiList = irf.fs.glob_files_dirs(dirPath, {filenameRegexp});   % OI = Object Info.
    oiList = oiList(~[oiList.isdir]);
    doList = {};             % DO = dataobj
    nFiles = numel(oiList);
    for i = 1:nFiles
        doList{end+1} = validate_file(oiList(i).fullPath);
    end
    if nFiles == 0
        % IMPLEMENTATION NOTE: Useful if the caller types a regex with no matches.
        % This could otherwise be mistaken for no validation warnings.
        warning('No files to validate.')
    end



    if nargout == 0
        % Do nothing
    elseif nargout == 1
        varargout = {doList};
    else
        % ASSERTION
        error('Multiple return values (more than one).')
    end
end



% 2020-11-12: Can only validate L2 datasets (yet).
%
%
% VARIABLE NAMING CONVENTION
% ==========================
% CDF_* : Refers to something "concrete" in the CDF file.
%
function Do = validate_file(filePath)
%
% PROPOSAL: Split validation into part applicable to all datasets, and part specifically for BICAS.
%    PRO: Can compare the compliance of master CDF files from other groups.
%    NOTE: Complicated(?) for semiconditional zvariable attributes for which choices have been made for BICAS datasets.
%
% PROPOSAL: Check that zVariables are correct w.r.t. sample/snapshot per record.
%    NOTE: Need list of zVariables which change with sample/snapshot per record: VDC, EDC, EAC, BIAS1/2/3.
%    PROPOSAL: Check that zVariables use the right dimension for snapshots.
%
% PROPOSAL: Check that sets of ZVs are analogous/almost identical (within a dataset).
%   Ex: IBIAS1-3
%   Ex: VDC, EDC, EAC
%   Ex: VDC/EDC/EAC_LABEL
%
% PROPOSAL: Some means of checking that large subsets of different CDF files are identical.
%    NOTE: Breaks present model of comparing CDF files independently.
%    TODO-DEC: Of a set of CDF files, compare which datasets with which?
%        PROPOSAL: Let one be special that serves as "master" for the other master files?
%    TODO-DEC: How define subsets that are be expected to be equal?
%        PROPOSAL: All attributes for a given zVariable.
%        PROPOSAL: Specific subset of global attributes.
%
% PROPOSAL: Test if FORMAT attribute covers actual skeleton-stored values in *_LABEL.
% PROPOSAL: Check that FORMAT has a correct format.



    PRINT_LATEST_GA_MODS = true;    % True/false
    % NOTE: All initials capitalized.
    RECEIVER_TEXT_MAP    = containers.Map(...
        {'LFR', 'TDS'}, ...
        {'Low Frequency Receiver', 'Time Domain Sampler'});
    MODE_TEXT_MAP        = containers.Map(...
        {'SBM1', 'SBM2', 'SURV', 'LFM'}, ...
        {'Selective Burst Mode 1', 'Selective Burst Mode 2', 'Survey Mode', 'Low Frequency Mode'});
    DATA_PRODUCT_TEXT_MAP = containers.Map(...
        {'CWF', 'SWF', 'RSWF'}, ...
        {'Continuous Waveform', 'Snapshot Waveform', 'Regular Snapshot Waveform'});



    
    
    fprintf('--------\nValidating %s\n', filePath)

    [~, fileBasename, filenameSuffix] = fileparts(filePath);
    filename = [fileBasename, filenameSuffix];

    Do         = dataobj(filePath);    % DO = dataobj
    Ga         = Do.GlobalAttributes;
    Zmd        = irf.cdf.get_zvs_metadata_struct(Do.Variables, Do.VariableAttributes);
    zvNameList = fieldnames(Zmd);
    
    
    
    if PRINT_LATEST_GA_MODS
        %=========================
        % Print latest MOD record
        %=========================
        % Useful for checking that all MODs have been updated.
        if ~isfield(Ga, 'MODS')
            fprintf('   File has no global attribute MODS. (Is it V01?)\n');
        else
            fprintf('   Last MODS = "%s"\n', Ga.MODS{end});
        end
    end
    
    
    
    %############################
    % Validate global attributes
    %############################
    
    %===================================
    % Validate set of global attributes
    %===================================
    missingDisallowedGlobAttributes = setxor(...
        bicas.tools.vbmc.const.EXACT_L2_GLOBAL_ATTRIBUTES, ...
        fieldnames(Ga));
    if ~isempty(missingDisallowedGlobAttributes)
        validation_warning_list(missingDisallowedGlobAttributes, ...
            'Found missing/disallowed global attributes:\n')
    end
    

    
    %=================================================
    % Retrieve various global attributes (later used)
    %=================================================
    CDF_Descriptor                 = Ga.Descriptor{1};       % E.g. "RPW-LFR-SURV-CWF-E>RPW Low Frequency Receiver Continuous Waveforms in survey mode. Electric component."
    CDF_DATASET_ID                 = Ga.Dataset_ID{1};       % E.g. "SOLO_L2_RPW-LFR-SURV-CWF-E"
    CDF_Skeleton_version           = Ga.Skeleton_version{1}; % E.g  "01"

    %========================================================
    % Derive various comparison values for global attributes
    %========================================================
    CDF_DATASET_ID_descriptor         = splitstr(CDF_DATASET_ID,            '_', 3, [3],     'Can not split DATASET_ID.');
    [receiver, mode, dataProduct]     = splitstr(CDF_DATASET_ID_descriptor, '-', 5, [2 3 4], 'Can not split DATASET_ID descriptor.');
    derived_CDF_DATASET_ID_descriptor = splitstr(CDF_Descriptor,            '>', 2, [1],     'Can not split Descriptor.');
    
    % Receivers     : LFR, TDS
    % Data products : CWF, SWF, RSWF
    % Modes         : SBM1, SBM2, SURV, LFM
    derived_CDF_Descriptor = sprintf('%s>RPW %s %s in %s. Electric component.', ...
        derived_CDF_DATASET_ID_descriptor, ...
        get_map_value(      RECEIVER_TEXT_MAP,     receiver), ...
        get_map_value(      DATA_PRODUCT_TEXT_MAP, dataProduct), ...
        lower(get_map_value(MODE_TEXT_MAP,         mode))   );    
    
    derived_CDF_TEXT = sprintf('This file contains RPW %s level 2 %s of electric data in %s.', ...
        receiver, ...
        lower(get_map_value(DATA_PRODUCT_TEXT_MAP, dataProduct)), ...
        lower(get_map_value(MODE_TEXT_MAP, mode)));
    
    derived_Logical_source_description = sprintf('Solar Orbiter Radio/Plasma Wave, %s L2 electric parameters', ...
        receiver);
    
    %===================================
    % Validate global attributes values
    %===================================
    validate_value(CDF_DATASET_ID_descriptor, derived_CDF_DATASET_ID_descriptor, 'DATASET_ID descriptor');
    validate_glob_attr(Ga, 'Descriptor', derived_CDF_Descriptor);
    
    % NOTE: Documentation seems to require version in SKELETON_PARENT, but in
    % reality no dataset has it.
    %derived_CDF_SKELETON_PARENT = sprintf('%s_V%s', CDF_DATASET_ID, CDF_Skeleton_version);
    derived_CDF_SKELETON_PARENT = sprintf('%s', CDF_DATASET_ID);
    validate_glob_attr(Ga, 'SKELETON_PARENT', derived_CDF_SKELETON_PARENT);
    
    validate_glob_attr(Ga, 'LEVEL', 'L2>Level 2 data processing');
    validate_glob_attr(Ga, 'TEXT',  derived_CDF_TEXT)
    
    derived_Logical_source = strrep(lower(CDF_DATASET_ID), '_l2_', '_L2_');
    validate_glob_attr(Ga, 'Logical_source',             derived_Logical_source)
    validate_glob_attr(Ga, 'Logical_source_description', derived_Logical_source_description);
    
    if ~ismember(CDF_DATASET_ID, {...
            'SOLO_L2_RPW-LFR-SBM1-CWF-E', ...
            'SOLO_L2_RPW-LFR-SBM2-CWF-E', ...
            'SOLO_L2_RPW-LFR-SURV-CWF-E', ...
            'SOLO_L2_RPW-LFR-SURV-SWF-E', ...
            'SOLO_L2_RPW-TDS-LFM-CWF-E', ...
            'SOLO_L2_RPW-TDS-LFM-RSWF-E'})
        validation_warning('Does not recognize DATASET_ID="%s"', CDF_DATASET_ID)
    end
    

    
    %#####################
    % Validate zVariables
    %#####################
    
    %============================
    % Validate set of zVariables
    %============================
    missingZvNames   = setdiff(bicas.tools.vbmc.const.MANDATORY_L2_ZV_NAMES, zvNameList);
    forbiddenZvNames = setdiff(zvNameList, ...
        union(...
            bicas.tools.vbmc.const.MANDATORY_L2_ZV_NAMES, ...
            bicas.tools.vbmc.const.SOMETIMES_L2_ZV_NAMES));
    if ~isempty(missingZvNames)
        validation_warning_list(missingZvNames, 'Can not find mandatory zVariables:\n')
    end
    if ~isempty(forbiddenZvNames)
        validation_warning_list(forbiddenZvNames, 'Found disallowed zVariables names:\n')
    end
    
    
    
    %===================================================================
    % Check Epoch
    % -----------
    % Seems values are always strings in cdfs, although
    % ROC-TST-GSE-NTT-00017-LES says it should be a number. Resolution,
    % Bin_location should probably be numbers.
    % Probably a deficiency in xlsx2skt.
    %===================================================================
    EPOCH_zVariableName = 'Epoch';   % Prescribed by the ROC-TST-GSE-NTT-00017-LES, iss2rev0.
    zVar_EPOCH_attributes = {...
        'FIELDNAM',           'Epoch'; ...
        'CATDESC',            'Default time'; ...
        'FILLVAL',            '9999-12-31T23:59:59.999999999'; ...
        'LABLAXIS',           'Epoch'; ...
        'UNITS',              'ns'; ...
        'VALIDMIN',           '2000-01-01T00:00:00.000000000'; ...
        'VALIDMAX',           '2050-12-31T23:59:59.999000000'; ...
        'SCALEMIN',           '1990-01-01T00:00:00.000000000'; ...
        'SCALEMAX',           '2050-12-31T23:59:59.999000000'; ...
        'VAR_TYPE',           'support_data'; ...
        'SCALETYP',           'linear'; ...      % Not misspelled according to ISTP/IACG Variable Attributes
        'TIME_BASE',          'J2000'; ...                    % Just taken from skeletons.
        'TIME_SCALE',         'Terrestrial Time'; ...         % Just taken from skeletons.
        'REFERENCE_POSITION', 'Spacecraft barycenter'; ...    % Just taken from skeletons.
        'Resolution',         '15258 ns'; ...                 % Just taken from skeletons.
        'Bin_location',       '0.5'; ...
        'VAR_NOTES',          'Primary time used as a reference in the file.'
    };
    for i = 1:size(zVar_EPOCH_attributes,1)
        validate_zv_attribute_value(...
            Zmd, ...
            EPOCH_zVariableName, ...
            zVar_EPOCH_attributes{i,1}, ...
            zVar_EPOCH_attributes{i,2})
    end
    

    
    %========================================
    % Validate "non-emptiness" of zVariables
    %========================================
    NONEMPTY_ZV_NAMES = {'VDC_LABEL', 'EDC_LABEL', 'EAC_LABEL'};
    for i = 1:numel(NONEMPTY_ZV_NAMES)
        validate_zv_presence(Do, NONEMPTY_ZV_NAMES{i}, 1)
    end


    
    %###############################
    % Validate zVariable attributes
    %###############################
    
    %======================================
    % Validate set of zVariable attributes
    %======================================
    % Likely obsolete comment.
    %   ROC-TST-GSESPC-00017-LES, iss01, rev 03 requires a number of zVar
    %       attributes, many of them conditionally.
    %   SRDB_PARAM_ID, SRDB_ENUM_ID - mandatory.
    %   LABLAXIS - conditional
    %   DISPLAY_TYPE - conditional
    %   FORMAT or FORM_PTR
    for iZv = 1:length(zvNameList)
        zvName         = zvNameList{iZv};
        zvAttrNameList = fieldnames(Zmd.(zvName).Attributes);
        
        missingZvAttrList   = setdiff(...
            bicas.tools.vbmc.const.MANDATORY_ZV_ATTRIBUTES, ...
            zvAttrNameList);        
        forbiddenZvAttrList = setdiff(zvAttrNameList, union(...
            bicas.tools.vbmc.const.MANDATORY_ZV_ATTRIBUTES, ...
            bicas.tools.vbmc.const.SOMETIMES_ZV_ATTRIBUTES));
        if ~isempty(missingZvAttrList)
            validation_warning_list(missingZvAttrList, 'zVariable "%s" misses attributes:\n', zvName)
        end
        if ~isempty(forbiddenZvAttrList)
            validation_warning_list(forbiddenZvAttrList, 'zVariable "%s" has forbidden attributes:\n', zvName)
        end
        
        % Check that FORMAT is uppercase.
        if isfield(Zmd.(zvName).Attributes, 'FORMAT')
            FORMAT = Zmd.(zvName).Attributes.FORMAT;
            validate_value(FORMAT, upper(FORMAT), sprintf('zVariable attribute "%s":FORMAT', zvName))
        end
    end
    
    % Validate DELTA_PLUS_MINUS:CATDESC value.
    validate_zv_attribute_value(...
        Zmd, 'DELTA_PLUS_MINUS', 'CATDESC', ...
        ['Time between sample timestamp and beginning/end of integration.', ...
        ' Total integration time is twice this value.'])

    %===========================================
    % Check FILLVAL
    % -------------
    % NOTE: Not to be confused with pad values.
    %===========================================
    ZV_FILLVAL_CHECK = {'VDC', 'EDC', 'EAC', 'IBIAS1', 'IBIAS2', 'IBIAS3'};
    ZV_FILLVAL = single(-1e31);
    for i = 1:length(ZV_FILLVAL_CHECK)
        
        zvName = ZV_FILLVAL_CHECK{i};
        validate_zv_presence(Do, zvName, 0)
        iZv = find(strcmp(Do.Variables(:,1), zvName));
        if numel(iZv) == 1
            validate_zv_attribute_value(Zmd, zvName, 'FILLVAL', ZV_FILLVAL, 0.0000)
        end
    end
    
    
    
    %==============================
    % Check filename (file system)
    %==============================
    % NOTE: Forces lower case filename suffix.    
    % NOTE: Does not compare with skeleton parent, but does indirectly. Is that appropriate?
    derivedFilename = sprintf('%s_V%s.cdf', CDF_DATASET_ID, CDF_Skeleton_version);
    if ~strcmp(filename, derivedFilename)
        validate_value(filename, derivedFilename, 'Filename (file system)')
    end
    
end



function validate_zv_attribute_value(ZvsMetadata, zvName, zvAttrName, comparisonValue, varargin)
    % TODO-DEC: How handle comparing different data types (MATLAB classes)?
    %   PROPOSAL: Act differently depending zVar data type (TT2000 or not).
    
    if ~isfield(ZvsMetadata, zvName)
        validation_warning('Can not find zVariable "%s"', zvName)
        return
    end
    if ~isfield(ZvsMetadata.(zvName).Attributes, zvAttrName)
        validation_warning('Can not find zVariable attribute "%s":"%s"', zvName, zvAttrName)
        return
    end
    
    value = ZvsMetadata.(zvName).Attributes.(zvAttrName);
    validate_value(value, comparisonValue, [zvName, ':', zvAttrName], varargin{:})
end



% ARGUMENTS
% =========
% Do               : dataobj object.
%                    NOTE: ZvsMetadata is not enough since may have to check the
%                    actual zVariable value.
% validateNonempty : Whether to validate that zVar is nonempty.
%
function validate_zv_presence(Do, zvName, validateNonempty)
    % PROPOSAL: Separate validate empty.
    %   PRO: Iterate over ZVs that should be present.
    
    if ~any(strcmp(Do.Variables(:,1), zvName))
        validation_warning('Can not find zVariable "%s".', zvName);
        return
    end

    if validateNonempty && isempty(Do.data.(zvName).data)
        validation_warning('Can not find any expected data/values for zVariable "%s".', zvName)
    end
    
    % Checks one variable attribute as a proxy for all variable attributes.
    zVarsWithAttribute = Do.VariableAttributes.FIELDNAM(:,1);
    if ~any(strcmp(zVarsWithAttribute, zvName))
        validation_warning('Can not find any FIELDNAM variable attribute for zVariable "%s".', zvName);
    end    
end



% Validate a specific value by comparing it to another value that should be
% identical.
% 
% NOTE: Special case for comparing zero.
% NOTE: Approximate comparison of numbers which can permit error even for
% permitted error zero (example: compare pad value -1e30 with -1e30). This is an
% unintended effect of a lack of precision.
%
% ARGUMENTS
% =========
% varName  : String description of value that is being validated.
% varargin : For comparing strings - No argument
%            For comparing numbers - Natural logarithm of permitted error,
%            abs(log(value/comparisonValue)).
%
function validate_value(value, comparisonValue, varName, varargin)
    msg = sprintf('%s does not match expected value.', varName);

    if ischar(value) && ischar(comparisonValue)
        % CASE: Comparing strings
        
        if numel(varargin) ~= 0
            error('Wrong number of arguments')
        end
        if ~strcmp(value, comparisonValue)
            validate_value___warning(msg, sprintf('"%s"', value), sprintf('"%s"', comparisonValue))
        end
        
    elseif isnumeric(value) && isnumeric(comparisonValue)
        %=========================
        % CASE: Comparing numbers
        %=========================
        
        if numel(varargin) ~= 1
            error('Wrong number of arguments when comparing numeric values.')
        end
        permittedError = varargin{1};
        
        if value == 0 && comparisonValue == 0
            return
        elseif value ~= 0 && comparisonValue ~= 0
            if abs(log(value/comparisonValue)) > permittedError   % Can not handle comparing zero.
                validate_value___warning(msg, sprintf('%f', value), sprintf('%f', comparisonValue))
            end
        else
            validate_value___warning(msg, sprintf('%f', value), sprintf('%f', comparisonValue))
        end
        
    else
        error('Can not compare these types of variables.')
    end
    
    
    
    %----------------------------------------------------------------------------------------------
    % Print standardized warning when comparing two values.
    % Printed line 1-N : Message string
    % Printed line N+1 : Value
    % Printed line N+2 : Comparison value
    %
    % NOTE: Does not automatically quote values so that the caller can chose
    % quoted/unquoted (e.g. strings/numbers).
    %
    function validate_value___warning(msg, value, comparisonValue)
        % NOTE: The value may be the filename, i.e. not something INSIDE the
        % file. The phrases should be chosen accordingly.
        % NOTE: Extra indentation in addition to what validation_warning adds.
        % NOTE: Example of long var_name: "ACQUISITION_TIME_UNITS(1,2)
        % (trimmed)" (37 characters)
        msg = [msg, sprintf('\n   Value            = %s',   value)];
        msg = [msg, sprintf('\n   Comparison value = %s\n', comparisonValue)];
        validation_warning(msg)
    end
    %----------------------------------------------------------------------------------------------
end



function validate_glob_attr(Ga, gaName, comparisonValue)
    
    if ~isfield(Ga, gaName)
        validation_warning('There is no global attribute "%s".', gaName)
        return
    end
    gaValue = Ga.(gaName);
    
    assert(isscalar(gaValue), 'Can not handle non-scalar global attribute value.')
    
    gaValue = gaValue{1};
    validate_value(gaValue, comparisonValue, sprintf('Global attribute "%s"', gaName))
end



function validation_warning_list(strList, varargin)
    titleMsg = sprintf(varargin{:});
    
    assert(titleMsg(end) == newline)
    
    % NOTE: Submitting one string without formatting.
    validation_warning([titleMsg, '   "', strjoin(strList, ['"', newline, '   "']), newline])
end



% Print warning on standardized format.
%
% varargin : Arguments to sprintf.
%            Can handle multiline messages.
%
function validation_warning(varargin)
    msg    = sprintf(varargin{:});
    msg    = irf.str.indent(msg, 3);
    msg(2) ='*';

    fwrite(1, msg)
end



% Get key in Map with default value rather than error if key does not exist.
function value = get_map_value(map, key)
    if map.isKey(key)
        value = map(key);
    else
        % Best to return something so that something is inserted into strings.
        % Absence can be ambiguous, e.g. when deriving Descriptor or TEXT.
        % NOTE: The caller might change the case.
        value = '<Unknown>';   
    end
end



% Split string using delimiter into substrings and return a specified subset of
% these. Verifies number of substrings. Useful for analyzing e.g. dataset IDs.
%
% ARGUMENTS
% =========
% str       : String to be split.
% delimiter : Delimiter used for splitting.
% nParts    : The expected number of substrings (error if wrong).
% iParts    : The substrings to be returned (numeric array).
%
function varargout = splitstr(str, delimiter, nParts, iParts, valWarningMsg)
    
    if nargout ~= length(iParts)
        error('Number of return values does not match input. This indicates a pure bug.')
    end

    parts = strsplit(str, delimiter);   % Returns ROW vector.
    if length(parts) ~= nParts
        validation_warning([valWarningMsg, '. str="', str, '"'])
        parts = [parts, cell(1, nParts-length(parts))];   % Extend ROW vector.
    end
    
    for iOut = 1:length(iParts)
        varargout{iOut} = parts{iParts(iOut)};
    end
end
