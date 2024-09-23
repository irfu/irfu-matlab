%
% Class which represents a parsable dataset filename on an ~official filenaming
% convention. Contains code for converting between (a) filename and (b)
% separate filename fields.
%
%
% RATIONALE
% =========
% Useful for
% * identifying and grouping datasets (based on filenames) for input to
%   BICAS during e.g. IRFU-local batch processing
% * creating filenames to be output from BICAS during e.g. IRFU-local batch
%   processing.
% * identifying datasets when comparing file modification dates for datasets
%   and QLI files.
% --
% Primarily meant for making it easy to
% * Distinguish between datasets and non-datasets.
% * Centralize the creation and parsing of filenames.
% * Centralize the understanding of filenaming conventions (as far as relevant).
% --
% Implementation as a class:
% * Makes it easy to support returning redundant values used for e.g. CDF GAs.
%   Ex: timeIntervalStr, versionStr
% * Avoids storing fields in a struct.
%
%
% NOTES
% =====
% * Meant to cover DSIs widely.
%   * Old DSIs beginning with ROC-SGSE, in order to cover old filenames.
%   * Datasets unrelated to BIAS processing, including non-RPW instruments.
%   * RCTs / CAL level.
%     NOTE: RCTs count as datasets in this context since they conform to the
%     same filenaming convention and are described in SOL-SGS-TN-0009, "Metadata
%     Definition for Solar Orbiter Science Data", 01/06, though using its own
%     level "CAL".
% * Recognizes -cdag (RPW-internal convention).
% * NOTE/BUG: Can not handle filenames where the DSI has mixed case in the
%             descriptor (excluding level).
%             Ex: solo_L1_swa-eas2-NM3D_20201027T000007-20201027T030817_V01.cdf
%
%
% VARIABLE NAMING CONVENTION
% ==========================
% DSICDAG = String containing the combination of DATASET_ID (DSI) and
%           (optionally) the suffix "-CDAG".
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef DatasetFilename
  % PROPOSAL: datasetId --> dsi
  % PROPOSAL: Require NaT when there is no time interval string.
  % PROPOSAL: Order the filenaming conventions.
  % PROPOSAL: Official abbreviation: DSFN = solo.adm.dsfn.DatasetFilename.
  %
  % PROBLEM: Risk of duplicating assertions on consistent fields. Which code
  %          should be responsible? Constructor? create/parse_dataset_filename?
  %
  % PROPOSAL: Field (string constant) for type of filenaming convention.
  %   NOTE: Only if multiple filenaming conventions in the same code.
  %
  % PROPOSAL: Property for file basename.
  %   PRO: Used by bicas.ga.get_output_dataset_GAs().
  %
  % PROPOSAL: Use field names identical to the terms used in specifications (RCS
  %           ICD, SOL-SGS-TN-0009).
  %   Ex: Datetime, descriptor
  %   CON: Terms do not follow own variable naming conventions.
  %
  % PROPOSAL: Forbid >2 char-length version strings beginning with zero.
  %
  %
  %
  % PROPOSAL: Separate classes for every separate filenaming scheme.
  %   CON: Harder to reuse similarities.
  %       Ex: DATASET_ID, version, time interval string
  %       CON-PROPOSAL: Use shared RE constants, functions.
  %   CON: Can not build composite function for simultaneously handling all
  %        filenaming conventions (as the current create/parse functions do).
  %     CON: Filenaming conventions are too dissimilar for doing that anyway: Have
  %          different variables.
  %       Ex: "CNE" filenaming convention has hash, has no time interval string.
  %       Ex: "LES" filenaming convention has hash.
  %       Ex: Inflight filenaming has isCdag (has no hash).
  %   CON: solo.adm.dsfn.parse_dataset_filename_many() would only apply to one
  %        filenaming convention.
  %     CON-PROPOSAL: Can call functions for multiple filenaming conventions.
  %       PRO: Functions would be made for making this easy.
  %     NOTE: solo.adm.dsfn.parse_dataset_filename_many() is only used by
  %           solo.adm.paths_to_DSMD_array().
  %     CON-PROPOSAL: Merge solo.adm.dsfn.parse_dataset_filename_many() into
  %                   solo.adm.paths_to_DSMD_array().
  %   PRO: Cleaner functions, classes (sets of properties).
  %       PRO: Consistent return format.
  %           CON: Not true for time vectors which reflect format of time
  %                interval string.
  %       PRO: Can easily have different requirements for different filenaming
  %            conventions.
  %           Ex: isCdag (only for official; not LES and CNES filenaming).
  %   PROPOSAL: Shared code for shared implementations.
  %       Ex: utility functions (used by top-level functions)
  %       Ex: Reg.exp. constants
  %   PROPOSAL: One top-level function that handles all filenaming conventions
  %             together. Should assign(?) string for identified convention.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % Dataset ID. NOTE: Always uppercase.
    datasetId

    % Logical. Whether or not the file is a CDAG
    % (whether dataset ID in filename is appended with "-CDAG"/"-cdag").
    isCdag

    % Start time (datetime).
    Dt1

    % End time (datetime).
    Dt2

    % String constant specifying the time interval format.
    timeIntervalFormat

    % Version number.
    versionNbr

    lesTestStr
    cneTestStr
  end

  properties(GetAccess=private, SetAccess=immutable)
    dsicdagUppercase
  end

  properties(Dependent)
    % The combination of DSI and optionally "-cdag" as found in
    % the filename, including case.
    % In practice meant to be interpreted as dataset glob.attr.
    % "Logical_source", which should include -CDAG when present (for
    % now at least).
    filenameDsiCdag

    timeIntervalStr

    % Version number as a string the way it is represented in
    % the filename. Excludes "V".
    versionStr

    filename
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % ARGUMENTS
    % =========
    % S
    %       Struct which de facto serves as a way to supply named arguments.
    function obj = DatasetFilename(S)
      assert(isstruct(S))

      %================================
      % Check legality and consistency
      %================================
      assert(ischar(S.datasetId))
      assert(islogical(S.isCdag))
      irf.dt.assert_UTC(S.Dt1)
      irf.dt.assert_UTC(S.Dt2)

      %=========================================================================
      % Lowest legal version number
      % ---------------------------
      % The lowest dataset version number should FORMALLY be 1.
      % """"The value of the version shall be consistent with the "version"
      % field in the filename and shall be an integer increment starting at
      % 01.""""
      % Source: SOL-SGS-TN-0009, 02/06, Section 3.2.2.1 "ISTP Global Attributes"
      % Table 3-16 "ISTP Global Attributes for Solar Orbiter Files", entry for
      % "Data_version".
      %
      % CAL files should have lowest version 1.
      % """"Initial version number shall be "01".""""
      % Source: ROC-PRO-PIP-ICD-00037-LES, 01/07, Section 4.2.4
      % "Data versioning",
      %
      % HOWEVER, there are datasets on SOAR with "V00" in the filename and these
      % are mirrored at IRFU.
      % Ex: solo_L1_swa-his-pha_20200310_V00.cdf
      %     solo_L1_swa-his-sensorrates_20200310_V00.cdf
      % NOTE: These files appear to have GA Data_version = "00" too.
      %=========================================================================
      assert(isnumeric(S.versionNbr) && (S.versionNbr >= 0))
      assert(ischar(S.timeIntervalFormat))
      assert(isempty(S.lesTestStr) || ischar(S.lesTestStr))
      assert(isempty(S.cneTestStr) || ischar(S.cneTestStr))
      % NOTE: Not asserting timeIntervalFormat since that is indirectly
      % checked below.

      obj.datasetId          = S.datasetId;
      obj.isCdag             = S.isCdag;
      obj.Dt1                = S.Dt1;
      obj.Dt2                = S.Dt2;
      obj.timeIntervalFormat = S.timeIntervalFormat;
      obj.versionNbr         = S.versionNbr;

      obj.lesTestStr         = S.lesTestStr;
      obj.cneTestStr         = S.cneTestStr;

      %===========================================
      % Check consistency between variables, etc.
      %===========================================
      if ~isempty(S.lesTestStr) && isempty(S.cneTestStr) && ~obj.isCdag ...
          && ismember(obj.timeIntervalFormat, {'DAY_TO_DAY', 'SECOND_TO_SECOND'})
        %=============================
        % "LES" filenaming convention
        %=============================
        obj.dsicdagUppercase = false;

      elseif isempty(S.lesTestStr) && ~isempty(S.cneTestStr) && ~obj.isCdag ...
          && strcmp(obj.timeIntervalFormat, {'NO_TIME_INTERVAL'})
        %=============================
        % "CNE" filenaming convention
        %=============================
        obj.dsicdagUppercase = true;

      elseif isempty(S.lesTestStr) && isempty(S.cneTestStr) ...
          && ismember(obj.timeIntervalFormat, {'DAY', 'DAY_TO_DAY', 'SECOND_TO_SECOND'})
        %=================================
        % In-flight filenaming convention
        %=================================
        obj.dsicdagUppercase = false;

      else

        error('Inconsistent struct argument')
      end

    end



  end    % methods(Access=public)



  %######################################
  %######################################
  % PUBLIC INSTANCE DEPENDENT PROPERTIES
  %######################################
  %######################################
  methods



    % dsicdagStr
    %       DSI+optional CDAG, but with the case to be used in the filename.
    function filenameDsiCdag = get.filenameDsiCdag(obj)
      if obj.isCdag
        cdagStr = '-CDAG';
      else
        cdagStr = '';
      end

      if obj.dsicdagUppercase
        cdagStr        = upper(cdagStr);
        filenameDsiStr = upper(obj.datasetId);
      else
        cdagStr        = lower(cdagStr);

        % NOTE: The DSI "level" is always uppercase. Can therefore not user
        % lower(datasetId).
        [sourceName, level, descriptor] = solo.adm.disassemble_DATASET_ID(obj.datasetId);
        filenameDsiStr = sprintf('%s_%s_%s', lower(sourceName), level, lower(descriptor));
      end

      filenameDsiCdag = [filenameDsiStr, cdagStr];
    end



    function timeIntervalStr = get.timeIntervalStr(obj)
      timeIntervalStr = solo.adm.dsfn.create_time_interval_str(...
        obj.Dt1, obj.Dt2, obj.timeIntervalFormat);
    end



    function versionStr = get.versionStr(obj)
      versionStr = sprintf('%02i', obj.versionNbr);
    end



    function filename = get.filename(obj)

      if ~isempty(obj.lesTestStr) ...
          && ismember(obj.timeIntervalFormat, {'DAY_TO_DAY', 'SECOND_TO_SECOND'})
        %=============================
        % "LES" filenaming convention
        %=============================
        filename = sprintf('%s_%s_V%s_%s.cdf', ...
          obj.filenameDsiCdag, obj.timeIntervalStr, obj.versionStr, obj.lesTestStr);

      elseif ~isempty(obj.cneTestStr) ...
          && strcmp(obj.timeIntervalFormat, {'NO_TIME_INTERVAL'})
        %=============================
        % "CNE" filenaming convention
        %=============================
        filename = sprintf('%s_%s_V%s.cdf', ...
          obj.filenameDsiCdag, obj.cneTestStr, obj.versionStr);

      elseif ismember(obj.timeIntervalFormat, {'DAY', 'DAY_TO_DAY', 'SECOND_TO_SECOND'})
        %===========================================
        % In-space filename convention (incl. CDAG)
        %===========================================
        filename   = sprintf('%s_%s_V%s.cdf', ...
          obj.filenameDsiCdag, obj.timeIntervalStr, obj.versionStr);

      else

        error('Illegal state.')
      end

    end



  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function Df = parse_filename(filename)

      function versionStr = version_RE_match_to_versionNbr(s)
        versionStr = str2double(s(2:end));
      end

      NO_MATCH_RETURN_VALUE = [];
      UNUSED_DT             = datetime('NaT', 'TimeZone', 'UTCLeapSeconds');



      % NOTE: Parse from the END.
      [~, trueBasename, n] = irf.str.read_token(filename, -1, '\.cdf');
      if n == -1
        Df = NO_MATCH_RETURN_VALUE;
        return
      end



      %=========================
      % Parse DATASET_ID + CDAG
      %=========================
      % IMPLEMENTATION NOTE: Could use separate "-CDAG" regexp but then have to
      % use negative lookahead to exclude "-CDAG" from the preceding DATASET_ID
      % (due to maximal munch). Could then use irf.str.read_token()
      % again for "-CDAG", but would then have to have separate cases of lowercase
      % and uppercase anyway. Might also be slower because of negative lookahead.
      %
      % NOTE: Lowercase DATASET_ID+CDAG always have uppercase dataset level.
      [filenameDsicdag, str, n] = irf.str.read_token(trueBasename, 1, ...
        '(SOLO|ROC-SGSE)_(HK|L1|L1R|L2|L3|CAL)_[A-Z0-2-]*', ...
        '(solo|roc-sgse)_(HK|L1|L1R|L2|L3|CAL)_[a-z0-2-]*');
      switch(n)
        case 1
          dsicdagUppercase = true;
          S.isCdag         = strcmp(filenameDsicdag(end-4:end), '-CDAG');
        case 2
          dsicdagUppercase = false;
          S.isCdag         = strcmp(filenameDsicdag(end-4:end), '-cdag');
        otherwise
          Df = NO_MATCH_RETURN_VALUE;
          return
      end
      if S.isCdag
        S.datasetId = upper(filenameDsicdag(1:end-5));
      else
        S.datasetId = upper(filenameDsicdag);
      end



      % Recurring regular expressions.
      VERSION_RE     = 'V[0-9][0-9]+';    % NOTE: Permits more than two digits.
      LES_TESTSTR_RE = 'les-[0-9a-f]{7,7}';
      CNE_TESTSTR_RE = '[0-9a-f]{7,7}_CNE';



      %===================================================================
      % NOTE: Checks for different naming conventions in assumed order of
      % decreasing likelyhood of being used.
      % ==> Potential speedup.
      %===================================================================

      % Non-rigorous reg.exp. for time interval string. A more rigorous check
      % is made in solo.adm.dsfn.parse_time_interval_str().
      TIME_INTERVAL_STR_RE = '[0-9T-]{8,31}';

      %==========================================
      % Standard in-flight filenaming convention
      % (incl. CDAG; tested for earlier)
      %==========================================
      [subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
        {'_', TIME_INTERVAL_STR_RE, '_', VERSION_RE}, ...
        'permit non-match');
      if perfectMatch && ~dsicdagUppercase

        [S.Dt1, S.Dt2, S.timeIntervalFormat] = ...
          solo.adm.dsfn.parse_time_interval_str(subStrCa{2});
        S.versionNbr      = version_RE_match_to_versionNbr(subStrCa{4});
        S.lesTestStr      = [];
        S.cneTestStr      = [];

        Df = solo.adm.dsfn.DatasetFilename(S);
        return
      end

      %=============================
      % "LES" filenaming convention
      %=============================
      [subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
        {'_', TIME_INTERVAL_STR_RE, '_', VERSION_RE, '_', LES_TESTSTR_RE}, ...
        'permit non-match');
      if perfectMatch && ~dsicdagUppercase && ~S.isCdag

        [S.Dt1, S.Dt2, S.timeIntervalFormat] = ...
          solo.adm.dsfn.parse_time_interval_str(subStrCa{2});
        S.versionNbr      = version_RE_match_to_versionNbr(subStrCa{4});
        S.lesTestStr      = subStrCa{6};
        S.cneTestStr      = [];

        Df = solo.adm.dsfn.DatasetFilename(S);
        return
      end

      %=============================
      % "CNE" filenaming convention
      %=============================
      [subStrCa, ~, perfectMatch] = irf.str.regexp_str_parts(str, ...
        {'_', CNE_TESTSTR_RE, '_', VERSION_RE}, ...
        'permit non-match');
      if perfectMatch && dsicdagUppercase && ~S.isCdag

        % S.timeIntervalStr    = NO_MATCH_RETURN_VALUE;

        S.Dt1                = UNUSED_DT;
        S.Dt2                = UNUSED_DT;
        S.timeIntervalFormat = 'NO_TIME_INTERVAL';
        %
        S.lesTestStr         = [];
        S.cneTestStr         = subStrCa{2};
        S.versionNbr         = version_RE_match_to_versionNbr(subStrCa{4});

        Df = solo.adm.dsfn.DatasetFilename(S);
        return
      end

      %============================================
      % CASE: Could not match filename to anything
      %============================================
      Df = NO_MATCH_RETURN_VALUE;

    end    % parse_filename()



  end    % methods(Static)



end
