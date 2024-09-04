%
% Functions (static methods) associated with bicas.proc.L1L2.cal.Cal finding,
% reading, and logging RCTs so that bicas.proc.L1L2.cal.Cal does not need to.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2021-08-16.
%
classdef findread
  % PROPOSAL: Abbreviation for "bias_rct_validity.json".
  %   BIAS, RCT, validity, JSON
  %   file
  %   BRVJ = BIAS RCT Validity JSON (also the words in the filename)
  %   BRJ  = BIAS RCT JSON
  %   BRJF = BIAS RCT JSON File
  %   BRVF = BIAS RCT Validity File



  %#####################
  %#####################
  % PRIVATE PROPERTIES
  %#####################
  %#####################
  properties(Access=private, Constant)

    READING_RCT_PATH_LL = 'info';
  end



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    function RctdCaMap = get_RctdCaMap(...
        useGactRct, nonBiasRcttid, rctDir, gact, zvcti, zv_BW, L)

      if useGactRct
        RctdCaMap = bicas.proc.L1L2.cal.rct.findread.find_read_RCTs_by_RCT_JSON_and_ZVCTI_GACT(...
          nonBiasRcttid, rctDir, ...
          gact, zvcti, zv_BW, ...
          L);
      else
        RctdCaMap = bicas.proc.L1L2.cal.rct.findread.find_read_RCTs_by_RCT_JSON_and_regexp(...
          nonBiasRcttid, rctDir, Bso, L);
      end
    end



    % Read the BIAS RCT JSON file and return the content.
    function [biasRctPath, DtValidityBegin, DtValidityEnd] = ...
        read_BIAS_RCT_JSON(rctDir, L)

      rctJsonPath = fullfile(rctDir, bicas.const.RCT_JSON_FILENAME);
      irf.assert.file_exists(rctJsonPath)

      L.logf('info', 'Reading "%s".', rctJsonPath)
      jsonStr    = fileread(rctJsonPath);
      JsonStruct = jsondecode(jsonStr);

      fnCa = fieldnames(JsonStruct);
      assert(length(fnCa) == 1, ...
        'File "%s" does not reference exactly one BIAS RCT as expected.', ...
        rctJsonPath)

      rctFilename   = fnCa{1};
      % NOTE: Renaming "start"-->"begin" since "begin" is more conventional.
      validityBegin = JsonStruct.(fnCa{1}).validity_start;
      validityEnd   = JsonStruct.(fnCa{1}).validity_end;

      %=====================================================
      % Correct the RCT filename returned from jsondecode()
      %=====================================================
      % IMPLEMENTATION NOTE: jsondecode() stores the filename in a struct field
      % name, but struct field names do not permit all characters (such as dash
      % and period) and replaces them with underscore instead. Must therefore
      % correct the filename string using knowledge of legal RCT filenames
      % (sigh...).
      %
      % Ex: solo_CAL_rpw-bias_20200210-20991231_V01.cdf
      %                 ^             ^            ^
      rctFilename(end-3) = '.';
      %
      iDsiDash = strfind(bicas.const.RCT_DSI, '-');
      assert(isscalar(iDsiDash))
      rctFilename(iDsiDash) = '-';
      %
      assert(strcmp(rctFilename(27), '_'))
      rctFilename(27) = '-';

      % Construct return values.
      biasRctPath = fullfile(rctDir, rctFilename);
      DtValidityBegin = datetime(validityBegin, 'TimeZone', 'UTCLeapSeconds');
      DtValidityEnd   = datetime(validityEnd,   'TimeZone', 'UTCLeapSeconds');
    end



    % (1) Load one BIAS RCT using the RCT JSON, and
    % (2) load one non-BIAS RCT by only using assumptions on filenames (and
    %     directory).
    %
    %
    % NOTES
    % =====
    % NOTE: Can be useful for manual experimentation with calibration of L1R
    %       (and L1) data.
    % NOTE: Necessary when processing L1-->L2 (unofficially) since L1 does
    %       not have ZVCTI+GACT.
    % NOTE: Will only load ONE RCT per RCTTID (since there is no potential RCT
    %       time dependence as there would be if using ZVCTI+GACT) and
    %       requires the caller to not use ZVCTI.
    %
    % RETURN VALUE
    % ============
    % RctdCaMap
    %       containers.Map. Can be used for bicas.proc.L1L2.cal.Cal
    %       constructor even if there is no ZVCTI+GACT.
    %       One key per specified RCTTID in argument rcttidCa.
    %       Exactly one RCT per RCTTID.
    %
    function RctdCaMap = find_read_RCTs_by_RCT_JSON_and_regexp(...
        nonBiasRcttid, rctDir, Bso, L)
      assert(ischar(nonBiasRcttid) & ~strcmp(nonBiasRcttid, 'BIAS'))

      RctdCaMap = containers.Map();

      %==========
      % BIAS RCT
      %==========
      [biasRctPath, ~, ~] = read_BIAS_RCT_JSON(rctDir, L);
      RctdCaMap('BIAS') = {...
          bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log(...
          'BIAS', biasRctPath, L) ...
          };

      %==============
      % Non-BIAS RCT
      %==============
      % Find path to RCT.
      settingKey     = bicas.proc.L1L2.cal.rct.RctData.RCTD_METADATA_MAP(...
        nonBiasRcttid).filenameRegexpSettingKey;

      filenameRegexp = Bso.get_fv(settingKey);
      filePath       = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        rctDir, filenameRegexp, L);

      % Read RCT file.
      RctdCa         = {...
        bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log(...
        nonBiasRcttid, filePath, L) ...
        };

      % NOTE: Placing all non-BIAS RCT data inside 1x1 cell arrays so that
      % they are stored analogously with when using L1R ZVCTI1+GACT.
      RctdCaMap(nonBiasRcttid) = RctdCa;
    end



    % (1) Load one BIAS RCT by regular expression.
    % (2) Load one or multiple non-BIAS RCT(s) for one selected RCTTID using
    %     ZVCTI, GACT, and ZV "BW".
    %
    %
    % IMPLEMENTATION NOTE
    % ===================
    % May load MULTIPLE RCTs with the same RCTTID, but will only load those RCTs
    % which are actually needed, as indicated by ZVCTI and ZV "BW". This is
    % necessary since GACT may reference unnecessary RCTs of types not
    % recognized by BICAS (LFR's ROC-SGSE_CAL_RCT-LFR-VHF_V01.cdf /2019-12-16),
    % and which are therefore unreadable by BICAS (BICAS will crash).
    %
    %
    % ARGUMENTS
    % =========
    % nonBiasRcttid
    %       RCTTID for one *non-BIAS* RCT.
    % gact
    %       1D cell array of strings.
    % zvcti
    % zv_BW
    %       Either
    %       (1) [] (as for TDS data), or
    %       (2) LFR input dataset zVariable BW.
    %
    %
    % RETURN VALUE
    % ============
    % RctdCaMap
    %       Returns containers.Map that can be used to initialize an instance
    %       of bicas.proc.L1L2.cal.Cal.
    %
    function RctdCaMap = find_read_RCTs_by_RCT_JSON_and_ZVCTI_GACT(...
        nonBiasRcttid, rctDir, gact, zvcti, zv_BW, L)

      assert(ischar(nonBiasRcttid) & ~strcmp(nonBiasRcttid, 'BIAS'))

      RctdCaMap = containers.Map();

      %==========
      % BIAS RCT
      %==========
      [biasRctPath, ~, ~] = bicas.proc.L1L2.cal.rct.findread.read_BIAS_RCT_JSON(rctDir, L);
      RctdCaMap('BIAS') = {...
          bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log(...
          'BIAS', biasRctPath, L) ...
          };

      %=========================================
      % Read potentially MULTIPLE NON-BIAS RCTs
      %=========================================
      NonBiasRctdCa = bicas.proc.L1L2.cal.rct.findread.find_read_RCTs_by_ZVCTI_GACT(...
        nonBiasRcttid, rctDir, ...
        gact, zvcti, ...
        zv_BW, L);
      RctdCaMap(nonBiasRcttid) = NonBiasRctdCa;
    end



    % Determine the path to the RCT that should be used according to
    % algorithm specified in the documentation(?). If there are multiple
    % matching candidates, choose the latest one as indicated by the
    % filename (last one in list of sorted strings, where the strings are
    % filenames).
    %
    %
    % IMPLEMENTATION NOTES
    % ====================
    % Useful to have this as separate functionality so that the chosen RCT
    % to use can be explicitly overridden via e.g. settings.
    % --
    % NOTE: Only public due to automatic testing.
    %
    function path = get_RCT_path_by_regexp(rctDir, filenameRegexp, L)
      % PROPOSAL: Better name.
      %   ~path, ~file, ~select
      %   find_RCT_by_regexp
      %   find_select_RCT_by_regexp
      %   NOTE: Does not read the file.
      %   NOTE: Cf. bicas.proc.L1L2.cal.rct.findread.find_read_RCTs_by_RCT_JSON_and_regexp()

      %=================================================
      % Find candidate files and select the correct one
      %=================================================
      dirObjectList = dir(rctDir);
      dirObjectList([dirObjectList.isdir]) = [];    % Eliminate directories.
      filenameList = {dirObjectList.name};
      % Eliminate non-matching filenames.
      filenameList(~irf.str.regexpf(filenameList, filenameRegexp)) = [];

      % ASSERTION / WARNING
      if numel(filenameList) == 0
        % ERROR
        error('BICAS:CannotFindRegexMatchingRCT', ...
          ['Can not find any calibration file that matches regular', ...
          ' expression "%s" in directory "%s".'], ...
          filenameRegexp, rctDir);
      end
      % CASE: There is at least one candidate file.

      filenameList = sort(filenameList);
      filename     = filenameList{end};    % NOTE: Selects one file out of many.
      path         = fullfile(rctDir, filename);

      if numel(filenameList) > 1
        % WARNING/INFO/NOTICE
        msg = sprintf(...
          ['Found multiple calibration files matching regular', ...
          ' expression "%s"\n in directory "%s".\n', ...
          'Selecting the latest one as indicated by', ...
          ' the filename (sorting strings): "%s".\n'], ...
          filenameRegexp, rctDir, filename);
        for i = 1:numel(filenameList)
          msg = [msg, sprintf('    %s\n', filenameList{i})];
        end
        L.log('debug', msg)
      end

      % IMPLEMENTATION NOTE: Not logging which calibration file is
      % selected, since this function is not supposed to actually load the
      % content.
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    % Potentially reads MULTIPLE NON-BIAS RCTs (for the same RCTDID) from
    % filenames indirectly specified by arguments "gact" and "zvcti".
    %
    % ARGUMENTS
    % =========
    % zv_BW
    %       LFR L1/L1R ZV "BW". If it does not exist (i.e. if processing TDS
    %       L1/L1R), then the caller should (!) create a fake one and submit
    %       it (normalize input).
    %
    function RctdCa = find_read_RCTs_by_ZVCTI_GACT(...
        nonBiasRcttid, rctDir, gact, zvcti, zv_BW, L)
      % PROPOSAL: Separate function for extracting filenames from ZVs.

      % ASSERTION
      assert(iscell(gact))
      nGactEntries = irf.assert.sizes(...
        gact,  [-1, 1], ...
        zvcti, [-2, 2], ...
        zv_BW, [-2, 1]);
      assert(all(ismember(zv_BW, [0,1])))

      % Obtain indices into GACT
      % ------------------------
      % NOTE: May exclude some values in "zvcti" due to zv_BW.
      % NOTE: GACT entry index (in MATLAB) is one greater than the value stored
      % in ZVCTI.
      iGactEntryArray = unique(zvcti(logical(zv_BW), 1)) + 1;



      % Pre-allocate array of RCTDs with the same RCTTID.
      RctdCa = cell(nGactEntries, 1);

      % IMPLEMENTATION NOTE: Iterate over those entries in
      % GACT that should be CONSIDERED, i.e. NOT all GACT entries.
      % May therefore legitimately leave some cells in cell array empty.
      for i = 1:numel(iGactEntryArray)
        iGactEntry         = iGactEntryArray(i);
        filePath           = fullfile(rctDir, gact{iGactEntry});
        RctdCa{iGactEntry} = bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log(...
          nonBiasRcttid, filePath, L);
      end

    end



    % For a given RCTTID+RCT file, do the following operations, customized for
    % the type of RCT:
    %   (1) read RCT file,
    %   (2) modify the content when required (e.g. extrapolate TFs), and
    %   (3) log the modified RCT content.
    % Effectively wraps the different RCT-reading functions/methods.
    %
    %
    % IMPLEMENTATION NOTES
    % ====================
    % This method exists to
    % (1) run shared code that should be run when reading any RCT (logging,
    %     modifying data),
    % (2) separate the logging from the RCT-reading code, so that external
    %     code can read RCTs without BICAS.
    %
    function Rctd = read_RCT_modify_log(rcttid, filePath, L)

      L.logf(bicas.proc.L1L2.cal.rct.findread.READING_RCT_PATH_LL, ...
        'Reading RCT (rcttid=%s): "%s"', rcttid, filePath)

      RctdMetadata = bicas.proc.L1L2.cal.rct.RctData.RCTD_METADATA_MAP(rcttid);

      % Call constructor(!) of specified class.
      Rctd = feval(RctdMetadata.className, filePath);

      Rctd.log_RCT(L);
    end



  end    % methods(Static, Access=private)



end
