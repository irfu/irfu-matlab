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
  %   BIAS, RCT, validity
  %   JSON
  %     CON: Implementation detail.
  %   file
  %   --
  %   BRJ  = BIAS RCT JSON
  %   BRJF = BIAS RCT JSON File
  %   BRVF = BIAS RCT Validity File -- IMPLEMENTED
  %   BRVJ = BIAS RCT Validity JSON (also the words in the filename)



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



    % Get an instance of RCTDC needed for a nominal instantiation of
    % bicas.proc.L1L2.cal.Cal.
    function Rctdc = get_nominal_RCTDC(...
        useGactRct, nonBiasRcttid, rctDir, ...
        gact, zvcti, zv_BW, tt2000Begin, tt2000End, L)
      % PROPOSAL: Better name.
      %   nominal
      % PROPOSAL: Make function instantiate bicas.proc.L1L2.cal.Cal.
      %   CON: Bad for testing.
      %     CON: Testing this function is impossible anyway.
      %   CON: Conceptually bad. Makes function "less generic".

      DtDataBegin = datetime(irf.cdf.TT2000_to_datevec(tt2000Begin), 'TimeZone', 'UTCLeapSeconds');
      DtDataEnd   = datetime(irf.cdf.TT2000_to_datevec(tt2000End),   'TimeZone', 'UTCLeapSeconds');

      Rctdc = bicas.proc.L1L2.cal.RctdCollection();

      %==========
      % BIAS RCT
      %==========
      [biasRctPath, brvfPath] = bicas.proc.L1L2.cal.rct.findread.get_BRVF_RCT_path(...
        rctDir, DtDataBegin, DtDataEnd);
      L.logf('info', 'Read "%s".', brvfPath)
      BiasRctd = bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log(...
        'BIAS', biasRctPath, L);
      Rctdc.add_RCTD('BIAS', {BiasRctd});

      %===============
      % Non-BIAS RCTs
      %===============
      if useGactRct
        NonBiasRctdCa = bicas.proc.L1L2.cal.rct.findread.find_read_RCTs_by_ZVCTI_GACT(...
          nonBiasRcttid, rctDir, ...
          gact, zvcti, ...
          zv_BW, L);

        Rctdc.add_RCTD(nonBiasRcttid, NonBiasRctdCa);
      else
        NonBiasRctd = bicas.proc.L1L2.cal.rct.findread.find_read_nonBIAS_RCTs_by_regexp(...
          nonBiasRcttid, rctDir, Bso, L);

        Rctdc.add_RCTD(nonBiasRcttid, {NonBiasRctd});
      end
    end



    % Load one non-BIAS RCT by only using assumptions on filenames (and
    % directory).
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
    %
    % RETURN VALUE
    % ============
    % NonBiasRctd
    %       NOTE: Exactly one RCT (not cell array).
    %
    function NonBiasRctd = find_read_nonBIAS_RCTs_by_regexp(...
        nonBiasRcttid, rctDir, Bso, L)

      assert(ischar(nonBiasRcttid) & ~strcmp(nonBiasRcttid, 'BIAS'))

      % Find path to RCT.
      settingKey     = bicas.proc.L1L2.cal.rct.RctDataImpl.RCTD_METADATA_MAP(...
        nonBiasRcttid).filenameRegexpSettingKey;
      filenameRegexp = Bso.get_fv(settingKey);
      filePath       = bicas.proc.L1L2.cal.rct.findread.get_RCT_path_by_regexp(...
        rctDir, filenameRegexp, L);

      % Read RCT file.
      NonBiasRctd    = bicas.proc.L1L2.cal.rct.findread.read_RCT_modify_log(...
        nonBiasRcttid, filePath, L);
    end



    function [biasRctPath, brvfPath] = get_BRVF_RCT_path(rctDir, DtDataBegin, DtDataEnd)
      [rctFilename, DtValidityBegin, DtValidityEnd, brvfPath] = ...
        bicas.proc.L1L2.cal.rct.findread.read_BRVF(rctDir);

      %-----------------------------------------------
      % ASSERT: RCT specified in BRVF covers the data
      %-----------------------------------------------
      % IMPLEMENTATION NOTE: Kind of a hypothetical test, since the BIAS RCT
      % should cover the entire mission.
      assert(DtValidityBegin <= DtDataBegin, ...
        'BIAS RCT validity according to "%s" begins at %s but does not cover the beginning of data at %s.', ...
        brvfPath, DtValidityBegin, DtDataBegin)
      assert(DtDataEnd      <= DtValidityEnd, ...
        'BIAS RCT validity according to "%s" ends at %s but does not cover the end of data and %s.', ...
        brvfPath, DtValidityEnd, DtDataEnd)



     % Add parent directory to RCT filename.
      biasRctPath = fullfile(rctDir, rctFilename);
    end



    % Read the BRVF and return the content.
    %
    % ARGUMENTS
    % =========
    % DtDataBegin, DtDataEnd
    %       UTC datetime objects for the beginning and end of data so that it
    %       can be checked against the BRVF validity time interval.
    function [rctFilename, DtValidityBegin, DtValidityEnd, brvfPath] = ...
        read_BRVF(rctDir)

      % PROPOSAL: Include checking against validity time interval.
      % PROPOSAL: Exclude adding parent directory to BIAS RCT filename.
      % PROPOSAL: ~Wrapper function which adds parent directory to BIAS RCT
      %           filename and checking of time interval.

      brvfPath = fullfile(rctDir, bicas.const.BRVF_FILENAME);
      irf.assert.file_exists(brvfPath)

      jsonStr    = fileread(brvfPath);
      JsonStruct = jsondecode(jsonStr);

      fnCa = fieldnames(JsonStruct);
      assert(length(fnCa) == 1, ...
        'File "%s" does not reference exactly one BIAS RCT as expected.', ...
        brvfPath)

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
      DtValidityBegin = datetime(validityBegin, 'TimeZone', 'UTCLeapSeconds');
      DtValidityEnd   = datetime(validityEnd,   'TimeZone', 'UTCLeapSeconds');
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



    % Load RCTs from filenames indirectly specified by arguments "gact" and
    % "zvcti".
    %
    % IMPLEMENTATION NOTE
    % ===================
    % May load MULTIPLE RCTs with the same RCTTID, but will only load those
    % RCTs which are actually needed, as indicated byGACT, ZVCTI and ZV "BW".
    % This is necessary since GACT may reference unnecessary RCTs of types not
    % recognized by BICAS (LFR's ROC-SGSE_CAL_RCT-LFR-VHF_V01.cdf /2019-12-16),
    % and which are therefore unreadable by BICAS (BICAS will crash).
    %
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
      %   PRO: Easier to test.

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

      % IMPLEMENTATION NOTE: Iterate over those entries in GACT that should be
      % CONSIDERED, i.e. NOT all GACT entries. May therefore legitimately leave
      % some cells in cell array empty.
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

      RctdMetadata = bicas.proc.L1L2.cal.rct.RctDataImpl.RCTD_METADATA_MAP(rcttid);

      % Call constructor(!) of specified class.
      Rctd = feval(RctdMetadata.className, filePath);

      Rctd.log_RCT(L);
    end



  end    % methods(Static, Access=private)



end
