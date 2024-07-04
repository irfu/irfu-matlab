%
% Abstract class of which instances of subclasses represent one RCT TYPE (*not*
% the data stored in an RCT). Subclasses should therefore only need to be
% instantiated once, in principle, and should NOT contain any actual RCT data.
%
% NOTE: BICAS may load multiple RCTs for the same RCT type.
%
% IMPLEMENTATION NOTES
% ====================
% * Subclasses effectively collect code associated with each RCT type so that
%   RCTs can processed without knowing the type (to some extent).
% * Class can not contain map to singleton objects of subclasses since MATLAB
%   prevents that due recursive definitions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) RctType
  % PROPOSAL: Classes for RCT data (not just RCT type).
  %   PRO: BIAS data has many fields.
  %   PRO: More well-defined data structs.
  %   PRO: Automatic assertions.
  %   CON: Structs are modified when bicas.proc.L1L2.cal.Cal uses them, i.e. one
  %        could just as well have classes for the format
  %        bicas.proc.L1L2.cal.Cal uses. ==> Too many classes.
  %     TODO-NI: Where does bicas.proc.L1L2.cal.Cal modify the structs?
  %       Does not RctType.modify_RCT_data() in subclasses do all modification
  %       in RctType.read_RCT_modify_log()? Has the question been obsoleted due
  %       to refactoring?
  %   PROPOSAL: Convert subclasses to stores of RCT data too.
  %       CON: Can have multiple non-BIAS RCTs loaded. Multiple instances of the
  %            same RCT type has no meaning.
  %
  % PROPOSAL: Use same code/function for reading calibration table, as for reading dataset (and master cdfs)?
  % PROPOSAL: Create general-purpose read_CDF function which handles indices correctly (1 vs many records).
  % PROPOSAL: Assert CDF skeleton/master version number.
  % PROPOSAL: Assert skeleton/master.
  %   PRO: Can give better error message when reading the wrong RCT.
  %
  % PROPOSAL: Assert/warn (depending on setting?) when CDF metadata imply that the RCT zVariables have the wrong units.
  % PROPOSAL: Use utility function for reading every zVariable.
  %   PROPOSAL: Assert units from zVar attributes.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Abstract, Constant, GetAccess=public)
    % Settings key for value that defines the regular expression that is
    % used for finding the corresponding RCT(s).
    filenameRegexpSettingKey
  end
  properties(GetAccess=public, Constant)

    % Minimum number of expected entries in tabulated transfer functions in
    % RCTs.
    TF_TABLE_MIN_LENGTH = 10;

    % LL = Log Level
    RCT_DATA_LL = 'debug';
  end



  %################################
  %################################
  % PUBLIC STATIC ABSTRACT METHODS
  %################################
  %################################
  methods(Static, Abstract)



    % Read RCT file.
    %
    %
    % DESIGN INTENT
    % =============
    % It is useful to be able to read an RCT into memory with as few
    % modifications as possible. Therefore, the returned data structures reflect
    % the content of the RCTs, but not necessarily on the form used by BICAS.
    % Changing the format of data should be done elsewhere, in particular
    % modifications of transfer functions, e.g. extrapolation, cut-offs,
    % inversions. For the same reason, the function should be independent of the
    % class instances (which format the data for BICAS to use).
    % --
    % NOTE: BIAS & LFR RCTs: contain FTFs which are not inverted in this code.
    %       TDS RCTs:        contain ITFs.
    % NOTE: Code still converts RCT TFs slightly:
    %   frequency      : Hz    --> rad/s
    %   phase+amplitude: degrees,dimensionless real value --> Z (complex number)
    %
    [RctData] = read_RCT(filePath);

    % Modify the data structure read by bicas.proc.L1L2.cal.rct.read_RCT()
    % to a data structure that BICAS can use.
    %
    % IMPLEMENTATION NOTE: There is a need to distinguish between (1) the
    % data structures in RCT files, which one may want to inspect manually,
    % or log, and should be quite analogous to the RCT file content, and (2)
    % the calibration data data structures which are convenient for BICAS to
    % use.
    [RctData] = modify_RCT_data(RctData);

    % Custom logging of modified RCT data.
    log_RCT(RctData, L);



  end    % methods(Static)



end
