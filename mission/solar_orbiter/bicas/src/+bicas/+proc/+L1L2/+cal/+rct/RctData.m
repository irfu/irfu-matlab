%
% Abstract class of which instances of subclasses represent one RCT TYPE (*not*
% the data stored in an RCT). Subclasses should therefore only need to be
% instantiated once, in principle, and should NOT contain any actual RCT data.
%
% NOTE: BICAS may load multiple RCTs for the same RCT type.
%
% IMPLEMENTATION NOTES
% ====================
% * Subclasses effectively collect code (static methods) associated with each
%   RCT type.
% * Instances of subclasses contain data loaded from one RCT.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) RctData
  % PROPOSAL: Classes for RCT data (not just RCT type).
  %   PRO: BIAS data has many fields.
  %   PRO: More well-defined data structs.
  %   PRO: Automatic assertions.
  %   CON: Structs are modified when bicas.proc.L1L2.cal.Cal uses them, i.e. one
  %        could just as well have classes for the format
  %        bicas.proc.L1L2.cal.Cal uses. ==> Too many classes.
  %     TODO-NI: Where does bicas.proc.L1L2.cal.Cal modify the structs?
  %       Does not RctData.modify_RCT_data() in subclasses do all modification
  %       in RctData.read_RCT_modify_log()? Has the question been obsoleted due
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
  %
  % PROPOSAL: Change name to something indicating a store of data.
  %   NOTE: RctData is a historical name.
  %     PROPOSAL: RCTD=RctData
  %       ~data, ~RCT




  %##########################
  %##########################
  % PUBLIC STATIC PROPERTIES
  %##########################
  %##########################
  properties(Constant, Access=public)

    % Map of singleton RCTT objects
    % -----------------------------
    % containers.Map: RCTTID --> struct containing information on every RCTT.
    % Its keys defines the set of RCTTID strings.
    RCTD_METADATA_MAP = bicas.proc.L1L2.cal.rct.RctData.init_RCTD_METADATA_MAP();
  end



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(GetAccess=public, Constant)

    % Minimum number of expected entries in tabulated transfer functions in
    % RCTs.
    TF_TABLE_MIN_LENGTH = 10;

    % LL = Log Level
    RCT_DATA_LL = 'debug';
  end
  properties(GetAccess=public, SetAccess=immutable)
    % Path to RCT file from which data was loaded.
    filePath
  end




  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = RctData(filePath)
        obj.filePath = filePath;
    end



  end    % methods(Access=public)



  %##################################
  %##################################
  % PUBLIC INSTANCE ABSTRACT METHODS
  %##################################
  %##################################
  methods(Access=public)



    % Custom logging of modified RCT data.
    log_RCT(obj, L);



  end    % methods(Access=public)



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
    [RctRawData] = read_RCT(filePath);



  end    % methods(Static, Abstract)



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)

    % Code to initialize hard-coded static constant RCTD_METADATA_MAP.
    %
    % IMPLEMENTATION NOTE: This data structure includes the filename reg.exp.
    % setting keys since it does not appear that MATLAB allows one to access a
    % "Constant instance field" of a class without instantiating it
    % (bicas.proc.L1L2.cal.rct.RctData subclasses). MATLAB does not have true
    % static variables (constant instance fields are the closest).
    %
    function RctdMetadataMap = init_RCTD_METADATA_MAP()
      RctdMetadataMap = containers.Map();

      % (1) Reference to class.
      % (2) Settings key for value that defines the regular expression that is
      %     used for finding the corresponding RCT(s).
      function S = info(className, filenameRegexpSettingKey)
        S = struct( ...
          'className',                className, ...
          'filenameRegexpSettingKey', filenameRegexpSettingKey);
      end

      RctdMetadataMap('BIAS')     = info('bicas.proc.L1L2.cal.rct.RctDataBias',    'PROCESSING.RCT_REGEXP.BIAS');
      RctdMetadataMap('LFR')      = info('bicas.proc.L1L2.cal.rct.RctDataLfr',     'PROCESSING.RCT_REGEXP.LFR');
      RctdMetadataMap('TDS-CWF')  = info('bicas.proc.L1L2.cal.rct.RctDataTdsCwf',  'PROCESSING.RCT_REGEXP.TDS-LFM-CWF');
      RctdMetadataMap('TDS-RSWF') = info('bicas.proc.L1L2.cal.rct.RctDataTdsRswf', 'PROCESSING.RCT_REGEXP.TDS-LFM-RSWF');
    end



  end    % methods(Static)



end
