%
% Abstract class of which instances of subclasses represent the (massaged,
% prepared) content of one RCT FILE. The *static* components of the subclasses
% also effectively represent one RCTTID each (not instances of subclasses).
%
% NOTE: Subclassing is indirect for testing purposes.
% RctDataAbstract
% -->RctDataTest
%   Mock object for tests.
% -->RctDataImpl
%   Shared superclass for classes which represent nominal RCTDs.
%   Contains some code unsuitable for tests (constructor reads RCT CDF).
%
% NOTE: BICAS may load multiple RCTs for the same RCTTID (LFR).
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) RctDataAbstract
  % PROPOSAL: Use same code/function for reading calibration table, as for
  %           reading dataset (and master CDFs)?
  % PROPOSAL: Create general-purpose read_CDF function which handles indices
  %           correctly (1 vs many records).
  % PROPOSAL: Assert CDF skeleton/master version number.
  % PROPOSAL: Assert skeleton/master.
  %   PRO: Can give better error message when reading the wrong RCT.
  %
  % PROPOSAL: Assert/warn (depending on setting?) when CDF metadata imply that
  %           the RCT zVariables have the wrong units.
  % PROPOSAL: Use utility function for reading every zVariable.
  %   PROPOSAL: Assert units from zVar attributes.
  %
  % PROPOSAL: Move RCTD_METADATA_MAP to bicas.const.



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
    % File name (not path) of RCT file from which data was loaded.
    fileName

    % Global attribute "Data_version" in RCT file as a string.
    % NOTE: Data_version is not always set correctly in RCT.
    %   Ex: SOLO_CAL_RCT-LFR-BIAS_V20190123171020.cdf
    %       SOLO_CAL_RPW-BIAS_V202011191204.cdf
    % NOTE: "scalarGa" refers to that GAs are asserted to be scalar (one GA
    %       entry).
    scalarGa_Data_version
    scalarGa_CAL_ENTITY_NAME
    scalarGa_CAL_ENTITY_AFFILIATION
    scalarGa_CAL_EQUIPMENT
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = RctDataAbstract(...
      filePath, ...
      scalarGa_Data_version, ...
      scalarGa_CAL_ENTITY_NAME, ...
      scalarGa_CAL_ENTITY_AFFILIATION, ...
      scalarGa_CAL_EQUIPMENT)

      obj.fileName                        = irf.fs.get_name(filePath);
      obj.scalarGa_Data_version           = scalarGa_Data_version;
      obj.scalarGa_CAL_ENTITY_NAME        = scalarGa_CAL_ENTITY_NAME;
      obj.scalarGa_CAL_ENTITY_AFFILIATION = scalarGa_CAL_ENTITY_AFFILIATION;
      obj.scalarGa_CAL_EQUIPMENT          = scalarGa_CAL_EQUIPMENT;
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
    % NOTE: Subclass code still converts the RCT TFs slightly:
    %   frequency      : Hz    --> rad/s
    %   phase+amplitude: degrees + dimensionless real value
    %                          --> Z (complex number)
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
    % (bicas.proc.L1L2.cal.rct.RctDataImpl subclasses). MATLAB does not have true
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

      RctdMetadataMap('BIAS')     = info('bicas.proc.L1L2.cal.rct.RctDataBias',    []);
      RctdMetadataMap('LFR')      = info('bicas.proc.L1L2.cal.rct.RctDataLfr',     'PROCESSING.RCT_REGEXP.LFR');
      RctdMetadataMap('TDS-CWF')  = info('bicas.proc.L1L2.cal.rct.RctDataTdsCwf',  'PROCESSING.RCT_REGEXP.TDS-LFM-CWF');
      RctdMetadataMap('TDS-RSWF') = info('bicas.proc.L1L2.cal.rct.RctDataTdsRswf', 'PROCESSING.RCT_REGEXP.TDS-LFM-RSWF');
    end



  end    % methods(Static)



end
