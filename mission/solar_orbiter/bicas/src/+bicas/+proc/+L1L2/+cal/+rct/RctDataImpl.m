%
% See superclass.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef(Abstract) RctDataImpl < bicas.proc.L1L2.cal.rct.RctDataAbstract
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



  %##########################
  %##########################
  % PUBLIC STATIC PROPERTIES
  %##########################
  %##########################
  properties(Constant, Access=public)

    % Map to metadata for all RCTD classes
    % ------------------------------------
    % containers.Map: RCTTID --> struct containing information on every RCTT.
    % Its keys defines the set of RCTTID strings.
    RCTD_METADATA_MAP = bicas.proc.L1L2.cal.rct.RctDataImpl.init_RCTD_METADATA_MAP();
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % NOTE: The constructor reads the RCT (file) on its own, just for obtaining
    %       relevant GAs. The RCT is reloaded again by subclasses for loading
    %       the bulk data.
    function obj = RctDataImpl(filePath)

      % Get specified GA from RCT.
      %
      % NOTE: GAs are asserted to be scalar (one entry), if found.
      %
      %
      % RETURN VALUE
      % ============
      % scalarGaValue
      %       [], if GA can not be found in the RCT.
      %       Not cell array.
      function scalarGaValue = get_scalar_GA(gaName)

        if isfield(Do.GlobalAttributes, gaName)
          % CASE: Found GA

          gaValue = Do.GlobalAttributes.(gaName);
          assert(iscell(gaValue))

          if numel(gaValue) ~= 1
            error(...
              'BICAS:FailedToReadInterpretRCT', ...
              ['Global attribute "%s" in RCT "%s" does not have exactly one', ...
              ' entry. Can therefore not interpret this.'], ...
              gaName, filePath)
          end

          scalarGaValue = gaValue{1};

        else

          % CASE: Did not find GA
          scalarGaValue = [];

        end
      end

      Do = dataobj(filePath);

      ga_Data_version = Do.GlobalAttributes.Data_version;
      assert(isscalar(ga_Data_version))
      assert(ischar(  ga_Data_version{1}))

      obj@bicas.proc.L1L2.cal.rct.RctDataAbstract(...
        filePath, ...
        ga_Data_version{1}, ...
        get_scalar_GA('CAL_ENTITY_NAME'), ...
        get_scalar_GA('CAL_ENTITY_AFFILIATION'), ...
        get_scalar_GA('CAL_EQUIPMENT'))
    end



  end    % methods(Access=public)



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
