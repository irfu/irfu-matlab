%
% DSMD = DataSet MetaData
%
% Container for metadata for one existing dataset (file).
%
%
% RATIONALE
% =========
% Useful for algorithms that have sets of (metadata on) datasets as input and/or
% output.
% --
% Ex: Select latest version
% Ex: Filter specific DATASET_ID
% Ex: Filter time range
% Ex: Find time-overlapping datasets of same DATASET_ID
%   Ex: For removing duplicates.
% Ex: Find time-overlapping datasets of different DATASET_ID
%   Ex: For plotting/analyzing data
%   Ex: For calling BICAS
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-05-28.
%
classdef DSMD
  % PROPOSAL: Better name,
  %     PROPOSAL: Long name: DatasetMetadata.
  %         PRO: Consistent with other class naming for BICAS.
  %         -model
  %         -metadata
  %         -file
  %
  % PROPOSAL: Make DSMD entirely immutable, including .datasetId. Modify
  %           solo.adm.convert_DSMD_DATASET_ID_to_SOLO() to create new instance instead.
  %
  % PROPOSAL: Add file modification timestamp.
  %   CON: Not desirable for DSMDs which do not correspond to preexisting files on disk.
  %       Ex: bicas.tools.batch.autocreate_input_BPCIs2___speedTest calling
  %           solo.adm.paths_to_DSMD_array() for just filenames.
  %       Ex: bicas.tools.batch.run_BICAS_all_passes() calling
  %           solo.adm.paths_to_DSMD_array() for just filenames.
  %   PROPOSAL: New class = path + modif. timestamp + DSMD(without path/modif.timestamp)
  %
  % PROPOSAL: Remove path from DSMD.
  %   PRO: path is not metadata.
  %   PRO: path is information that is often "correlated" with the metadata in algorithms, but it is not used by algorithms.
  %       Ex: Filter DSMDs ==> Filter paths
  %       Ex: Group  DSMDs ==> Group paths
  %       PROPOSAL: Algorithms should not return DSMDs but indices into the
  %           argument DSMD array.
  %           Ex: Group DSMDs should return groups of indices instead.
  %   PROPOSAL: New class = path, modification date + DSMD(without path and modification date)
  %         DSMDF=DatasetMetadataFile
  %         DSMDM=DatasetMetadataModel
  %     PROPOSAL: New class inherits from DSMD.
  %     PROPOSAL: New class owns one DSMD (composition).
  %         PRO: May want to extract non-file DSMDs from (array of) file DSMDs.
  %
  % PROPOSAL: Add static function for extracting CA of dataset filenames from
  %           array of DSMDs.
  %
  % TODO-DEC: What constitutes equality between DSMDs?
  %           Same DSI and time interval?
  %           isCdag? Version number?
  %           Path? Filename?
  %
  % PROPOSAL: Datetime objects should have uppercase initial
  %     PRO: Consistent with (BICAS) naming conventions.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess = public)
    % NOTE: Mutable because of
    % solo.adm.convert_DSMD_DATASET_ID_to_SOLO() (which
    % modifies DSMDs; not to be confused with
    % solo.adm.convert_DATASET_ID_to_SOLO).
    datasetId
  end
  properties(SetAccess = immutable)
    path

    versionNbr
    isCdag

    % Beginning & end of (presumed) time covered by file. datetime.
    dt1
    dt2
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = DSMD(path, datasetId, versionNbr, isCdag, dt1, dt2)
      assert(ischar(path))
      assert(ischar(datasetId))
      assert(isnumeric(versionNbr))
      assert(isscalar(isCdag))
      assert(islogical(isCdag))
      assert(isa(dt1, 'datetime'))
      assert(isa(dt2, 'datetime'))
      assert(strcmp(dt1.TimeZone, 'UTCLeapSeconds'), 'dt1 is not UTC.')
      assert(strcmp(dt2.TimeZone, 'UTCLeapSeconds'), 'dt2 is not UTC.')
      assert(~isnat(dt1))
      assert(~isnat(dt2))

      obj.path       = path;
      obj.datasetId  = datasetId;
      obj.versionNbr = versionNbr;
      obj.isCdag     = isCdag;
      obj.dt1 = dt1;
      obj.dt2 = dt2;
    end



  end    % methods(Access=public)


end
