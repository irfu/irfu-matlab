%
% Collection of RCTDs as used by bicas.proc.L1L2.cal.Cal.
%
% Stores column cell arrays of RCTDs, one array per RCTTID. For non-BIAS, array
% corresponds to GACT. Index iRct-1 corresponds to GACT and ZVCTI(:,1) when
% those are used. Permits empty cell elements to avoid loading RCTs which are
% not used by BICAS processing(!)
%
% IMPLEMENTATION NOTE: Permits (1) adding more than the nominal number of
% RCTTIDs (1x BIAS + 1x non-BIAS) to make the bicas.proc.L1L2.cal.Cal work for
% all types of signals at the same time which could be useful for manual
% experimentation with calibraiton, (2) for historical reasons. Ideally,
% bicas.proc.L1L2.cal.Cal should probably be split up into multiple classes
% after which this class probably makes no sense.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef RctdCollection
  % PROPOSAL: Automatic test code.
  % PROPOSAL: Assert scalar BIAS RctdCa.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(Access=private)
    RctdCaMap
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = RctdCollection()
      obj.RctdCaMap = containers.Map();
    end



    function add_RCTD(obj, rcttid, RctdCa)
      assert(~obj.RctdCaMap.isKey(rcttid))
      assert(ismember(...
        rcttid, ...
        bicas.proc.L1L2.cal.rct.RctData.RCTD_METADATA_MAP.keys))

      assert(iscell(RctdCa) && iscolumn(RctdCa))
      for iRctd = 1:numel(RctdCa)
        Rctd = RctdCa{iRctd};

        % NOTE: Deliberately permits empty cell array elements (avoid loading
        % RCTs specified in GACT but for data which does not exist in any BICAS
        % input dataset).
        assert(isempty(Rctd) | isa(Rctd, 'bicas.proc.L1L2.cal.rct.RctData'))
      end

      obj.RctdCaMap(rcttid) = RctdCa;
    end



    function RctdCa = get_RCTD_CA(obj, rcttid)
      RctdCa = obj.RctdCaMap(rcttid);
    end



    % Convert stored RCTDs (excluding empty cell array elements) to plain cell
    % array of RCTDs. This is useful for constructing GAs.
    function RctdCa = get_global_RCTD_CA(obj)
      % IMPLEMENTATION NOTE: It appears that MATLAB does not permit one to
      % create an empty, typed array of instances of suclasses to an abstract
      % superclass, unless the abstract superclass does not inherit from
      % matlab.mixin.Heterogeneous, which seems ugly. Therefore returning cell
      % array instead.
      %
      % >> bicas.proc.L1L2.cal.rct.RctData.empty(0, 1)
      % Error using bicas.proc.L1L2.cal.rct.RctData.empty
      % Abstract classes cannot be instantiated. Class
      % 'bicas.proc.L1L2.cal.rct.RctData' defines abstract methods and/or
      % properties.

      RctdCa = cell(0, 1);

      RctdCaCa = obj.RctdCaMap.values;
      for i = 1:numel(RctdCaCa)
        for j = 1:numel(RctdCaCa{i})
          Rctd = RctdCaCa{i}{j};
          if ~isempty(Rctd)
            RctdCa{end+1, 1} = Rctd;
          end
        end
      end
    end



  end    % methods(Access=public)



end
