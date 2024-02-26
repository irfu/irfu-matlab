%
% Find all possible BPCIs for a given set of DSMDs and SWMs.
%
% NOTE: One can limit the number of SWMs by trimming SwmArray.
% NOTE: Requires DSMDs to not overlap in time for each DSI separately.
%       Therefore typically wants to only supply the latest versions of
%       datasets.
%       Can therefore not handle SOLO_L2_RPW-LFR-SBM1-CWF-E.
%
%
% ARGUMENTS
% =========
% See bicas.tools.batch.autocreate_one_SWM_BPCI().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-01-23.
%
function BpciArray = autocreate_many_BPCIs(...
  DsmdArray, SwmArray, createOutputPathFh)

% ASSERTIONS
solo.adm.assert_no_time_overlap(DsmdArray)
assert(isa(SwmArray, 'bicas.swm.SoftwareMode'))

BpciArray = bicas.tools.batch.BicasProcessingCallInfo.empty(0,1);
for iSwm = 1:numel(SwmArray)

  Swm          = SwmArray(iSwm);
  dsiList      = {Swm.inputsList.dsi};
  dsmdGroupsCa = solo.adm.find_overlapping_DSMD_groups(DsmdArray, dsiList);

  nGroups = numel(dsmdGroupsCa);
  for iGrp = 1:nGroups
    % t2 = tic();
    Bpci = bicas.tools.batch.autocreate_one_SWM_BPCI(...
      Swm, ...
      dsmdGroupsCa{iGrp}, ...
      createOutputPathFh);
    % fprintf('SPEED: autocreate_many_BPCIs(): autocreate_one_SWM_BPCI(): %.1f [s] wall time\n', toc(t2))

    BpciArray = [BpciArray; Bpci];
  end
end
end
