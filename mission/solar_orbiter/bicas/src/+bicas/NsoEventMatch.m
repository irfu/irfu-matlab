%
% Return value from bicas.NsoTable.get_NSO_events_timestamps().
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef NsoEventMatch
  % PROPOSAL: Better name.
  %   ~NSO event
  %   ~extracted
  %   ~matching, match
  %   ~QRCB
  %   ExtractedNsoEvent



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    qrcbAr
    qrcid
    iNsoEvent
  end    % properties(SetAccess=immutable)



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = NsoEventMatch(qrcbAr, qrcid, iNsoEvent)
      assert(iscolumn(qrcbAr)    & islogical(qrcbAr))
      assert(isscalar(qrcid)     & isstring(qrcid))
      assert(isscalar(iNsoEvent) & isnumeric(iNsoEvent))

      obj.qrcbAr    = qrcbAr;
      obj.qrcid     = qrcid;
      obj.iNsoEvent = iNsoEvent;
    end



  end    % methods(Access=public)



end
