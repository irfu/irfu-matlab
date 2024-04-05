%
% Collection of constants.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef const



  %###################
  %###################
  % STATIC PROPERTIES
  %###################
  %###################
  properties(Constant)

    % Return data structure describing in which data directory structure data
    % for different DSIs is assumed to come from.
    %
    % NOTE: Specifies DSIs (i.e. excluding CDAG).
    SOURCE_DSI_DICT = solo.qli.batch.const.get_DSIs()

  end



  %########################
  % PRIVATE STATIC METHODS
  %########################
  methods(Static, Access=private)



    function SourceDsiDict = get_DSIs()

      SourceDsiDict = dictionary();

      % SourceDsiDict('IRFU') = {
      %   'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS', ...
      %   'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS'...
      % };
      SourceDsiDict('LESIA') = {{
        'SOLO_L2_RPW-TNR-SURV', ...
        'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS', ...
        'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS'...
      }};
      SourceDsiDict('SOAR') =  {{
        'SOLO_L2_MAG-RTN-NORMAL', ...
        'SOLO_L2_MAG-RTN-NORMAL-1-MINUTE', ...
        'SOLO_L2_SWA-PAS-EFLUX', ...
        'SOLO_L2_SWA-PAS-GRND-MOM', ...
      }};

      % ASSERTION: Unique DSIs
      % NOTE: Can not assert misspelled DSIs.
      DsiCa = [SourceDsiDict.values{:}];
      irf.assert.castring_set(DsiCa)
    end



  end



end
