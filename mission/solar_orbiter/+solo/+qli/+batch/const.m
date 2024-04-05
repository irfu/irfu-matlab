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

    %=======================================================================
    % DSIs (i.e. excluding CDAG) of filenames which will be searched for in
    % the respective logs.
    %=======================================================================
    % IRFU_LOGS_DSI_CA = {
    %   'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS', ...
    %   'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS'...
    % }
    LESIA_LOGS_DSI_CA = {
      'SOLO_L2_RPW-TNR-SURV', ...
      'SOLO_L3_RPW-BIA-EFIELD-10-SECONDS', ...
      'SOLO_L3_RPW-BIA-DENSITY-10-SECONDS'...
      }
    SOAR_LOGS_DSI_CA = {
      'SOLO_L2_MAG-RTN-NORMAL', ...
      'SOLO_L2_MAG-RTN-NORMAL-1-MINUTE', ...
      'SOLO_L2_SWA-PAS-EFLUX', ...
      'SOLO_L2_SWA-PAS-GRND-MOM', ...
      }

  end



end
