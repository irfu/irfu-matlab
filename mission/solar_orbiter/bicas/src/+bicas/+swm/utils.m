%
% Miscellaneous utility functions in the form of static methods.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Assert that string contains human-readable text.
    function assert_text(str)
      % Require at least one whitespace to prevent confusing value with
      % other values.
      irf.assert.castring_regexp(str, '.* .*')

      irf.assert.castring_regexp(str, '[^<>]*')
    end



    function assert_SWM_CLI_option(swmCliOption)
      irf.assert.castring_regexp(...
        swmCliOption, ...
        bicas.const.SWM_CLI_OPTION_REGEX)
    end



    % NOTE: Really refers to "option body".
    function assert_SIP_CLI_option(sipCliOptionBody)
      irf.assert.castring_regexp(...
        sipCliOptionBody, ...
        bicas.const.SIP_CLI_OPTION_BODY_REGEX)
    end



    % NOTE: Wrapper around global counterpart.
    function assert_DSI(dsi)
      bicas.assert_BICAS_DSI(dsi)

      % ASSERTION: Only using SOLO_* DSIs.
      [sourceName, ~, ~] = solo.adm.disassemble_DATASET_ID(dsi);
      assert(strcmp(sourceName, 'SOLO'))
    end



  end    % methods(Static)



end
