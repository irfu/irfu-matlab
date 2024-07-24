%
% Class for defining argument value to bicas.utils.cli.parse_CLI_options().
%
% Immutable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef CliOptionConfig
  % PROPOSAL: Better name
  %   PRO: "CLI" is unnecessary given the location in package "cli".
  % PROPOSAL: Abbreviation
  %   PROPOSAL: COPC (not COC)



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % The option header, including any prefix (e.g. dash) expressed as a
    % regular expression.
    optionHeaderRegexp

    % Default=0. "Interpretation Priority". In case of multiple reg.exp.
    % matches, the priority determines which interpretation is used. If multiple
    % options have the same priority, then assertion error.
    interprPriority

    % String specifying the number of times the option may occur.
    % Permitted alternatives (strings):
    % '0-1'   = Option must occur once or never.
    % '1'     = Option must occur exactly once.
    % '0-inf' = Option may occur any number of times (zero or more).
    %
    % IMPLEMENTATION NOTE: This option exists so that multiple
    % optionHeaderRegexp can be allowed to overlap in their coverage.
    % Regexps can not express negation ("match all of this, except
    % this") which can create problems and this tries to mitigate
    % that.
    occurrenceRequirement

    % The number of option values that must follow the option header.
    nValues
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = CliOptionConfig(...
        optionHeaderRegexp, occurrenceRequirement, nValues, interprPriority)

      irf.assert.castring(optionHeaderRegexp)
      irf.assert.castring(occurrenceRequirement)
      assert(isnumeric(nValues) && nValues >= 0)
      assert(isnumeric(interprPriority))

      obj.optionHeaderRegexp    = optionHeaderRegexp;
      obj.occurrenceRequirement = occurrenceRequirement;
      obj.nValues               = nValues;
      obj.interprPriority       = interprPriority;
    end



  end    % methods(Access=public)



end
