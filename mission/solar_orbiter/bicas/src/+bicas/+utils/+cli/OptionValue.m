%
% Class for defining return value from bicas.utils.cli.parse_CLI_options().
%
% Immutable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef OptionValue
  % PROPOSAL: Better name.
  %   NOTE: Has also used term "option occurrence".
  %     ~Occurrence, finding



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable)
    % The location (index) of the option header for which this instance contains
    % information.
    iOptionHeaderCliArgument

    optionHeader
    optionValuesCa
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = OptionValue(...
        iOptionHeaderCliArgument, optionHeader, optionValuesCa)

      assert(isnumeric(iOptionHeaderCliArgument) && iOptionHeaderCliArgument >= 1)
      irf.assert.castring(optionHeader)
      assert(iscell(optionValuesCa) & iscolumn(optionValuesCa))

      obj.iOptionHeaderCliArgument = iOptionHeaderCliArgument;
      obj.optionHeader             = optionHeader;
      obj.optionValuesCa           = optionValuesCa;
    end



  end    % methods(Access=public)



end
