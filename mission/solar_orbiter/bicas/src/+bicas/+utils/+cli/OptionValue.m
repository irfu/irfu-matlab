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
  %
  % PROPOSAL: Rename optionValues --> optionValuesCa
  %   PRO: Is CA.



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
    optionValues
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = OptionValue(...
        iOptionHeaderCliArgument, optionHeader, optionValues)

      assert(isnumeric(iOptionHeaderCliArgument) && iOptionHeaderCliArgument >= 1)
      irf.assert.castring(optionHeader)
      assert(iscell(optionValues) & iscolumn(optionValues))

      obj.iOptionHeaderCliArgument = iOptionHeaderCliArgument;
      obj.optionHeader             = optionHeader;
      obj.optionValues             = optionValues;
    end



  end    % methods(Access=public)



end
