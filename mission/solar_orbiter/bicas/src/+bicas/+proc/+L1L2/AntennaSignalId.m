%
% Immutable class which instances represent an ASID with metadata.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef AntennaSignalId
  % PROPOSAL: Rename string constants for diffs: V-->E
  %   CON: BIAS specification, Table 4 uses abbreviation "V12" etc. for diffs.
  %
  % PROPOSAL: Rename category-->categoryId



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=public)
    % Corresponding SSID.
    ssid

    % String constant that represents the type of signal (single/diff, DC/AC).
    category

    % 1x1, or 1x2 numeric array with components representing antennas.
    % Its exact interpretation depends on "category".
    % Since comparisons must be well defined, permitted values must be well
    % defined to ensure uniqueness:
    % (1) whether it is a column or row vector,
    % (2) the order of antennas in diffs (e.g. E12 vs E21).
    antennas
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Constructor
    function obj = AntennaSignalId(ssid, category, antennas)
      assert(isa(ssid, 'uint8'))
      assert(isnumeric(antennas))
      % NOTE: Assertion permits empty value, []. Assert vector length later.
      assert(all(ismember(antennas, [1,2,3])))

      if isequal(size(antennas), [1,1])
        % CASE: single,  1x1
        % No assertion

      elseif isequal(size(antennas), [1,2])
        % CASE: diff, 1x2 row vector

        % NOTE: Implicitly checks that antennas are different.
        assert(antennas(1) < antennas(2))

      else
        error(...
          'BICAS:Assertion:IllegalArgument', ...
          'Trying to define illegal AntennaSignalId.')
      end

      % ASSERTION: category
      assert(isstring(category))

      % ASSERTIONS: category, antennas
      nAntennas = numel(antennas);
      switch(category)
        case "DC_SINGLE"
          assert(nAntennas == 1)
        case "DC_DIFF"
          assert(nAntennas == 2)
        case "AC_DIFF"
          assert(nAntennas == 2)
        otherwise
          % ASSERTION
          error(...
            'BICAS:Assertion:IllegalArgument', ...
            'Illegal argument category="%s".', category)
      end

      % Assign object.
      obj.ssid     = ssid;
      obj.antennas = antennas;
      obj.category = category;
    end



    function isAc = is_AC(obj)
      isAc = strcmp(obj.category, 'AC_DIFF');
    end



    function isDiff = is_diff(obj)
      isDiff = numel(obj.antennas) == 2;
    end



  end    % methods(Access=public)



end
