%
% Immutable class which instances represent an ASR channel ID (DC single/DC
% diff/AC diff). This may or may not refer to
% (1) a physical source signal, or
% (2) how a signal (from antennas) should be represented in a dataset, since
%     output datasets organize output data sorted as if they were physical
%     signals (though the actual values/samples might actually be e.g. GND,
%     2.5V Ref, or unknown).
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
    function obj = AntennaSignalId(category, antennas)
      % ASSERTIONS: antennas
      assert(isnumeric(antennas))
      % NOTE: Assertion permits empty value, []. Assert vector length
      % later.
      assert(all(ismember(antennas, [1,2,3])))

      if isequal(size(antennas), [1,1])
        % CASE: single (antenna)
        % No assertion

      elseif isequal(size(antennas), [1,2])   % Row vector
        % CASE: diff

        % NOTE: Implicitly checks that antennas are different.
        assert(antennas(1) < antennas(2))

      else
        error(...
          'BICAS:Assertion:IllegalArgument', ...
          'Trying to define illegal AntennaSignalId.')
      end

      % ASSERTION: category
      irf.assert.castring(category)

      % ASSERTIONS: category, antennas
      nAntennas = numel(antennas);
      switch(category)
        case 'DC_SINGLE'
          assert(nAntennas == 1)
        case 'DC_DIFF'
          assert(nAntennas == 2)
        case 'AC_DIFF'
          assert(nAntennas == 2)
        otherwise
          % ASSERTION
          error(...
            'BICAS:Assertion:IllegalArgument', ...
            'Illegal argument category="%s".', category)
      end

      % Assign object.
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
