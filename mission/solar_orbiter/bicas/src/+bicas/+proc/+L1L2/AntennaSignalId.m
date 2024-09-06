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
% NOTE: The constant field names used in bicas.proc.L1L2.AntennaSignalId.C
% auto-propagate(!) to constant field names in bicas.proc.L1L2.SignalSourceId.C
% and bicas.proc.L1L2.SignalDestinationId.C.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef AntennaSignalId
  % PROPOSAL: Rename string constants for diffs: V-->E
  %   CON: BIAS specification, Table 4 uses abbreviation "V12" etc. for diffs.
  %
  % PROPOSAL: Rename "s" to something more standard.
  %   id
  %       CON: Does not imply string(?). The object is also named "ID".
  % PROPOSAL: Abolish "s" by using ~map-like object which accepts ASIDs as
  %           keys and using map for implementing SameRowsMap.
  %   PROPOSAL: dictionary?
  % PROPOSAL: Rename category-->categoryId
  %
  % PROPOSAL: Not use term ASR?
  %
  % PROPOSAL: Store ASID constants in dictionary.
  %   PROPOSAL: get_derived_ASR_constants() should return dictionary.
  % PROPOSAL: Store mapping string-->ASID objects separately from collections
  %           of string constants and object constants (ALL_ASID_NAMES_CA,
  %           ALL_ASID_CA).
  %   PROPOSAL: Use dictionary. Obtain sets of strings/objects via methods.
  %
  % PROPOSAL: Private ASID/SSID/SDID constructors. Only instantiate objects for
  %           unique ASID/SSID/SDID ONCE.
  %
  % PROPOSAL: Separate class file for storing ASID/SSID/SDID constants.
  %   PRO: Shorter path to constants.
  %   PRO: Better overview.



  properties(Access=public, Constant)
    % IMPLEMENTATION NOTE: Defining one constant struct, which contains
    % multiple constants as fields. Advantages:
    % (1) Some constants need to be defined using other constants (in this
    %     class) and MATLAB then requires that one uses the full qualifiers,
    %     i.e. bicas.proc.L1L2.AntennaSignalId.* which makes the code very
    %     long. One can also *not* use "import".
    % (2) Makes it possible to access constants through a variable copy of
    %     this constant rather than using the long qualifiers.
    C = bicas.proc.L1L2.AntennaSignalId.init_const()
  end



  properties(SetAccess=immutable, GetAccess=public)
    % String constant that represents unique object. The name is deliberately
    % short since objects themselves are IDs and should preferably be used
    % except when one has to convert to string. The only reason this exists
    % is to be used as a replacement for objects as keys in containers.Map
    % objects. containers.Map does not accept objects as keys.
    s

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



  methods(Access=public)



    % Constructor
    function obj = AntennaSignalId(name, category, antennas)
      assert(isstring(name))

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
      obj.s        = name;
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



  methods(Access=public, Static)



    % Function for quickly generating struct of constants with the same names
    % (field names) as the ASID constants (bicas.proc.L1L2.AntennaSignalId.C),
    % and derived from the corresponding ASID constants. Can be used by other
    % classes that define immutable objects using ASIDs.
    function C = get_derived_ASR_constants(fh)
      ASID = bicas.proc.L1L2.AntennaSignalId.C;
      C    = struct();

      for asidNameCa = ASID.ALL_ASID_NAMES_CA'
        asidName     = asidNameCa{1};
        C.(asidName) = fh(ASID.(asidName));
      end
    end



  end



  methods(Access=private, Static)



    function C = init_const()
      C = struct();

      function add(name, asidCategory, asidAntennas)
        C.(name) = bicas.proc.L1L2.AntennaSignalId(...
          name, asidCategory, asidAntennas);
      end

      % =====================================
      % Add every possible unique ASID object
      % =====================================
      add("DC_V1",  'DC_SINGLE', [1   ]);
      add("DC_V2",  'DC_SINGLE', [2   ]);
      add("DC_V3",  'DC_SINGLE', [3   ]);

      add("DC_V12", 'DC_DIFF',   [1, 2]);
      add("DC_V13", 'DC_DIFF',   [1, 3]);
      add("DC_V23", 'DC_DIFF',   [2, 3]);

      add("AC_V12", 'AC_DIFF',   [1, 2]);
      add("AC_V13", 'AC_DIFF',   [1, 3]);
      add("AC_V23", 'AC_DIFF',   [2, 3]);

      % =======================================
      % Create lists of all unique ASID objects
      % =======================================
      ALL_ASID_NAMES_CA = {};    % All ASID names.
      ALL_ASID_CA       = {};    % All ASID objects.
      for fnCa = fieldnames(C)'
        ALL_ASID_NAMES_CA{end+1, 1} = string(fnCa{1});
        ALL_ASID_CA{      end+1, 1} =     C.(fnCa{1});
      end
      C.ALL_ASID_NAMES_CA = ALL_ASID_NAMES_CA;
      C.ALL_ASID_CA       = ALL_ASID_CA;
    end



  end



end
