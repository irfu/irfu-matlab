%
% Mutable class that stores values for GA "MODS" for multiple output datasets.
% Intended for making it easy to build the corresponding hard-coded constants
% (1) with a lot of overlap between dataset IDs and entries, and
% (2) that conforms to certain format specified by the RCS ICD.
%
% Contains map DSID-->bicas.ga.mods.DsidEntry.
%
% MUTABLE. HANDLE CLASS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Database < handle
  % PROPOSAL: Abolish specifying DSIDs in constructor. Simply add DsidEntry when
  %           needed.
  %   CON: Useful for asserting that only valid DSIDs are specified later.
  % PROPOSAL: Method for returning the latest BICAS version.
  %   PRO: Can use for asserting that it equals the current BICAS version.
  %   PROBLEM: Version string (x.y.z) is stored as string in code.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private, GetAccess=private)
    % Map DSID-->GMDE
    DsidGmdeMap
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    % Constructor. Creates database which is initialized with empty
    % bicas.ga.mods.DsidEntry objects for every specified DSID.
    function obj = Database(dsidCa)
      assert(iscolumn(dsidCa))

      obj.DsidGmdeMap = containers.Map('KeyType', 'char', 'ValueType', 'Any');

      for i = 1:numel(dsidCa)
        dsid = dsidCa{i};

        % ASSERTION: "dsid" is an unused key
        % ---------------------------------
        % IMPLEMENTATION NOTE: In principle overkill since later code
        % effectively contains the same assertion but without proper
        % error message), but it is actually useful when configuring
        % hardcoded values manually.
        if obj.DsidGmdeMap.isKey(dsid)
          error('BICAS:Assertion', ...
            'Argument "dsidCa" contains key dsid="%s" at least twice.', dsid)
        end

        % NOTE: Effectively (additional) assertion on that "dsid" is a
        % valid key.
        obj.DsidGmdeMap(dsid) = bicas.ga.mods.DsidEntry();
      end

    end



    % Add one GMVE to multiple DSIDs.
    function add_GMVE(obj, dsidCa, Gmve)
      assert(iscolumn(dsidCa))
      irf.assert.castring_set(dsidCa)
      assert(isa(Gmve, 'bicas.ga.mods.VersionEntry'))

      for i = 1:numel(dsidCa)
        dsid = dsidCa{i};
        Gmde = obj.DsidGmdeMap(dsid);
        Gmde.add_GMVE(Gmve)
      end
    end



    % Return cell array of strings to be used as value GA MODS for the
    % specified DSID.
    function gaModsStrCa = get_MODS_strings_CA(obj, dsid)
      Gmde        = obj.DsidGmdeMap(dsid);
      gaModsStrCa = Gmde.get_MODS_strings_CA();
    end



  end    % methods(Access=public)



end
