%
% Mutable class that stores values for GA "MODS" for multiple output datasets.
% Intended for making it easy to build the corresponding hard-coded constants
% (1) with a lot of overlap between dataset IDs and entries, and
% (2) that conforms to certain format specified by the RCS ICD.
%
% Contains map DSI-->bicas.gamods.DsiEntry.
%
% MUTABLE. HANDLE CLASS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Database < handle
% PROPOSAL: Abolish specifying DSIs in constructor. Simply add DsiEntry when
%           needed.
%   CON: Useful for asserting that only valid DSIs are specified later.
% PROPOSAL: Method for returning the latest BICAS version.
%   PRO: Can use for asserting that it equals the current BICAS version.
%   PROBLEM: Version string (x.y.z) is stored as string in code.

    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=private)
        % Map DSI-->GMDE
        DsiGmdeMap
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        % Constructor. Creates database which is initialized with empty
        % bicas.gamods.DsiEntry for specified DSIs.
        function obj = Database(dsiCa)
            obj.DsiGmdeMap = containers.Map('KeyType', 'char', 'ValueType', 'Any');

            for i = 1:numel(dsiCa)
                dsi = dsiCa{i};

                % ASSERTION: dsi is an unused key
                % -------------------------------
                % IMPLEMENTATION NOTE: In principle overkill since later code
                % effectively contains the same assertion but without proper
                % error message), but it is actually useful when configuring
                % hardcoded values manually.
                if obj.DsiGmdeMap.isKey(dsi)
                    error('BICAS:Assertion', ...
                        'Map already has a key dsi="%s".', dsi)
                end

                % NOTE: Effectively (additional) assertion on that "dsi" is a
                % valid key.
                obj.DsiGmdeMap(dsi) = bicas.gamods.DsiEntry();
            end

        end



        % Add one GMVE to multiple DSIs.
        function add_GMVE(obj, dsiCa, Gmve)
            irf.assert.castring_set(dsiCa)
            assert(isa(Gmve, 'bicas.gamods.VersionEntry'))

            for i = 1:numel(dsiCa)
                dsi = dsiCa{i};
                Gmde = obj.DsiGmdeMap(dsi);
                Gmde.add_GMVE(Gmve)
            end
        end



        % Return cell array of strings to be used as value GA MODS for the
        % specified DSI.
        function gaModsStrCa = get_MODS_strings_CA(obj, dsi)
            Gmde        = obj.DsiGmdeMap(dsi);
            gaModsStrCa = Gmde.get_MODS_strings_CA();
        end



    end    % methods(Access=public)



end
