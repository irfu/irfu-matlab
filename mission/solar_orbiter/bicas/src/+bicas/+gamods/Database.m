%
% Mutable class that stores values for GA "MODS" for all output datasets.
% Intended to make it easy to build corresponding hard-coded constant
% (1) with a lot of overlap between dataset IDs and entries, and
% (2) that conforms to certain format.
%
% MUTABLE.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Database < handle
    % PROPOSAL: Abolish specifying DSIs in constructor. Simply add DsiEntry when
    %           needed.
    %   CON: Useful for asserting that only valid DSIs are specified later.

    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=private)
        DsiMap
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = Database(dsiCa)
            obj.DsiMap = containers.Map('KeyType', 'char', 'ValueType', 'Any');

            for i = 1:numel(dsiCa)
                dsi = dsiCa{i};

                % ASSERTION: dsi is an unused key
                % -------------------------------
                % IMPLEMENTATION NOTE: In principle overkill since later code
                % effectively contains the same assertion but without proper
                % error message), but it is actually useful when configuring
                % hardcoded values manually.
                if obj.DsiMap.isKey(dsi)
                    error('BICAS:Assertion', ...
                        'Map already has a key dsi="%s".', dsi)
                end

                % NOTE: Effectively (additional) assertion on that "dsi" is a
                % valid key.
                obj.DsiMap(dsi) = bicas.gamods.DsiEntry();
            end

        end



        function add_version_entry(obj, dsiCa, ve)
            irf.assert.castring_set(dsiCa)           
            assert(isa(ve, 'bicas.gamods.VersionEntry'))

            for i = 1:numel(dsiCa)
                dsi = dsiCa{i};
                de = obj.DsiMap(dsi);
                de.add_version_entry(ve)
            end
        end



        function gaModsStrCa = get_MODS_strings_CA(obj, dsi)
            de = obj.DsiMap(dsi);
            gaModsStrCa = de.get_MODS_strings_CA();
        end



    end    % methods(Access=public)

end
