%
% Class that represents the GA "MODS" information for one DSI.
% Contains list of instances of bicas.gamods.VersionEntry.
% 
%
% MUTABLE. HANDLE CLASS.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef DsiEntry < handle
% PROPOSAL: Assert than all comment strings (separate ones in VersionEntry)
%           are unique.

    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=private)
        GmveAr
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        function obj = DsiEntry()
            obj.GmveAr = bicas.gamods.VersionEntry.empty(0, 1);
        end



        % Add GMVE to GMDE.
        %
        %
        % IMPLEMENTATION NOTE: Does not assert not re-using earlier BICAS
        % version since
        % (1) there have been L3 deliveries made with the same BICAS version
        %     (despite minor BICAS modifications), and
        % (2) one can can not easily assert sorting of BICAS version strings.
        %
        %
        % ARGUMENTS
        % =========
        % GmveNew
        %       Date strings must be equal to or later than the LAST stored GMVE
        %       if there is one.
        %       NOTE: There is no constraint on the BICAS version number(!).
        function add_GMVE(obj, GmveNew)
            
            assert(isa(GmveNew, 'bicas.gamods.VersionEntry'))

            if isempty(obj.GmveAr)
                obj.GmveAr(1, 1) = GmveNew;
            else
                GmveLast = obj.GmveAr(end);

                % Obtain list of all date strings, except the last one.
                nonlastDateStrCa = arrayfun(@(ca) ca.dateStr, ...
                    obj.GmveAr(1:end-1), 'UniformOutput', false);

                assert(~ismember(GmveNew.dateStr, nonlastDateStrCa), ...
                    'The date string Gmve.dateStr="%s" has already been used by a non-last GMVE in this GMDE.', ...
                    GmveNew.dateStr)
                if strcmp(GmveNew.dateStr, GmveLast.dateStr)
                    if strcmp(GmveNew.bicasVersionStr, GmveLast.bicasVersionStr)
                        % NOTE: Overwrite last GMVE in array
                        obj.GmveAr(end, 1) = GmveLast + GmveNew;
                    else
                        obj.GmveAr(end+1, 1) = GmveNew;
                    end

                elseif issorted({GmveLast.dateStr, GmveNew.dateStr})
                    % NOTE: issorted() does not check for equality for strings,
                    % but we know that the strings are not equal at this point.
                    obj.GmveAr(end+1, 1) = GmveNew;

                else
                    error('BICAS:Assertion:IllegalArgument', ...
                        'The date string Gmve.dateStr="%s" has already been used in this GMDE.', ...
                        Gmve.dateStr)
                end
            end
        end



        % Return cell array of strings to be used as value GA MODS for the
        % specified DSI.
        function gaModsStrCa = get_MODS_strings_CA(obj)
            gaModsStrCa = arrayfun(...
                @(Gmve) (Gmve.get_str()), obj.GmveAr, ...
                'UniformOutput', false);
        end



    end



end
