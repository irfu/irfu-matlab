%
% Class that represents the data in one GA "MODS" entry for one unique
% combination of DSI and dataset version. One such entry can then be "applied"
% to multiple DSIs via bicas.gamods.Database.
%
% IMMUTABLE.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef VersionEntry   % < handle

    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=private, GetAccess=public)
        dateStr
        bicasVersionStr
        commentsCa
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)

        function obj = VersionEntry(dateStr, bicasVersionStr, commentsCa)
            % ASSERTIONS
            irf.assert.castring_regexp(dateStr, ...
                '20[1-9][0-9]-[0-1][0-9]-[0-3][0-9]')
            % NOTE: Without initial "V".
            irf.assert.castring_regexp(bicasVersionStr, '[0-9]+.[0-9]+.[0-9]+')
            bicas.gamods.VersionEntry.assert_commentsCa(commentsCa)

            obj.dateStr         = dateStr;
            obj.bicasVersionStr = bicasVersionStr;
            obj.commentsCa      = commentsCa(:);
        end



        function obj = add_comments(obj, commentsCa)
            bicas.gamods.VersionEntry.assert_commentsCa(commentsCa)
            obj.commentsCa = [obj.commentsCa; commentsCa(:)];
        end


        function s = get_str(obj)
            assert(obj.is_valid())

            commentsStr = strjoin(obj.commentsCa, ' | ');
            % NOTE: Add "V" before/to BICAS version string.
            s = sprintf('%s -- V%s -- %s', ...
                obj.dateStr, obj.bicasVersionStr, commentsStr);
        end



        function is_valid = is_valid(obj)
            is_valid = ~isempty(obj.commentsCa);
        end



    end



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)

        function assert_commentsCa(commentsCa)
            irf.assert.castring_set(commentsCa)

            for i = 1:numel(commentsCa)
                s = commentsCa{i};
                irf.assert.castring_regexp(s, '[-<=_.,()&:''/ a-zA-Z0-9]+')

                % Check that comments ends with period.
                % NOTE: Besides for consistency, this is useful for checking
                % that comment strings hardcoded over multiple rows (inside a
                % cell array) are not accidentally split up into multiple
                % strings (one per row).
                assert(s(end) == '.')
            end
        end

    end

end
