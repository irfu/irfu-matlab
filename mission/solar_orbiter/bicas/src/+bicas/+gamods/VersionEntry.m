%
% Class that represents the data in one GA "MODS" entry for one unique
% combination of (1) DSI and (2) dataset version. One such entry can then be
% "applied" to multiple DSIs via bicas.gamods.Database. One entry contains a
% date (of a BICAS version), a BICAS version number, and a list of comments. One
% entry does NOT contain the DSI or dataset version. That is for the owner of
% the object to store.
%
% IMMUTABLE.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef VersionEntry



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
      % NOTE: Version string without initial "V".
      irf.assert.castring_regexp(bicasVersionStr, '[0-9]+.[0-9]+.[0-9]+')
      bicas.gamods.VersionEntry.assert_commentsCa(commentsCa)

      obj.dateStr         = dateStr;
      obj.bicasVersionStr = bicasVersionStr;
      obj.commentsCa      = commentsCa(:);
    end



    % NOTE: Does not modify the object, but returns a modified object(!).
    function obj = add_comments(obj, commentsCa)
      obj = bicas.gamods.VersionEntry(...
        obj.dateStr, obj.bicasVersionStr, ...
        [obj.commentsCa; commentsCa(:)]);
    end



    function s = get_str(obj)
      commentsStr = strjoin(obj.commentsCa, ' | ');
      % NOTE: Add "V" before/to BICAS version string.
      s = sprintf('%s -- V%s -- %s', ...
        obj.dateStr, obj.bicasVersionStr, commentsStr);
    end



    % Merge two GMVEs with the same date strings and BICAS version number.
    %
    % NOTE: Syntactic sugar.
    %
    function obj = plus(obj1, obj2)
      assert(strcmp(obj1.dateStr,         obj2.dateStr))
      assert(strcmp(obj1.bicasVersionStr, obj2.bicasVersionStr))

      obj = bicas.gamods.VersionEntry(obj1.dateStr, obj1.bicasVersionStr, ...
        [obj1.commentsCa; obj2.commentsCa]);
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
      assert(~isempty(commentsCa))

      for i = 1:numel(commentsCa)
        s = commentsCa{i};
        irf.assert.castring_regexp(s, '[-<=_.,()&:''/ a-zA-Z0-9]+')

        % Check that comments ends with period.
        % NOTE: Besides for consistency, this is useful for checking
        % that comment strings hardcoded over multiple rows (inside a
        % cell array) are not accidentally split up into multiple
        % strings (one per row).
        assert(s(end) == '.', 'Comment %i does not end with period.', i)
      end
    end



  end



end
