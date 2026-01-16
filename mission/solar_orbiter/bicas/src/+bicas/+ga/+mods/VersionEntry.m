%
% Class that represents the data in one GA "MODS" entry for one unique
% combination of (1) DSID and (2) dataset version. One such entry can then be
% "applied" to multiple DSIDs via bicas.ga.mods.Database. One entry contains
% (1) a date (of a BICAS version),
% (2) a BICAS version number, and
% (3) a list of comments which apply to that BICAS version.
% One entry does NOT contain the DSID or dataset version. That is for the owner
% of the object to store.
%
% IMMUTABLE.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef VersionEntry
  % PROPOSAL: BICAS version date should be obtained from table which translates
  %           BICAS version-->BICAS version date, and should therefore not be
  %           stored in the class.
  %   PROBLEM: Might conflict with possibility of having the multiple BICAS
  %            versions with the same date.
  %     Ex: Quick bugfixes. Has happened once(?).



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=private, GetAccess=public)
    dateStr
    bicasVersionStr
    commentsAr
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = VersionEntry(dateStr, bicasVersionStr, commentsAr)
      % ASSERTIONS
      assert(isstring(dateStr))
      irf.assert.castring_regexp(char(dateStr), ...
        '20[1-9][0-9]-[0-1][0-9]-[0-3][0-9]')

      assert(isstring(bicasVersionStr))
      % NOTE: Version string without initial "V".
      irf.assert.castring_regexp(char(bicasVersionStr), '[0-9]+.[0-9]+.[0-9]+')

      bicas.ga.mods.VersionEntry.assert_commentsAr(commentsAr)

      obj.dateStr         = dateStr;
      obj.bicasVersionStr = bicasVersionStr;
      obj.commentsAr      = commentsAr(:);
    end



    % NOTE: Does not modify the object, but returns a modified object(!).
    function obj = add_comments(obj, commentsAr)
      assert(isstring(commentsAr))
      assert(iscolumn(commentsAr))

      obj = bicas.ga.mods.VersionEntry(...
        obj.dateStr, obj.bicasVersionStr, ...
        [obj.commentsAr; commentsAr]);
    end



    % NOTE: Returns char array (no string object).
    function s = get_str(obj)
      commentsStr = strjoin(obj.commentsAr, ' | ');

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

      obj = bicas.ga.mods.VersionEntry(obj1.dateStr, obj1.bicasVersionStr, ...
        [obj1.commentsAr; obj2.commentsAr]);
    end



  end



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)



    function assert_commentsAr(commentsAr)
      assert(isstring(commentsAr))
      assert(iscolumn(commentsAr))
      assert(numel(unique(commentsAr)) == numel(commentsAr))
      assert(~isempty(commentsAr))

      for i = 1:numel(commentsAr)
        str = commentsAr{i};

        % Check that comment only use expected characters
        % -----------------------------------------------
        % NOTE: No "|" since it is used to separate comments within the same
        %       string.
        % PROPOSAL: Check against more than one whitespace in a row.
        irf.assert.castring_regexp(char(str), '[-<=_.,()&*:''/ a-zA-Z0-9]+')

        % Check that comment does not contain more than one whitespace in a row
        % ---------------------------------------------------------------------
        % This may happen if code does not correctly merge hardcoded strings
        % over multiple rows.
        assert(~contains(str, "  "))

        % Check that comment does not contain more than one minus in a row
        % ----------------------------------------------------------------
        % Double minus is used to separate different sections of a MODS entry.
        assert(~contains(str, "--"))

        % Check that comment begins with permitted character
        % --------------------------------------------------
        % NOTE: Excludes special characters and whitespace.
        % NOTE: Besides for consistency, this is useful for checking
        % that comment strings hardcoded over multiple rows (inside a
        % cell array) are not accidentally split up into multiple
        % strings (one per row).
        firstCharStr = extract(str, 1);
        irf.assert.castring_regexp(char(firstCharStr), '[a-zA-Z0-9]')

        % Check that comments ends with period
        % ------------------------------------
        % NOTE: Besides for consistency, this is useful for checking
        % that comment strings hardcoded over multiple rows (inside a
        % cell array) are not accidentally split up into multiple
        % strings (one per row).
        lastCharStr = extract(str, strlength(str));
        assert(lastCharStr == ".", 'Comment %i does not end with period.', i)
      end
    end



  end



end
