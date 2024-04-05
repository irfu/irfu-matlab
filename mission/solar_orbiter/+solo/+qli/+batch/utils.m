%
% Miscellaneous functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils   % < handle



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Function useful for automated tests.
    function write_file(filePath, rowsCa)
      % PROPOSAL: Convert to generic function.
      % TODO-NI: Is there no equivalent function already?!

      fileId = fopen(filePath, 'w');
      for i = 1:numel(rowsCa)
        fprintf(fileId, '%s\n', rowsCa{i});
      end
      fclose(fileId);
    end



  end    % methods(Static)



  %########################
  %########################
  % PRIVATE STATIC METHODS
  %########################
  %########################
  methods(Static, Access=private)
  end    % methods(Static, Access=private)



end
