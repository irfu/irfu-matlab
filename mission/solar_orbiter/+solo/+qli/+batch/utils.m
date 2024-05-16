%
% Miscellaneous functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils



  %#######################
  %#######################
  % PUBLIC STATIC METHODS
  %#######################
  %#######################
  methods(Static)



    % Function useful for automated tests.
    function write_file(filePath, rowsCa)
      % PROPOSAL: Convert to generic function.
      %   NOTE:
      %     irf.fs.write_file()     writes byte array.
      %     irf.fs.read_file()      reads  byte array.
      %     irf.fs.read_text_file() reads rows.
      %   CON: Not generic enough.
      %     CON: Could be made more generic if backwards-compatible extensions
      %          if needed.

      fileId = fopen(filePath, 'w');
      for i = 1:numel(rowsCa)
        fprintf(fileId, '%s\n', rowsCa{i});
      end
      fclose(fileId);
    end



  end    % methods(Static)



end
