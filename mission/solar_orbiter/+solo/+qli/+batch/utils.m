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



    function Config = read_config_file(configFilePath)
      % ====================================
      % RATIONALE FOR THE CONFIG FILE FORMAT
      % ====================================
      % Absolute path to IRF logo (not relative)
      %   PRO: Does not need to convert relative --> absolute path.
      %     PRO: A file-reading function should not make such conversion.
      %       PRO: Makes testing easier w.r.t. that field.
      %   PRO: Can put IRF logo outside of git repo.
      % Excluding output directory.
      %   PRO: It is useful to modify this value manually when testing setup.
      % Excluding output log directory.
      %   PRO: It is useful to modify this value manually when testing setup.

      irf.log('n', sprintf('Loading config file "%s"', configFilePath))



      bytesArray = irf.fs.read_file(configFilePath);
      % NOTE: String (char array) must be row array.
      s          = char(bytesArray');
      JsonConfig = jsondecode(s);

      assert(...
        numel(fieldnames(JsonConfig)) == 6, ...
        'Config file "%s" has the wrong number of keys.', configFilePath)

      % ========================================================
      % Convert JSON file contents to struct suitable for MATLAB
      % ========================================================
      % NOTE: Includes using different naming conventions in JSON file and
      % MATLAB code.
      % NOTE: Is part of the assertion on the format (check exact set of
      % keys/field names).

      Config = [];
      Config.vhtDir        = JsonConfig.vhtDir;
      Config.irfLogoPath   = JsonConfig.irfLogoPath;
      Config.datasetDirsCa = JsonConfig.datasetDirs;
      Config.soloDbDirPath = JsonConfig.solo_db_init_local_file_db;
      Config.fmdQliDir     = JsonConfig.fmdQliDir;

      Config.LogFileDirPatternDict = dictionary();
      fnCa = fieldnames(JsonConfig.logFileDirPatterns);
      for i = 1:numel(fnCa)
        fn = fnCa{i};
        Config.LogFileDirPatternDict(fn) = JsonConfig.logFileDirPatterns.(fn);
      end
    end



  end    % methods(Static)



end
