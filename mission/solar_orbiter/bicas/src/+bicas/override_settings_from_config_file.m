%
% Override values in SETTINGS object with values from config file.
%
% NOTE: This function exists mostly so that it can be used by BICAS-external
% code that constructs SETTINGS objects.
%
% NOTE: No return value since SETTINGS is a handle object.
%
% NOTE: Hardcoded settings value source string.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-05-18.
%
function override_settings_from_config_file(configFile, SETTINGS, L)
    
    rowList                 = EJ_library.fs.read_text_file(...
        configFile, '(\r\n|\r|\n)');
    
    ConfigFileSettingsVsMap = bicas.interpret_config_file(rowList, L);
    
    SETTINGS.override_values_from_strings(...
        ConfigFileSettingsVsMap, 'configuration file');
    
end
