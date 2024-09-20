function default_anomaly_handling___MTEST
L = bicas.Logger('HUMAN_READABLE', false);

settingValueList = {'ERROR', 'WARNING', 'illegal value'};

for i = 1:numel(settingValueList)
  disp('=======================================================================================================')
  try

    bicas.default_anomaly_handling(L, ...
      settingValueList{i}, 'ANOMALY_HANDLING_SETTING', 'ERROR_WARNING_ILLEGAL_SETTING', ...
      'Anomaly XYZ has been detected.', 'BICAS:default_anomaly_handling___MTEST')
  catch Exc
    fprintf('Exc.message = "%s"\n', Exc.message);
  end
end

disp('=======================================================================================================')
disp('=======================================================================================================')
disp('=======================================================================================================')


bicas.default_anomaly_handling(L, ...
  settingValueList{i}, 'ANOMALY_HANDLING_SETTING', 'other', ...
  'Anomaly XYZ has been detected.', 'BICAS:default_anomaly_handling___MTEST')
end
