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



    % ~Generic dictionary utility function.
    %
    % Merge arbitrary number of dictionaries into one. Only keep the greatest
    % value when there are key collisions.
    %
    % NOTE: Does not care about the type of keys, only whether they are
    % identical to each other or not.
    %
    function MergedDict = merge_dictionaries_max(DictCa, keyType, valueType)
      assert(iscell(DictCa))

      MergedDict = dictionary(keyType, valueType);

      for iDict = 1:numel(DictCa)
        InputDict  = DictCa{iDict};

        keyArray   = InputDict.keys;
        valueArray = InputDict.values;

        for iKey = 1:numel(keyArray)
          MergedDict = solo.qli.batch.utils.dictionary_set_value_max(...
            MergedDict, keyArray(iKey), valueArray(iKey));
        end
      end
    end



    % ~Generic dictionary utility function.
    %
    % Set one dictionary key value. In case of key collision, use the smaller of
    % the old and new value.
    %
    function Dict = dictionary_set_value_min(Dict, key, newValue)
      if Dict.isKey(key)
        oldValue = Dict(key);
        if oldValue > newValue
          Dict(key) = newValue;
        end
      else
        Dict(key) = newValue;
      end
    end



    % ~Generic dictionary utility function.
    %
    % Set one dictionary key value. In case of key collision, use the larger of
    % the old and new value.
    %
    function Dict = dictionary_set_value_max(Dict, key, newValue)
      if Dict.isKey(key)
        oldValue = Dict(key);
        if oldValue < newValue
          Dict(key) = newValue;
        end
      else
        Dict(key) = newValue;
      end
    end



  end    % methods(Static)



end
