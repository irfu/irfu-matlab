%
% Implementation of superclass for automated tests only. See superclass.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef FileSystemReaderTest < solo.qli.batch.FileSystemReaderAbstract
  % PROPOSAL: Automatic test code.



  %#####################
  %#####################
  % INSTANCE PROPERTIES
  %#####################
  %#####################
  properties(SetAccess=immutable, GetAccess=private)
    ReturnValuesDict
  end



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function obj = FileSystemReaderTest(ReturnValuesDict)      
      obj.ReturnValuesDict = ReturnValuesDict;

      for i = 1:ReturnValuesDict.numEntries
        key = ReturnValuesDict.keys{i};
        assert(iscell(key))

        value = ReturnValuesDict({key});
        assert(iscell(value))
        assert(numel(value{1}) == 2)
        filePathCa = value{1}{1};
        FmdDt      = value{1}{2};
        assert(iscell(filePathCa))
        assert(isa(FmdDt, 'datetime'))
      end
      
    end



    function [pathsCa, fmdSdnArray] = get_file_paths_FMD_SDNs(obj, dirsCa)
      value       = obj.ReturnValuesDict({dirsCa});
      value       = value{1};

      pathsCa     = value{1};
      fmdSdnArray = value{2};

      % Normalize to column arrays.
      pathsCa     = pathsCa(:);
      fmdSdnArray = fmdSdnArray(:);
    end



  end    % methods(Access=public)



end
