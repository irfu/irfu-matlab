%
% Nominal implementation of abstract superclass.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef FileSystemReaderImplementation < solo.qli.batch.FileSystemReaderAbstract
  % PROPOSAL: Automatic test code.



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Access=public)



    function [pathsCa, fmdDtArray] = get_file_paths_FMD_SDNs(obj, dirsCa)

      [pathsCa, FsoiArray] = bicas.tools.batch.get_file_paths(dirsCa);
      if ~isempty(FsoiArray)
        fmdDtArray = datetime([FsoiArray.datenum], 'ConvertFrom', 'datenum');
      else
        fmdDtArray = solo.qli.const.EMPTY_DT_ARRAY;
      end
    end



  end    % methods(Access=public)



end
