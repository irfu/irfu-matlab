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



    function [pathsCa, fmdSdnArray] = get_file_paths_FMD_SDNs(obj, dirsCa)

      [pathsCa, FsoiArray] = bicas.tools.batch.get_file_paths(dirsCa);
      fmdSdnArray = [FsoiArray.datenum];
    end



  end    % methods(Access=public)



end
