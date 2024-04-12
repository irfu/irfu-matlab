%
% Class which exists to facilitate automated tests. An instance provides access
% to methods (so far only one) for reading (some things) from the file system.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef FileSystemReaderAbstract   % < handle



  %#########################
  %#########################
  % PUBLIC INSTANCE METHODS
  %#########################
  %#########################
  methods(Abstract)



    [pathsCa, fmdSdnArray] = get_file_paths_FMD_SDNs(obj, dirsCa);



  end



end
