%
% Miscellaneous functions.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef utils   % < handle
    % PROPOSAL: Automatic test code.



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
    end    % properties(SetAccess=immutable)



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)
    end    % methods(Access=public)



    %##########################
    %##########################
    % PRIVATE INSTANCE METHODS
    %##########################
    %##########################
    methods(Access=private)
    end    % methods(Access=private)



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)



      function write_file(filePath, rowsCa)
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
