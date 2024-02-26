%
% Implementation of abstract class for nominal use.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef BicasProcessingAccessImpl < bicas.tools.batch.BicasProcessingAccessAbstract



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)



        % OVERRIDE
        function [varargout] = bicas_main(obj, varargin)
            for i = 1:numel(varargin)
                assert(ischar(varargin{i}))
            end

            [varargout{1:nargout}] = bicas.main(varargin{:});
        end



    end    % methods(Access=public)



end
