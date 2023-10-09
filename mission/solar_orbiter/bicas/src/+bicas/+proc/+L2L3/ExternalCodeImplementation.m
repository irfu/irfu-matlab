%
% Implementation of abstract superclass that actually proceeds call to nominal
% external code.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ExternalCodeImplementation < bicas.proc.L2L3.ExternalCodeAbstract
    % PROPOSAL: Automatic test code.
    %   PROPOSAL: Code that just calls the external functions to make sure that
    %             the interfaces seems to be consistent.

    %##########################################################
    %##########################################################
    % PUBLIC INSTANCE METHODS THAT OVERRIDE SUPERCLASS METHODS
    %##########################################################
    %##########################################################
    methods(Access=public)

        function varargout = vdccal(obj, varargin)
            [varargout{1:nargout}] = solo.vdccal(varargin{:});
        end

        function varargout = psp2ne(obj, varargin)
            [varargout{1:nargout}] = solo.psp2ne(varargin{:});
        end

    end    % methods(Access=public)

end
