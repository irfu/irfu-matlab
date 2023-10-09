%
% Abstract class which, through a subclasses, provides access to functions
% provide "science calculations" and calibrations by scientists. One can thus
% not write conventional test code for these. This class exists so that the
% dependence on this code can be exchanged for custom implementations for
% tests.
%
% NOTE: The interface must follow the interface of the underlying function,
% although it could technically be allowed to do some minor modifications (e.g.
% alter the order of arguments or return values)
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ExternalCodeAbstract
    % PROPOSAL: Name?
    %   ~external
    %   ~science
    %   ~science calculations
    %   ~science calibration
    %   ~science processing
    %   ~science data
    %   NOTE: Should fit with naming of subclasses.
    %
    % PROPOSAL: How implement interface?
    %   PROPOSAL: nargout = func(varargin)
    %       PRO: Avoids repeating function definitions.
    %       CON: Must manually check that test code uses same interface as real
    %            code.
    %   PROPOSAL: Name all return values and arguments explicitly.


    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Abstract)

        varargout = vdccal(varargin);

        varargout = psp2ne(varargin);

    end

end
