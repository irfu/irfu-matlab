%
% Implementation of abstract superclass for tests. The constructor sets function
% return values.
%
% NOTE 2023-11-03: Not yet used by any tests.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ExternalCodeTest < bicas.proc.L2L3.ExternalCodeAbstract
    % PROPOSAL: Automatic test code for class itself.
    % PROPOSAL: How set data to be returned from functions?
    %   PROPOSAL: Specify cell arrays of return data.
    %       CON: Difficult to distinguish which RV is which.
    %   PROPOSAL: Struct argument.
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        VdccalRv
        Psp2neRv
    end


    
    methods
        %function obj = ExternalCodeTest(vdccalOutputCa, psp2neOuotputCa)
        function obj = ExternalCodeTest(VdccalRv, Psp2neRv)
            obj.VdccalRv = VdccalRv;
            obj.Psp2neRv = Psp2neRv;
        end
    end
    

    
    %##########################################################
    %##########################################################
    % PUBLIC INSTANCE METHODS THAT OVERRIDE SUPERCLASS METHODS
    %##########################################################
    %##########################################################
    methods(Access=public)

        function varargout = vdccal(obj, varargin)
            % function [DCE_SRF_out, PSP_out, ScPot_out, codeVerStr, matVerStr] = vdccal(VDC_inp, EDC_inp, calFilename)
            assert(nargin  == 3)
            assert(nargout == 5)
            varargout = {...
                obj.VdccalRv.DCE_SRF_out, ...
                obj.VdccalRv.PSP_out, ...
                obj.VdccalRv.ScPot_out, ...
                obj.VdccalRv.codeVerStr, ...
                obj.VdccalRv.matVerStr ...
            };
        end
        
        function varargout = psp2ne(obj, varargin)
            % function [NeScp, NeScpQualityBit, codeVerStr] = psp2ne(PSP)
            assert(nargin  == 1)
            assert(nargout == 3)
            varargout = { ...
                obj.Psp2neRv.NeScp, ...
                obj.Psp2neRv.NeScpQualityBit, ...
                obj.Psp2neRv.codeVerStr, ...
            };
        end

    end    % methods(Access=public)



end
