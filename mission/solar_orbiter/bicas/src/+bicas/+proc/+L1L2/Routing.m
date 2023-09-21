%
% Class that represents a particular routing of signals via a particular BLTS:
% Where the physical signal comes from and how it should be stored in the
% datasets objects.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Routing   % < handle



    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        % bicas.proc.L1L2.PhysicalSignalSrcDest.
        % Where the physical signal in the BLTS ultimately comes from. This is
        % used to determine how the signal should be calibrated.
        src
        
        % bicas.proc.L1L2.PhysicalSignalSrcDest.
        % How the BLTS should be stored in the datasets.
        dest
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)

        
        
        function obj = Routing(src, varargin)
            assert(isa(src,  'bicas.proc.L1L2.PhysicalSignalSrcDest'))
            obj.src  = src;
            
            switch numel(varargin)
                case 0
                    obj.dest = src;
                case 1
                    dest = varargin{1};
                    assert(isa(dest, 'bicas.proc.L1L2.PhysicalSignalSrcDest'))
                    obj.dest = dest;
                otherwise
                    error('BICAS:Assertion:IllegalArgument', ...
                        'Illegal number of extra arguments.')
            end
        end



    end    % methods(Access=public)
        
end
