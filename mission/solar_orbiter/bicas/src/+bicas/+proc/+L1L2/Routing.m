%
% Immutable class that represents a particular routing of signals via a
% particular BLTS: (1) where a physical signal comes from, and (2) how it should
% be stored in the dataset object (if at all).
%
% Immutable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef Routing   % < handle
    % PROPOSAL: src --> ssid, dest --> sdid



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
            
            assert(obj.dest.is_ASR())
        end



    end    % methods(Access=public)



end
