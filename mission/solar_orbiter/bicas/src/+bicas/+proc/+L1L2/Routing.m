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



    %###################
    %###################
    % STATIC PROPERTIES
    %###################
    %###################
    properties(GetAccess=public, Constant)
        C = bicas.proc.L1L2.Routing.init_const();
    end
    
    
    
    %#####################
    %#####################
    % INSTANCE PROPERTIES
    %#####################
    %#####################
    properties(SetAccess=immutable)
        % bicas.proc.L1L2.SignalSourceId.
        % Where the physical signal in the BLTS ultimately comes from. This is
        % used to determine how the signal should be calibrated.
        src
        
        % bicas.proc.L1L2.SignalDestinationId.
        % How the BLTS should be stored in the datasets.
        dest
    end



    %#########################
    %#########################
    % PUBLIC INSTANCE METHODS
    %#########################
    %#########################
    methods(Access=public)

        
        
        function obj = Routing(ssid, varargin)
            assert(isa(ssid, 'bicas.proc.L1L2.SignalSourceId'))
            obj.src = ssid;
            
            switch numel(varargin)
                case 0
                    assert(ssid.is_ASR())
                    sdid = bicas.proc.L1L2.SignalDestinationId(ssid.value);
                case 1
                    sdid = varargin{1};
                    assert(isa(sdid, 'bicas.proc.L1L2.SignalDestinationId'))
                otherwise
                    error('BICAS:Assertion:IllegalArgument', ...
                        'Illegal number of extra arguments.')
            end
            obj.dest = sdid;
            
        end



    end    % methods(Access=public)



    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Access=private, Static)
        
        
        
        function R = init_const()
            % PROPOSAL: Distinguish between different "channels" for 2.5V Ref
            %           and GND in the source (SSID).
            
            SSID = bicas.proc.L1L2.SignalSourceId.C;
            SDID = bicas.proc.L1L2.SignalDestinationId.C;
            R = bicas.proc.L1L2.AntennaSignalId.get_derived_ASR_constants(...
                @(asid) (bicas.proc.L1L2.Routing(...
                    bicas.proc.L1L2.SignalSourceId(asid))));
                
            R.REF25V_TO_DC_V1    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V1);
            R.REF25V_TO_DC_V2    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V2);
            R.REF25V_TO_DC_V3    = bicas.proc.L1L2.Routing(SSID.REF25V,  SDID.DC_V3);
            R.GND_TO_DC_V1       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V1);
            R.GND_TO_DC_V2       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V2);
            R.GND_TO_DC_V3       = bicas.proc.L1L2.Routing(SSID.GND,     SDID.DC_V3);
            R.UNKNOWN_TO_NOWHERE = bicas.proc.L1L2.Routing(SSID.UNKNOWN, SDID.NOWHERE);
        end
        
        
        
    end



end
