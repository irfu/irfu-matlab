%
% Simple class which instances represent either of two things:
% (1) SIGNAL SOURCE: where a particular BLTS ultimately comes from, i.e.
%     ** an ASR (DC single/DC diff/AC diff)
%     ** "2.5 V Ref"
%     ** "GND"
%     ** that its origin is unknown (mux mode unknown)
% OR
% (2) ~SIGNAL STORAGE: how the BLTS should be stored in the dataset (since
%     output datasets are only designed to store data measured data as ASRs),
%     i.e.
%     ** an ASR (DC single/DC diff/AC diff)
%     ** nowhere (mux mode unknown).
% NOTE: One instance of this class represents either one of the two above
% alternatives. "src_dest" should thus be interpreted as "source OR dest".
%
% Immutable.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2019-11-19
%
classdef PhysicalSignalSrcDest
    % NOTE: Only(?) the demultiplexer creates its own instances.
    % PROPOSAL: Make publically un-instantiable. Define a fixed set of legal instances accessible via constants.
    %   PRO: Limiting the number of instantiations speeds up the code. Having assertions is then no downside.
    %   PRO: Can move constants from bicas.proc.L1L2.demuxer to this class.
    %   CON: Unpractical to access via names (long qualifiers, including the class name).
    %       CON: bicas.proc.L1L2.demuxer only uses them for its routings, but
    %            can import them for setting the routings.
    %
    % PROPOSAL: Need acronym for all physical signal sources which is a superset of ASR/AS ID.
    %
    % PROPOSAL: Separate classes for
    %           (1) SignalSourceId
    %               DC single, DC diff, AC diff, GND, 2.5V Ref, Unknown
    %               Implement using AntennaSignalId.
    %           (2) SignalDestinationId
    %               DC single, DC diff, AC diff, Nowhere
    %               Implement using AntennaSignalId.
    %           (3) AntennaSignalId
    %               DC single, DC diff, AC diff
    %   CON: Class implementation would be very similar. Code duplication.
    %       CON: Not conceptually similar.
    %   NOTE: "Signal source"      is almost superset of "signal destination".
    %   NOTE: "Signal source"      is        superset of ASR.
    %   NOTE: "Signal destination" is almost superset of ASR/dataset representation.
    %   NOTE: "Unknown": Always a source.
    %   NOTE: "Nowhere": Always a destination. Is signal storage but not ASR?!
    %   NOTE: Only AntennaSignalId can be the key in a future class for storing
    %         all ASRs (samples).
    %   --
    %   CON-PROPOSAL: One class, but methods for determining which type of object.
    %       PRO: Useful for assertions.
    %   PROPOSAL: Classes can use other class: Composition.
    %   PROPOSAL: All inherit from same superclass: Signal ID



    properties(SetAccess=immutable, GetAccess=public)
        % String constant
        category
        
        % 0x0, 1x1, or 1x2 numeric array with components representing antennas.
        % Its exact interpretation depends on "category".
        % Row/column vector important for comparisons. Therefore well defined.
        antennas
    end
    
    
    
    methods(Access=public)



        % Constructor
        function obj = PhysicalSignalSrcDest(category, antennas)
            
            % ASSERTIONS: antennas
            assert(isnumeric(antennas))
            % NOTE: OK for empty value, [].
            assert(all(ismember(antennas, [1,2,3])))
            if isequal(size(antennas), [0,0])
                % CASE: No antennas
                
                % Do nothing
                
            elseif isequal(size(antennas), [1,1])
                % CASE: single
                
                % Do nothing
                
            elseif isequal(size(antennas), [1,2])
                % CASE: diff
                
                % Implicitly checks that antennas are different.
                assert(antennas(1) < antennas(2))
                
            else
                error(...
                    'BICAS:Assertion:IllegalArgument', ...
                    'Trying to define illegal PhysicalSignalSrcDest.')
            end
            
            % ASSERTION: category
            irf.assert.castring(category)

            % ASSERTIONS: category, antennas
            nAntennas = numel(antennas);
            switch(category)
                case 'DC single'
                    assert(nAntennas == 1)
                case 'DC diff'
                    assert(nAntennas == 2)
                case 'AC diff'
                    assert(nAntennas == 2)
                case '2.5V Ref'
                    assert(nAntennas == 0)
                case 'GND'
                    assert(nAntennas == 0)
                case 'Unknown'
                    assert(nAntennas == 0)
                    % Represents that the source of the BLTS is unknown.
                case 'Nowhere'
                    assert(nAntennas == 0)
                    % Represents that the BLTS should be routed to nowhere.
                otherwise
                    % ASSERTION
                    error(...
                        'BICAS:Assertion:IllegalArgument', ...
                        'Illegal argument category="%s".', category)
            end
            
            % Assign object.
            obj.antennas = antennas;
            obj.category = category;
        end



        function isAsr = is_ASR(obj)
            isAsr = ismember(obj.category, {'DC single', 'DC diff', 'AC diff'});
        end
        
        
        
        function isAc = is_AC(obj)
            isAc = strcmp(obj.category, 'AC diff');
        end



    end    % methods(Access=public)
    
end
