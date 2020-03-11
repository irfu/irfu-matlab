%
% Simple class which instances can represent either of two things
% (1) where a particular BLTS comes from, i.e.
%     ** an ASR, DC single/DC diff/DC AC
%     ** "2.5 V Ref"
%     ** "GND"
%     ** that its origin is unknown (mux mode unknown)
% or
% (2) how the BLTS should be stored in the dataset (since output datasets are only designed to store data measured data
% as ASRs), i.e.
%     ** an ASR, DC single/DC diff/DC AC
%     ** nowhere (mux mode unknown).
%
% Immutable.
%
%
% DEFINITIONS
% ===========
% See bicas.calib.
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2019-11-19
%
classdef BLTS_src_dest %< handle
    % NOTE: Only(?) the demultiplexer creates its own instances.
    % PROPOSAL: Make publically un-instantiable. Define a fixed set of legal instances accessible via constants.
    %   PRO: Limiting the number of instantiations speeds up the code. Having assertions is then no downside.
    %   CON: Unpractical to access via names instead of numbers.
    %   CON-PROPOSAL: Have the demultiplexer.main use its own pre-defined constants.
    %       NOTE/CON: Slightly impractical for routings which depend on the latching relay.
    %
    % PROPOSAL: Other name that does not reference BLTS.
    %   See bicas.calib BOGIQ.
    %   PRO: Reference to BLTS is confusing.
    %   PROPOSAL: Need acronym for all physical signal sources which is a superset of ASR/AS ID.
    %
    % PROPOSAL: Separate classes for (1) BLTS physical signal source, and (2) BLTS signal representation in dataset.
    %   CON: (2) is subset of (1).
    %   CON-PROPOSAL: Method for whether object represents a destination in dataset.
    %       PRO: Useful for assertions.



    properties(SetAccess=immutable, GetAccess=public)
        % String constant
        category
        
        % 0x0, 1x1, or 1,2 numeric array with components representing antennas.
        % Its exact interpretation depends on "category".
        % Row/column vector important for comparisons. Therefore well defined.
        antennas
    end
    
    
    
    methods(Access=public)
        
        
        
        function obj = BLTS_src_dest(category, antennas)
            
            % ASSERTIONS: antenna
            assert(isnumeric(antennas))
            assert(all(ismember(antennas, [1,2,3])))    % NOTE: OK for empty "antennas".
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
                error('BICAS:BLTS_src_dest:Assertion:IllegalArgument', 'Trying to define illegal BLTS_src_dest.')
            end
            
            % ASSERTION: category
            EJ_library.assert.castring(category)

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
                    error('BICAS:BLTS_src_dest:Assertion:IllegalArgument', 'Illegal argument category="%s".', category)
            end
            
            % Assign object.
            obj.antennas = antennas;
            obj.category = category;
        end



    end    % methods(Access=public)
    
end
