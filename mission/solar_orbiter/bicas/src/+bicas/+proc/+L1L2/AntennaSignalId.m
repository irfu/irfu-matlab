%
% Immutable class which instances represent an ASR channel ID (DC single/DC diff/AC
% diff). This may or may not refer to a physical source signal, or how a signal
% (from antennas) should be represented in a dataset.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef AntennaSignalId



    properties(Access=public, Constant)
        % IMPLEMENTATION NOTE: Defining one constant struct, which contains
        % multiple constants as fields. Advantages:
        % (1) Some constants need to be defined using other constants (in this
        %     class) and MATLAB then requires that one uses the full qualifiers,
        %     i.e. bicas.proc.L1L2.AntennaSignalId.* which makes the code very
        %     long. One can also *not* use "import".
        % (2) Makes it possible to access constants through a variable copy of
        %     this constant rather than using the long qualifiers.
        C = bicas.proc.L1L2.AntennaSignalId.init_const()
    end



    properties(SetAccess=immutable, GetAccess=public)
        % String constant that represents unique object. Name is deliberately
        % short since objects themselves are IDs and should preferably be used
        % except when one has to convert to string. The only reason this exists
        % is to be used as a replacement for objects as keys in containers.Map
        % objects. containers.Map does not accept objects as keys.
        s
        
        % String constant that represents the type of signal (DC/AC,
        % single/diff).
        category

        % 1x1, or 1x2 numeric array with components representing antennas.
        % Its exact interpretation depends on "category".
        % Row/column vector important for comparisons. Therefore well defined.
        antennas
    end



    methods(Access=public)



        % Constructor
        function obj = AntennaSignalId(name, category, antennas)
            assert(ischar(name))
            
            % ASSERTIONS: antennas
            assert(isnumeric(antennas))
            % NOTE: Assertion permits empty value, []. Assert vector length
            % later.
            assert(all(ismember(antennas, [1,2,3])))

            if isequal(size(antennas), [1,1])
                % CASE: single (antenna)
                % No assertion

            elseif isequal(size(antennas), [1,2])   % Row vector
                % CASE: diff

                % NOTE: Implicitly checks that antennas are different.
                assert(antennas(1) < antennas(2))

            else
                error(...
                    'BICAS:Assertion:IllegalArgument', ...
                    'Trying to define illegal AntennaSignalId.')
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
                otherwise
                    % ASSERTION
                    error(...
                        'BICAS:Assertion:IllegalArgument', ...
                        'Illegal argument category="%s".', category)
            end

            % Assign object.
            obj.s        = name;
            obj.antennas = antennas;
            obj.category = category;
        end



        function isAc = is_AC(obj)
            isAc = strcmp(obj.category, 'AC diff');
        end



    end    % methods(Access=public)



    methods(Access=public, Static)

        % Function for quickly generating struct of constants with the same
        % names (field names) as the ASID constants, and derived from the
        % corresponding ASID constants.
        function C = get_derived_ASR_constants(fh)
            ASID = bicas.proc.L1L2.AntennaSignalId.C;
            C    = struct;

            for asidNameCa = ASID.ALL_ASID_NAMES_CA'
                asidName = asidNameCa{1};
                C.(asidName) = fh(ASID.(asidName));
            end
        end

    end



    methods(Access=private, Static)



        function C = init_const()
            C = struct();
            
            function add(name, asidCategory, asidAntennas)
                C.(name) = bicas.proc.L1L2.AntennaSignalId(name, asidCategory, asidAntennas);
            end
            
            add('DC_V1',  'DC single', [1  ]);
            add('DC_V2',  'DC single', [2  ]);
            add('DC_V3',  'DC single', [3  ]);
            
            add('DC_V12', 'DC diff',   [1, 2]);
            add('DC_V13', 'DC diff',   [1, 3]);
            add('DC_V23', 'DC diff',   [2, 3]);
            
            add('AC_V12', 'AC diff',   [1, 2]);
            add('AC_V13', 'AC diff',   [1, 3]);
            add('AC_V23', 'AC diff',   [2, 3]);
        
            % Create lists.
            ALL_ASID_NAMES_CA = {};    % All ASID names.
            ALL_ASID_CA       = {};    % All ASID objects.
            for fnCa = fieldnames(C)'
                ALL_ASID_NAMES_CA{end+1, 1} = fnCa{1};
                ALL_ASID_CA{   end+1, 1} = C.(fnCa{1});
            end
            C.ALL_ASID_NAMES_CA = ALL_ASID_NAMES_CA;
            C.ALL_ASID_CA    = ALL_ASID_CA;
        end
        
        
        
    end



end
