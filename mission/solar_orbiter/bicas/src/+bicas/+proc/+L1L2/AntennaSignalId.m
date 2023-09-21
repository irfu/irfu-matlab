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
        C = struct( ...
            'DC_V1',  bicas.proc.L1L2.AntennaSignalId('DC single', [1  ]), ...
            'DC_V2',  bicas.proc.L1L2.AntennaSignalId('DC single', [2  ]), ...
            'DC_V3',  bicas.proc.L1L2.AntennaSignalId('DC single', [3  ]), ...
            ...
            'DC_V12', bicas.proc.L1L2.AntennaSignalId('DC diff',   [1,2]), ...
            'DC_V13', bicas.proc.L1L2.AntennaSignalId('DC diff',   [1,3]), ...
            'DC_V23', bicas.proc.L1L2.AntennaSignalId('DC diff',   [2,3]), ...
            ...
            'AC_V12', bicas.proc.L1L2.AntennaSignalId('AC diff',   [1,2]), ...
            'AC_V13', bicas.proc.L1L2.AntennaSignalId('AC diff',   [1,3]), ...
            'AC_V23', bicas.proc.L1L2.AntennaSignalId('AC diff',   [2,3]) ...
        );
    end



    properties(SetAccess=immutable, GetAccess=public)
        % String constant
        category

        % 1x1, or 1x2 numeric array with components representing antennas.
        % Its exact interpretation depends on "category".
        % Row/column vector important for comparisons. Therefore well defined.
        antennas
    end



    methods(Access=public)



        % Constructor
        function obj = AntennaSignalId(category, antennas)

            % ASSERTIONS: antennas
            assert(isnumeric(antennas))
            % NOTE: OK for empty value, [].
            assert(all(ismember(antennas, [1,2,3])))

            if isequal(size(antennas), [1,1])
                % CASE: single
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

            for asidNameCa = fieldnames(ASID)'
                asidName = asidNameCa{1};
                C.(asidName) = fh(ASID.(asidName));
            end
        end

    end

end
