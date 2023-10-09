%
% Collection of functions for calling code that is external to BICAS.
%
% Class should at the very minimum provide wrappers around external functions
% that (1) convert input & output, (2) validate the output.
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
%
classdef ext
    % PROPOSAL: Automatic test code.
    % PROPOSAL: Move constants to bicas.const.
    %
    % PROPOSAL: Instantiable class for raw wrappers around external code.
    %   PRO: Can switch out external code. Good for testing.
    
    
    
    properties(GetAccess=private, Constant)
        % Regular expression for the format of version strings from
        % BICAS-external code.
        % Equivalent to: yyyy-mm-ddThh:mm:ss
        CODE_VER_STR_REGEXP = ...
            '[0-9]{4}-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]';
    end



    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static, Access=public)
        
        
        
        % Indirectly call BICAS-external code to calculate
        % (1) EFIELD, SCPOT (solo.vdccal), and from that
        % (2) DENSITY       (solo.psp2ne).
        function R = calc_EFIELD_SCPOT_DENSITY(LfrCwfZv, SETTINGS)
            
            QUALITY_FLAG_minForUse = SETTINGS.get_fv(...
                'PROCESSING.L2_TO_L3.ZV_QUALITY_FLAG_MIN');
            
            % =================================
            % Call wrapper around solo.vdccal()
            % =================================
            R1 = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT(LfrCwfZv, QUALITY_FLAG_minForUse);
            
            % ===================================
            % Caller wrapper around solo.psp2ne()
            % ===================================
            [NeScpTs, NeScpQualityBitFpa, psp2neCodeVerStr] = ...
                bicas.proc.L2L3.ext.calc_DENSITY(R1.PspTs);


            assert(strcmp(R1.PspTs.units,   'V'))
            assert(strcmp(R1.ScpotTs.units, 'V'))
            assert(strcmp(NeScpTs.units,    'cm^-3'))
            
            R = [];
            R.pspVolt            = R1.PspTs.data;
            R.scpotVolt          = R1.ScpotTs.data;
            R.edcSrfMvpm         = R1.EdcSrfTs.data;
            R.vdccalCodeVerStr   = R1.vdccalCodeVerStr;
            R.vdccalMatVerStr    = R1.vdccalMatVerStr;
            R.bNotUsed           = R1.bNotUsed;
            R.neScpCm3           = NeScpTs.data;
            % NOTE: Ignoring return value NeScpQualityBit(Ts) for now. Value is
            %       expected to be used by BICAS later.
            R.NeScpQualityBitFpa = NeScpQualityBitFpa;
            R.psp2neCodeVerStr   = psp2neCodeVerStr;
        end



    end    % methods(Access=public)


        
    %########################
    %########################
    % PRIVATE STATIC METHODS
    %########################
    %########################
    methods(Static, Access=private)



        % Calculate both
        %   (1) ELECTRIC FIELD, and
        %   (2) SPACECRAFT POTENTIALS
        % via the same BICAS-external code solo.vdccal() (still inside
        % irfu-matlab).
        %
        % Largely a wrapper around solo.vdccal().
        %
        % NOTE: Needs to be careful with the units, and incompatible updates to
        % solo.vdccal() without the knowledge of the BICAS author. Therefore
        % uses extra assertions to detect such changes.
        %
        % RETURN VALUE
        % ============
        % R
        %       Struct with multiple variables.
        %       NOTE: Return values are packaged as a struct to provide named
        %       return values and avoid confusing similar return results with
        %       each other.
        %
        function R = calc_EFIELD_SCPOT(Zv, QUALITY_FLAG_minForUse)

            irf.assert.struct(Zv, {'Epoch', 'VDC', 'EDC', 'QUALITY_FLAG_Fpa'}, {})



            %======================================================
            % Create input variables for solo.vdccal()
            % ----------------------------------------
            % Set records to NaN for QUALITY_FLAG below threshold.
            %======================================================
            % NOTE: Unclear how treat QUALITY_FLAG=FV.
            bNotUsedFpa         = Zv.QUALITY_FLAG_Fpa < QUALITY_FLAG_minForUse;
            bNotUsed            = bNotUsedFpa.get_data(false);   % FV = false wise?
            Zv.VDC(bNotUsed, :) = NaN;
            Zv.EDC(bNotUsed, :) = NaN;
            %
            % NOTE: Should TSeries objects really use TensorOrder=1 and
            % repres={x,y,z}?!! VDC and EDC are not time series of vectors, but
            % fo three scalars. Probably does not matter. solo.vdccal() does
            % indeed use VDC.x, EDC.x etc.
            VdcTs = TSeries(...
                EpochTT(Zv.Epoch), Zv.VDC, ...
                'TensorOrder', 1, ...
                'repres',      {'x', 'y', 'z'});
            EdcTs = TSeries(...
                EpochTT(Zv.Epoch), Zv.EDC, ...
                'TensorOrder', 1, ...
                'repres',      {'x', 'y', 'z'});



            %##########################
            % CALL BICAS-EXTERNAL CODE
            %##########################
            % NOTE: Not specifying calibration file.
            % ==> Use current official calibration file, hardcoded in
            %     solo.vdccal(), that should be used for official datasets.
            [EdcSrfTs, PspTs, ScpotTs, vdccalCodeVerStr, vdccalMatVerStr] ...
                = solo.vdccal(VdcTs, EdcTs, []);
            clear VdcTs EdcTs
            %##########################



            % ASSERTIONS: Check solo.vdccal() return values.
            irf.assert.sizes(...
                Zv.Epoch,      [-1, 1], ...
                EdcSrfTs.data, [-1, 3], ...
                PspTs.data,    [-1, 1], ...
                ScpotTs.data,  [-1, 1]);
            assert(strcmp(EdcSrfTs.units,            'mV/m'))
            assert(strcmp(EdcSrfTs.coordinateSystem, 'SRF' ))
            assert(strcmp(PspTs.units,               'V'))
            assert(strcmp(ScpotTs.units,             'V'))
            assert(all(Zv.Epoch == EdcSrfTs.time.ttns))
            assert(all(Zv.Epoch ==    PspTs.time.ttns))
            assert(all(Zv.Epoch ==  ScpotTs.time.ttns))
            irf.assert.castring(vdccalMatVerStr)
            assert(~isempty(vdccalMatVerStr), ...
                ['solo.vdccal() returns an empty vdccalMatVerStr', ...
                ' (string representing the version of the corresponding', ...
                ' .mat file). BICAS therefore needs to be updated.'])
            irf.assert.castring_regexp(vdccalCodeVerStr, bicas.proc.L2L3.ext.CODE_VER_STR_REGEXP)



            %===================================================================
            % Normalize the representation of E-field X-component
            % ---------------------------------------------------
            % Set E_x = NaN, but ONLY if assertion deems that the corresponding
            % information is missing.
            %
            % IMPLEMENTATION NOTE: solo.vdccal() sets antenna 1 values to be
            % zero, if its input data is non-fill value/NaN, but NaN if fill
            % value. Must therefore check for both zero and NaN.
            % Ex: Dataset 2020-08-01
            %===================================================================
            % IMPLEMENTATION NOTE: ismember() does not work for NaN.
            assert(all(EdcSrfTs.data(:, 1) == 0 | isnan(EdcSrfTs.data(:, 1))), ...
                ['EDC for antenna 1 returned from', ...
                ' solo.vdccal() is neither zero nor NaN and can therefore', ...
                ' not be assumed to be unknown anymore.', ...
                ' Verify that this is correct solo.vdccal() behaviour and', ...
                ' (if correct) then update BICAS to handle this.'])
            EdcSrfTs.data(:, 1) = NaN;



            % Build return struct.
            R = [];
            R.PspTs            = PspTs;
            R.ScpotTs          = ScpotTs;
            R.EdcSrfTs         = EdcSrfTs;
            R.vdccalCodeVerStr = vdccalCodeVerStr;
            R.vdccalMatVerStr  = vdccalMatVerStr;
            R.bNotUsed         = bNotUsed;
        end



        % Calculate DENSITY via a BICAS-external code solo.psp2ne() (still
        % inside irfu-matlab).
        %
        % Essentially a wrapper around solo.psp2ne().
        % 
        % NOTE: One needs to be careful with units and incompatible updates to
        % solo.vdccal() without the knowledge of the BICAS author. Therefore
        % uses extra assertions to detect such changes.
        %
        % NOTE: Empirically, some return values are NaN.
        % NOTE: Shortening "SCP" = SCPOT comes from the return variable name in
        % solo.psp2ne().
        %
        % IMPLEMENTATION NOTE: Does not need to check QUALITY_FLAG limit since
        % relies on PSP values for which this has already been done.
        %
        function [NeScpTs, NeScpQualityBitFpa, psp2neCodeVerStr] = calc_DENSITY(PspTs)
            
            %##########################
            % CALL BICAS-EXTERNAL CODE
            %##########################
            [NeScpTs, NeScpQualityBitTs, psp2neCodeVerStr] = solo.psp2ne(PspTs);
            %##########################
            
            %===============================================
            % ASSERTIONS: Check solo.psp2ne() return values
            %===============================================
            irf.assert.sizes(...
                PspTs.data,             [-1, 1], ...
                NeScpTs.data,           [-1, 1], ...
                NeScpQualityBitTs.data, [-1, 1] ...
            );
            assert(all(PspTs.time == NeScpTs.time          ))
            assert(all(PspTs.time == NeScpQualityBitTs.time))

            assert(all( (NeScpTs.data > 0) | isnan(NeScpTs.data)), ...
                'solo.psp2ne() returned non-positive (non-NaN) plasma density.')
            assert(strcmp(NeScpTs.units, 'cm^-3'))

            % NOTE: Not permitting NaN quality bit. Unsure if that is the
            %       best behaviour.
            assert(...
                all(ismember(NeScpQualityBitTs.data, [0, 1])), ...
                'solo.psp2ne() returned illegal NeScpTsQualityBitTs. Contains values which are not 0 or 1.')

            irf.assert.castring_regexp(psp2neCodeVerStr, bicas.proc.L2L3.ext.CODE_VER_STR_REGEXP)

            % ==================================
            % Convert NeScpQualityBitTs --> FPA
            % ==================================
            NeScpQualityBitFpa = bicas.utils.FillPositionsArray.floatNan2logical(...
                NeScpQualityBitTs.data);
        end



    end    % methods(Static, Access=private)



end
