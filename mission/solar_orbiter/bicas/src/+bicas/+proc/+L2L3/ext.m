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
  %   NOTE: Should take advantage of bicas.proc.L2L3.ExternalCodeAbstract.
  % PROPOSAL: Move constants to bicas.const.



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
    function R = calc_EFIELD_SCPOT_DENSITY(Excd, Zv, A)
      arguments
        Excd
        Zv.tt2000
        Zv.VDC_Fpa
        Zv.EDC_Fpa
        A.L
      end
      assert(isa(Excd, 'bicas.proc.L2L3.ExternalCodeAbstract'))

      % =================================
      % Call wrapper around solo.vdccal()
      % =================================
      R1 = bicas.proc.L2L3.ext.calc_EFIELD_SCPOT(...
        Excd, ...
        tt2000 =Zv.tt2000, ...
        VDC_Fpa=Zv.VDC_Fpa, ...
        EDC_Fpa=Zv.EDC_Fpa);

      % =================================
      % Call wrapper around solo.psp2ne()
      % =================================
      % NOTE: The name "NeScpQualityBit" is used by solo.psp2ne() and
      % refers to its other return value "NeScp", i.e. "Scp" only refers
      % to the data the density is based on, but the quality bit only
      % refers to density (and not to SCPOT).
      [NeScpTs, NeScpQualityBitFpa, psp2neCodeVerStr] = ...
        bicas.proc.L2L3.ext.calc_DENSITY(R1.PspTs, Excd, A.L);

      assert(strcmp(R1.PspTs.units,   'V'))
      assert(strcmp(R1.ScpotTs.units, 'V'))
      assert(strcmp(NeScpTs.units,    'cm^-3'))

      %==============================
      % Package function return data
      %==============================
      R = [];
      R.PspVoltFpa         = bicas.utils.FPArray(R1.PspTs.data,    'FILL_VALUE', NaN);
      R.ScpotVoltFpa       = bicas.utils.FPArray(R1.ScpotTs.data,  'FILL_VALUE', NaN);
      R.EdcSrfMvpmFpa      = bicas.utils.FPArray(R1.EdcSrfTs.data, 'FILL_VALUE', NaN);
      R.vdccalCodeVerStr   = R1.vdccalCodeVerStr;
      R.vdccalMatVerStr    = R1.vdccalMatVerStr;
      R.NeScpCm3Fpa        = bicas.utils.FPArray(NeScpTs.data, 'FILL_VALUE', NaN);
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
    function R = calc_EFIELD_SCPOT(Excd, Zv)
      arguments
        Excd
        Zv.tt2000
        Zv.VDC_Fpa
        Zv.EDC_Fpa
      end



      %==========================================
      % Create input variables for solo.vdccal()
      %==========================================
      % NOTE: Should TSeries objects really use TensorOrder=1 and
      % repres={x,y,z}?!! VDC and EDC are not time series of vectors, but
      % of three scalars. Probably does not matter. solo.vdccal() does
      % indeed use VDC.x, EDC.x etc.
      VdcTs = TSeries(...
        EpochTT(Zv.tt2000), Zv.VDC_Fpa.array(single(NaN)), ...
        'TensorOrder', 1, ...
        'repres',      {'x', 'y', 'z'});
      EdcTs = TSeries(...
        EpochTT(Zv.tt2000), Zv.EDC_Fpa.array(single(NaN)), ...
        'TensorOrder', 1, ...
        'repres',      {'x', 'y', 'z'});



      %#################################################################
      % CALL BICAS-EXTERNAL CODE
      %#################################################################
      % NOTE: Not specifying calibration file.
      % ==> Use current official calibration file, hardcoded in
      %     solo.vdccal(), that should be used for official datasets.
      [EdcSrfTs, PspTs, ScpotTs, vdccalCodeVerStr, vdccalMatVerStr] = ...
        Excd.vdccal(VdcTs, EdcTs, []);
      clear VdcTs EdcTs
      %#################################################################



      % ASSERTIONS: Check solo.vdccal() return values.
      irf.assert.sizes(...
        Zv.tt2000,     [-1, 1], ...
        EdcSrfTs.data, [-1, 3], ...
        PspTs.data,    [-1, 1], ...
        ScpotTs.data,  [-1, 1]);
      assert(strcmp(EdcSrfTs.units,            'mV/m'))
      assert(strcmp(EdcSrfTs.coordinateSystem, 'SRF' ))
      assert(strcmp(PspTs.units,               'V'))
      assert(strcmp(ScpotTs.units,             'V'))
      assert(all(Zv.tt2000 == EdcSrfTs.time.ttns))
      assert(all(Zv.tt2000 ==    PspTs.time.ttns))
      assert(all(Zv.tt2000 ==  ScpotTs.time.ttns))
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
      % IMPLEMENTATION NOTE: solo.vdccal() sets EdcSrfTs X component to ZERO,
      % if its input data is non-fill value/non-NaN, and NaN if fill value/NaN.
      % Must therefore check for both zero and NaN.
      %     Ex: Dataset 2020-08-01
      % --
      % NOTE: The X component can never be a measurement value since RPW can
      % not measure E field in the X direction.
      % --
      % TODO-DEC: Is solo.vdccal() returning zero for the X component a
      % solo.vdccal() bug? The value is unknown, rather than assumed to be
      % zero(?).
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
    function [NeScpTs, NeScpQualityBitFpa, psp2neCodeVerStr] = ...
        calc_DENSITY(PspTs, Excd, L)

      %##################################################################
      % CALL BICAS-EXTERNAL CODE
      %##################################################################
      [NeScpTs, NeScpQualityBitTs, psp2neCodeVerStr] = Excd.psp2ne(PspTs);
      %##################################################################

      %===============================================
      % ASSERTIONS: Check solo.psp2ne() return values
      %===============================================
      irf.assert.sizes(...
        PspTs.data,             [-1, 1], ...   % Implicitly checks Epoch's size.
        NeScpTs.data,           [-1, 1], ...
        NeScpQualityBitTs.data, [-1, 1] ...
        );
      assert(all(PspTs.time == NeScpTs.time          ))
      assert(all(PspTs.time == NeScpQualityBitTs.time))

      assert(isfloat(NeScpTs.data))
      % if ~all( (NeScpTs.data > 0) | isnan(NeScpTs.data) )
      %   errorMsg = 'solo.psp2ne() returned non-positive (non-NaN) plasma density.';
      if ~all( (NeScpTs.data >= 0) | isnan(NeScpTs.data) )
        errorMsg = 'solo.psp2ne() returned negative (non-NaN) plasma density.';
        % IMPLEMENTATION NOTE: The real check should probably be to assert
        % positive density values, but this has been TEMPORARILY changed to
        % non-negative values while awaiting an update to solo.psp2.ne() from
        % Jordi Boldu. psp2ne() output for 2023-09-06 and 2023-09-07 currently
        % includes density=0 values.
        % /Erik P G Johansson 2025-10-06
        nZero     = numel(find(      NeScpTs.data == 0));
        nNegative = numel(find(      NeScpTs.data <  0));
        nNan      = numel(find(isnan(NeScpTs.data)));
        nAll      = numel(NeScpTs.data);
        L.log( 'error', errorMsg)
        L.logf('error', '    #Zeroes          = %i (%f%%)', nZero,     100*nZero    /nAll)
        L.logf('error', '    #Negative values = %i (%f%%)', nNegative, 100*nNegative/nAll)
        L.logf('error', '    #NaN             = %i (%f%%)', nNan,      100*nNan     /nAll)
        error(errorMsg)
      end
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
      NeScpQualityBitFpa = bicas.utils.FPArray.floatNan2logical(...
        NeScpQualityBitTs.data);
    end



  end    % methods(Static, Access=private)



end
