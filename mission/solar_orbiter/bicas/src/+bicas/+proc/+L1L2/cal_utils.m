%
% Collection of utility functions used by bicas.proc.L1L2.Cal to reduce its
% size. Only meant to contain static methods.
%
% Selected functions in bicas.proc.L1L2.Cal are meant to be moved here.
% 
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2020-11-05.
%
classdef cal_utils
    % PROPOSAL: More automatic test code.
    %
    % PROPOSAL: Make functions that are only used internally private.
    %   Ex: TF_LF_constant_abs_Z()  (the only one?)



    methods(Static)
        
        
        
        % Given a sequence of Epoch values, determine for each value which
        % calibration time index should be used. The caller will have to decide
        % which sequence of data that should be calibrated together (e.g. if
        % calibration time changes in the middle of CWF), and which Epoch values
        % should be used to determine calibration time (e.g. first Epoch value
        % for a snapshot determines entire snapshot).
        %
        %
        % ARGUMENTS AND RETURN VALUES
        % ===========================
        % Epoch
        %       Column vector with Epoch values.
        % CalibEpochList
        %       Column vector of monotonically increasing timestamps ("Epoch
        %       format"). In practice intended to be Bias.epochL or Bias.epochH.
        % iCalib
        %       Array. iCalibList(i) = calibration time index for Epoch(i).
        %
        function [iCalib] = get_calibration_time(Epoch, CalibEpochList)
            
            % ASSERTIONS
            bicas.utils.assert_ZV_Epoch(Epoch)
            bicas.utils.assert_ZV_Epoch(CalibEpochList)
            % IMPLEMENTATION NOTE: Does not work if CalibEpochList is empty,
            % since discretize behaves differently for scalar second argument.
            assert(~isempty(CalibEpochList))
            
            % IMPLEMENTATION NOTE: "discretize" by itself returns NaN for Epoch
            % values outside the outermost edges. Therefore (1) must add upper
            % edge "Inf", (2) asserts non-Nan afterwards.
            % IMPLEMENTATION NOTE: "discretize" behaves differently for scalar
            % second argument. Adding edges at infinity hides this problem. If
            % one does not add infinities and uses a scalar edge list, then one
            % has to treat those cases manually.
            iCalib = discretize(Epoch, [CalibEpochList; Inf], 'IncludedEdge', 'left');
            assert(all(~isnan(iCalib(:))), ...
                'BICAS:SWMProcessing', ...
                ['Can not derive which calibration data to', ...
                ' use for all specified timestamps.'])
        end



        % Modify TABULATED transfer functions, if needed.
        % NOTE: Converts tabulated TF --> Tabulated TF
        %
        % Extrapolate to 0 Hz, if needed, using the value for the lowest
        % tabulated frequency.
        %
        % NOTE: Does NOT remove high frequencies.
        %
        function ModifTabTf = extrapolate_tabulated_TF_to_zero_Hz(TabTf)
            % ASSERTIONS
            assert(TabTf.omegaRps(1) > 0)
            assert(isa(TabTf, 'irf.utils.tabulated_transform'))

            % NOTE: Can not just use the lowest-frequency Z value for 0 Hz since
            % it has to be real (not complex).
            Z1       = TabTf.Z(1);
            signZ0   = sign(real(Z1));
            assert(signZ0 ~= 0, ...
                'BICAS:FailedToReadInterpretRCT:Assertion', ...
                ['Can not extrapolate tabulated inverse transfer function', ...
                ' (ITF) to zero Hz due to ambiguity. real(Z(1)) = 0.'])
            Z0       = abs(Z1) * signZ0;   % Z value at 0 Hz.
            
            omegaRps = [0;  TabTf.omegaRps(:)];
            Z        = [Z0; TabTf.Z(:)       ];
            
            ModifTabTf = irf.utils.tabulated_transform(omegaRps, Z);
        end

        
        
        % ~Modify TF to have constant gain for all frequencies below specified
        % frequency (LF=Low Frequency).
        %
        %
        % NOTE: Special case for tf(0 rad/s) being non-finite.
        %
        %
        % ARGUMENTS
        % =========
        % omegaRps
        %       Arbitrary size. Numeric.
        % omegaLimitRps
        % zLimit
        %       Absolute value is used as gain for frequencies lower than
        %       omegaLimitRps. Primarily intended to be equal to be equal to
        %       tf(omegaLimitRps), but does not have to be.
        %       NOTE: Allowed to be NaN.
        %       NOTE: Function could determine this value but asks the caller to
        %       evaluate it. Since this function is meant to be used in
        %       anonymous functions/function handles (e.g. as a wrapper around a
        %       tf), this avoids calling argument tf every time this function is
        %       called.
        % 
        function Z = TF_LF_constant_abs_Z(tf, omegaRps, omegaLimitRps, zLimit)
            
            % ASSERTIONS
            assert(isa(tf, 'function_handle'))
            assert(isscalar(omegaLimitRps) && ~isnan(omegaLimitRps))
            % IMPLEMENTATION NOTE: Must be able to handle zLimit==NaN in order
            % to be applied to that NaN (BIAS) TF function.
            assert(isscalar(zLimit))
            
            
            
            % NOTE: May evaluate 1/0 at 0 Hz (for e.g. BIAS AC TF), but that
            % should be overwritten afterwards.
            Z = tf(omegaRps);
            b = omegaRps < (omegaLimitRps);
            
            % Handling of special case
            % ========================
            % Z(0 Hz) non-finite (e.g. for BIAS AC ITF).
            % ==> Phase is undetermined.
            % ==> Can not set to non-zero gain with same phase as before.
            % IMPLEMENTATION NOTE: Identifies indices before normalizing, just
            % to be sure that condition stems from input data, not normalization
            % bugs.
            b2 = ~isfinite(Z(b)) & omegaRps(b)==0;
            
            Z(b) = Z(b) ./ abs(Z(b)) * abs(zLimit);
            
            Z(b(b2)) = 0;
        end
        
        
        
        % EXPERIMENTAL
        %
        % Interpolate tabulated TF.
        %
        % Special interpolation for TFs that potentially includes:
        %   ** Separate interpolation for abs(Z) and arg(Z) (unwrapped)
        %   ** Moving average before interpolation
        %   ** Linear interpolation (x)or splines
        %
        %
        % RETURN VALUE
        % ============
        % Z
        %       One per component per omegaPRps. NaN when omegaPRps is outside
        %       the range of the tabulated TF.
        %
%         function Zp = interpolate_TF(omegaRps, Z, omegaEvalRps)
%             % NOTE: spline() extrapolates by default. interp1() does not
%             % (returns NaN).
%             %
%             % PROPOSAL: Interpolate over log(omega)
%             %   CON: Can not do for omega=0.
%             
%             absZ = abs(Z);
%             argZ = unwrap(angle(Z));
%             
%             assert(all(absZ >= 0))
% %             assert(all((min(omegaRps) <= omegaEvalRps) & (omegaEvalRps <= max(omegaRps))), ...
% %                 'Trying to extrapolate outside frequency range of tabulated transfer function.')
%             
%             switch(2)
%                 case 1
%                     Zp = interp1(omegaRps, Z, omegaEvalRps, 'linear');
%                 case 2
%                     bInRange = (min(omegaRps) <= omegaEvalRps) & (omegaEvalRps <= max(omegaRps));
%                     
%                     %absZ = smooth(absZ, 2);
%                     %argZ = smooth(argZ, 2);
%                     
%                     % NOTE:
%                     % ** interp1() returns NaN outside of tabulated range (by
%                     %    default).
%                     % ** spline() extrapolates outside of tabulated range (at
%                     %    least by default).
%                     
%                     absZp = interp1(omegaRps, absZ, omegaEvalRps, 'linear');
%                     %absZp = spline(omegaRps, absZ, omegaPRps);
%                     
%                     argZp = interp1(omegaRps, argZ, omegaEvalRps, 'linear');                   
%                     %argZp = spline(omegaRps, argZ, omegaPRps);
%                     
%                     Zp = absZp .* exp(1i*argZp);
%                     
%                     % Remove extrapolation (from e.g. spline()).
%                     Zp(~bInRange) = NaN;
%             end
%         end
        
        
        
        % Manual test code for bicas.proc.L1L2.cal_utils.interpolate_TF().
        %
%         function interpolate_TF___MTEST()
%             % IMPORTANT NOTE: Tranposing with apostrophe in MATLAB also complex
%             % conjugates the values. Use transpose() to NOT change the complex
%             % values.
%                         
%             close all
%             
%             %Z      = [1,1+1i, 2i, -2+2i, -3i];
%             %omega  = 1:numel(Z);
%             
%             omega  = 1:10;
%             Z = exp((1i+0.2)*omega);
%             
%             %omegaP = [-1:0.1:5];
%             omegaP = [min(omega):0.1:max(omega)];
%             %omegaP = [-2+min(omega):0.1:max(omega)+2];
%             
%             
%             
%             Zp = bicas.proc.L1L2.cal_utils.interpolate_TF(omega, Z, omegaP);
%             
%             %=====================
%             % Plot input & output
%             %=====================
%             subplot(1,3, 1)
%             plot(Z, 'o')
%             hold on
%             plot(Zp, 'linewidth', 3)
%             plot(1i*eps(0), '*')    % Can't (?) plot origin without trick.
%             title('Zp (complex plane)')
%             grid on
%             axis square
%             
%             subplot(1,3, 2)
%             plot(omega, abs(Z), 'o')
%             hold on
%             plot(omegaP, abs(Zp))
%             grid on
%             xlabel('omega')
%             ylabel('|Zp|')
%             
%             subplot(1,3, 3)
%             plot(omega, unwrap(angle(Z)), 'o')
%             hold on
%             plot(omegaP, unwrap(angle(Zp)))
%             grid on
%             xlabel('omega')
%             ylabel('unwrapped arg(Zp)')
%             
%             % Print result
%             %transpose(Zp)
%         end

        
        
        % Evaluate a tabulated transfer function.
        %
        % NOTE: This function is effectively meant to specify how tabulated
        % transfer functions should be interpreted w.r.t. interpolation.
        %
        % NOTE: LFR & TDS RSWF tabulated transfer functions do not strictly
        % cover the frequency range needed for data. Must extrapolate somehow.
        %
        %
        % ARGUMENTS
        % =========
        % TabTf
        % valueOutsideTable
        %       Value used outside of tabulated frequencies.
        %
        %
        % RETURN VALUE
        % ============
        % Z
        %       TF Z values corresponding to omegaRps. NaN for omegaRps values
        %       outside of the range of TabTf, return NaN.
        %
        function Z = eval_tabulated_TF(TabTf, omegaRps, valueOutsideTable)
            % OLD COMMENTS: """"Intended specifically for INVERSE transfer
            % functions. Therefore setting Z=0 for frequencies lower than the
            % table covers.""""
            
            % PROPOSAL: valueOutsideTable only applies within some specified margins (not to infinity).
            % PROPOSAL: Automatic test code. 
            
            assert(isa(TabTf, 'irf.utils.tabulated_transform'))
            assert(isfinite(valueOutsideTable))
            
            % NOTE: interp1 returns NaN for values outside range.
            Z = interp1(TabTf.omegaRps, TabTf.Z, omegaRps, 'linear');
            % CASE: Z == NaN for omegaRps not covered by tabulated TF.
            
            % Set to zero (overwrite) for values above highest tabulated
            % frequency.            
            %bUseTabTf = (omegaRps <= Tf.omegaRps(end));
            %Z(~bUseTabTf) = 0;
            
            Z(~isfinite(Z)) = valueOutsideTable;
            
%             if 0
%                 % ASSERTION
%                 if ~all(isfinite(Z))
%                     % IMPLEMENTATION NOTE: Experience shows that it is useful to
%                     % have an extended error message confirming that the requested
%                     % frequence range is outside the tabulated one, and by how much.
%                     errorMsg = sprintf(...
%                         ['Can not evaluate tabulated transfer function for', ...
%                         ' frequencies OUTSIDE the range of tabulated frequencies.\n', ...
%                         'Range of frequencies for which there are tabulated Z values:\n', ...
%                         '    min(TabulatedTf.omegaRps) = %g\n', ...
%                         '    max(TabulatedTf.omegaRps) = %g\n', ...
%                         'Range of frequencies for which evaluation (interpolation) of Z was attempted:\n', ...
%                         '    min(omegaRps)     = %g\n', ...
%                         '    max(omegaRps)     = %g\n'], ...
%                         min(Tf.omegaRps), ...
%                         max(Tf.omegaRps), ...
%                         min(omegaRps), ...
%                         max(omegaRps));
%                     
%                     error('BICAS:Assertion', errorMsg)
%                 end
%             end

        end
        
        
        
        function itf = create_LFR_BIAS_ITF(...
                itfLfr, itfBias, isAc, acConstGainLowFreqRps)
            % PROPOSAL: Re-purpose into function only for combining BIAS and
            % non-BIAS TFs.
            
            assert(isscalar(isAc), islogical(isAc))
            
            itf = @(omegaRps) (TF_product(omegaRps));            
            
            if isAc()
                % NOTE: Modifies combined LFR+BIAS TF.                
                
                zLimit = itf(acConstGainLowFreqRps);

                itf = @(omegaRps) (bicas.proc.L1L2.cal_utils.TF_LF_constant_abs_Z(...
                    itf, omegaRps, acConstGainLowFreqRps, zLimit));
            end
            
            
            %###################################################################
            % IMPLEMENTATION NOTE: In principle, this function is quite
            % unnecessary for multiplying TFs, but it is useful for putting
            % breakpoints in when debugging TFs which are built from layers of
            % anonymous functions and function handles.
            function Z = TF_product(omegaRps)
                Z = itfLfr(omegaRps) ...
                    .* ...
                    itfBias(omegaRps);
            end
        end
        
        
        
        function log_TF_tabulated(logLevel, tfName, Tf, L)
            % PROPOSAL: Somehow prevent printing unnecessary trailing zeros.
            
            assert(isa(Tf, 'irf.utils.tabulated_transform'))
            
            assert(numel(tfName) <= 38, ...
                'String argument "tfName" is too long. numel(tfName)=%i.', ...
                numel(tfName))
            L.logf(logLevel, '%-38s f={%g--%g} [Hz]', ...
                tfName, ...
                Tf.omegaRps(1  )/(2*pi), ...
                Tf.omegaRps(end)/(2*pi));
        end
        
        
        
        % ARGUMENTS
        % =========
        % freqHzArray  : Array of frequencies for which the TF Z value should be
        %                logged.
        % TfFunchandle : Z(omegaRps).
        %                NOTE: rad/s (omegaRps) as opposed to Hz (freqHzArray).
        %
        function log_TF_function_handle(...
                logLevel, tfName, tfUnit, freqHzArray, tfFuncHandle, L)
            assert(isa(tfFuncHandle, 'function_handle'))
            
            zArray = tfFuncHandle(freqHzArray * 2*pi);
            for i=1:numel(freqHzArray)
                freqHz = freqHzArray(i);
                Z      = zArray(i);
                
                inverseZValueStr = sprintf('1/%10.5f', 1/abs(Z));
                
                %======================================================================================================
                % NOTE 2020-04-30: Execution at ROC fails due to not finding
                % function "phase" for unknown reason.
                % ----------------------------------------------------------
                % Exception.identifier = "MATLAB:UndefinedFunction"
                % Exception.message    = "Undefined function 'phase' for input arguments of type 'double'."
                % Matching MATLAB error message identifier parts (error types derived from Exception1.identifier):
                %     UntranslatableErrorMsgId : Error occurred, but code can not translate the error's MATLAB message
                %     identifier into any of BICAS's internal standard error codes.
                % MATLAB call stack:
                %     row  969, calib.m,                    calib.log_ITF_Z
                %     row  303, calib.m,                    calib.calib
                %     row   68, execute_SWM.m,          execute_SWM
                %     row  455, main.m,                     main_without_error_handling
                %     row  116, main.m,                     main
                % --------------------------------------------------------------------------------------
                % See also
                % https://se.mathworks.com/matlabcentral/answers/408657-which-toolbox-is-phase-in
                % """"phase() as a routine by itself is part of the System
                % Identification Toolbox, in the "obsolete" category. phase() is
                % also a method of the newer iddata() class from the System
                % Identification Toolbox. But what you probably want is angle()
                % followed by unwrap(), which is part of basic MATLAB.""""
                %
                % Therefore using function "angle" instead of "phase".
                %======================================================================================================
                if i ~= 1
                    % Do not print name except on first row.
                    tfName = '';
                end
                
                % Check that string is not too long for neat printouts (all
                % calls from BICAS).
                assert(numel(tfName) <= 46, ...
                    'String argument "tfName" is too long. numel(tfName)=%i.', ...
                    numel(tfName))
                % NOTE: Does not print unwrapped phase. Print ~"not unwrapped",
                % ~"wrapped"?
                L.logf(logLevel, ...
                    '%-46s %7.2f [Hz]: abs(Z)=%8.5f=%12s [%s], phase(Z)=% 6.1f [deg]', ...
                    tfName, freqHz, abs(Z), inverseZValueStr, ...
                    tfUnit, rad2deg(angle(Z)))
            end    % for
        end
        
        
        
        % Utility function for creating string representing 1D vector.
        % Ex: '(3.1416, 2.7183, 1.6180)'
        function s = vector_string(pattern, v)
            assert(~isempty(v))
            irf.assert.vector(v)
            s = sprintf(...
                '(%s)', ...
                strjoin(irf.str.sprintf_many(pattern, v), ', '));
        end
        
        
        
        function assert_iBlts(iBlts)
            assert(ismember(iBlts, [1:5]), ...
                'BICAS:IllegalArgument:Assertion', ...
                'Illegal value iBlts=%g', iBlts)
        end
        
        
        
        function assert_iLsf(iLsf)
            assert(ismember(iLsf,  [1:4]), ...
                'BICAS:IllegalArgument:Assertion', ...
                'Illegal value iLsf=%g.', iLsf)
        end
        
        
        
    end    % methods(Static)
    
    
    
end    % classdef cal_utils
