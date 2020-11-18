%
% Collection of utility functions used by bicas.calib to reduce its size. Only
% meant to contain static methods.
%
% Selected functions in bicas.calib are meant to be moved here.
% 
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2020-11-05.
%
classdef calib_utils
    % PROPOSAL: Automatic test code.



    methods(Static)
        
        
        
        % Modify TABULATED transfer functions, if needed.
        % NOTE: Converts Tabulated TF --> Tabulated TF
        %
        % Extrapolate to 0 Hz, if needed.
        %
        % NOTE: Does NOT remove high frequencies.
        %
        function ModifTabTf = extrapolate_tabulated_TF_to_zero_Hz(TabTf)
            % ASSERTIONS
            assert(TabTf.omegaRps(1) > 0)
            assert(isa(TabTf, 'EJ_library.utils.tabulated_transform'))

            % NOTE: Can not just use the lowest-frequency Z value for 0 Hz since
            % it has to be real (not complex).
            Z1       = TabTf.Z(1);
            signZ0   = sign(real(Z1));
            assert(signZ0 ~= 0, ...
                'BICAS:calib:extrapolate_tabulated_TF_to_zero_Hz:FailedToReadInterpretRCT:Assertion', ...
                ['Can not extrapolate tabulated inverse transfer function', ...
                ' (ITF) to zero Hz due to ambiguity. real(Z(1)) = 0.'])
            Z0       = abs(Z1) * signZ0;   % Z value at 0 Hz.
            
            omegaRps = [0;  TabTf.omegaRps(:)];
            Z        = [Z0; TabTf.Z(:)       ];
            
            ModifTabTf = EJ_library.utils.tabulated_transform(omegaRps, Z);
        end

        
        
        % EXPERIMENTAL. NOT USED YET(?)
        %
        % Modify TF to have constant abs(Z) for all frequencies below specified
        % frequency (LF=Low Frequency).
        %
        % NOTE: Function could determine this value but asks the caller for it
        % to avoid calling the tf every time this function is called (e.g. if it
        % is used as a wrapper around a tf).
        % 
        function Z = TF_LF_constant_abs_Z(tf, omegaRps, fLimitHz, zLimit)
            
            assert(isa(tf, 'function_handle'))
            %zLimit = tf(fLimitHz * 2*pi);     % TEMP
            assert(isfinite(zLimit))
            
            Z = tf(omegaRps);
            b = omegaRps < (fLimitHz * 2*pi);
            
            Z(b) = Z(b) ./ abs(Z(b)) * abs(zLimit);

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
        function Zp = interpolate_TF(omegaRps, Z, omegaEvalRps)
            % NOTE: spline() extrapolates by default. interp1() does not (returns
            % NaN).
            %
            % PROPOSAL: Interpolate over log(omega)
            %   CON: Can not do for omega=0.
            
            absZ = abs(Z);
            argZ = unwrap(angle(Z));
            
            assert(all(absZ >= 0))
%             assert(all((min(omegaRps) <= omegaEvalRps) & (omegaEvalRps <= max(omegaRps))), ...
%                 'Trying to extrapolate outside frequency range of tabulated transfer function.')
            
            switch(2)
                case 1
                    Zp = interp1(omegaRps, Z, omegaEvalRps, 'linear');
                case 2
                    bInRange = (min(omegaRps) <= omegaEvalRps) & (omegaEvalRps <= max(omegaRps));
                    
                    %absZ = smooth(absZ, 2);
                    %argZ = smooth(argZ, 2);
                    
                    % NOTE:
                    % ** interp1() returns NaN outside of tabulated range (by
                    %    default).
                    % ** spline() extrapolates outside of tabulated range (at
                    %    least by default).
                    
                    absZp = interp1(omegaRps, absZ, omegaEvalRps, 'linear');
                    %absZp = spline(omegaRps, absZ, omegaPRps);
                    
                    argZp = interp1(omegaRps, argZ, omegaEvalRps, 'linear');                   
                    %argZp = spline(omegaRps, argZ, omegaPRps);
                    
                    Zp = absZp .* exp(1i*argZp);
                    
                    % Remove extrapolation (from e.g. spline()).
                    Zp(~bInRange) = NaN;
            end
        end
        
        
        
        % Manual test code for bicas.calib_utils.interpolate_TF().
        %
        function interpolate_TF___MTEST()
            % IMPORTANT NOTE: Tranposing with apostrophe in MATLAB also complex
            % conjugates the values. Use transpose() to NOT change the complex
            % values.
                        
            close all
            
            %Z      = [1,1+1i, 2i, -2+2i, -3i];
            %omega  = 1:numel(Z);
            
            omega  = 1:10;
            Z = exp((1i+0.2)*omega);
            
            %omegaP = [-1:0.1:5];
            omegaP = [min(omega):0.1:max(omega)];
            %omegaP = [-2+min(omega):0.1:max(omega)+2];
            
            
            
            Zp = bicas.calib_utils.interpolate_TF(omega, Z, omegaP);
            
            %=====================
            % Plot input & output
            %=====================
            subplot(1,3, 1)
            plot(Z, 'o')
            hold on
            plot(Zp, 'linewidth', 3)
            plot(1i*eps(0), '*')    % Can't (?) plot origin without trick.
            title('Zp (complex plane)')
            grid on
            axis square
            
            subplot(1,3, 2)
            plot(omega, abs(Z), 'o')
            hold on
            plot(omegaP, abs(Zp))
            grid on
            xlabel('omega')
            ylabel('|Zp|')
            
            subplot(1,3, 3)
            plot(omega, unwrap(angle(Z)), 'o')
            hold on
            plot(omegaP, unwrap(angle(Zp)))
            grid on
            xlabel('omega')
            ylabel('unwrapped arg(Zp)')
            
            % Print result
            %transpose(Zp)
        end

        
        
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
            
            assert(isa(TabTf, 'EJ_library.utils.tabulated_transform'))
            assert(isfinite(valueOutsideTable))
            
            if 1
                % NOTE: interp1 returns NaN for values outside range.
                Z = interp1(TabTf.omegaRps, TabTf.Z, omegaRps, 'linear');
            else
                % EXPERIMENTAL
                Z = bicas.calib_utils.interpolate_TF(TabTf.omegaRps, TabTf.Z, omegaRps);
            end
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
        
        
        
        % Function for multiplying two TFs.
        %
        % RATIONALE
        % =========
        % In principle, this function is quite unnecessary for multiplying TFs,
        % but it is useful for putting breakpoints in when debugging TFs which
        % are built from layers of anonymous functions and function handles. It
        % is also therefore that this function does NOT return a function
        % pointer.
        function Z = multiply_TFs(omegaRps, tf1, tf2)
            % PROPOSAL: Re-purpose into function only for combining BIAS and
            % non-BIAS TFs.
            
            Z = tf1(omegaRps) ...
                .* ...
                tf2(omegaRps);
        end
        
        
        
        function log_TF_tabulated(logLevel, tfName, Tf, L)
            % PROPOSAL: Somehow prevent printing unnecessary trailing zeros.
            
            assert(isa(Tf, 'EJ_library.utils.tabulated_transform'))
            
            assert(numel(tfName) <= 38, ...
                'String argument "tfName" is too long. numel(tfName)=%i.', numel(tfName))
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
        function log_TF_function_handle(logLevel, tfName, tfUnit, freqHzArray, tfFuncHandle, L)
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
                %     row   68, execute_sw_mode.m,          execute_sw_mode
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
                    'String argument "tfName" is too long. numel(tfName)=%i.', numel(tfName))
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
            EJ_library.assert.vector(v)
            s = sprintf('(%s)', strjoin(EJ_library.str.sprintf_many(pattern, v), ', '));
        end
        
        
        
        function assert_iBlts(iBlts)
            assert(ismember(iBlts, [1:5]), ...
                'BICAS:calib:IllegalArgument:Assertion', ...
                'Illegal value iBlts=%g', iBlts)
        end
        
        
        
        function assert_iLsf(iLsf)
            assert(ismember(iLsf,  [1:4]), ...
                'BICAS:calib:IllegalArgument:Assertion', ...
                'Illegal value iLsf=%g.', iLsf)
        end
        
        
        
    end    % methods(Static)
    
    
    
end    % classdef calib_utils



            % TEST
%             if 0
%                 
%                 OMEGA_LIMIT_RPS = 100 * 2*pi;
%                 zLimit = itfIvpt(OMEGA_LIMIT_RPS);
%                 
%                 if isfinite(zLimit)
%                     omegaHz = logspace(log10(0.1), log10(1e3), 1e3);
%                     omegaRps = omegaHz * 2*pi;
%                     
%                     Z = itfIvpt(omegaRps);
%                     
%                     close all
%                     subplot(4,1, 1)
%                     loglog(omegaHz, abs(Z))
%                     subplot(4,1, 2)
%                     semilogx(omegaHz, angle(Z))
% 
%                     itfIvpt = @(omegaRps) (bicas.calib.TF_w_constant_abs_Z_LF(itfIvpt, omegaRps, OMEGA_LIMIT_RPS, zLimit));
%                     Z = itfIvpt(omegaRps);
%                     
%                     subplot(4,1, 3)
%                     loglog(omegaHz, abs(Z))
%                     subplot(4,1, 4)
%                     semilogx(omegaHz, angle(Z))
%                     
%                     false;
%                 end
%                 
%             end
            
