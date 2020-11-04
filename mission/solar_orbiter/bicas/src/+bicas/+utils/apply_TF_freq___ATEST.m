%
% Semi-automatic test code for function "bicas.utils.apply_TF_freq".
%
%
% Author: Erik P G Johansson, IRF, Uppsala, Sweden
% First created 2017-02-13
%
function apply_TF_freq___ATEST
    %
    % NOTE: Want to test even & odd N, dt<>1.
    
    % BOGIQ:
    % ------
    % PROPOSAL: Utility Z(omega) for generating single frequency change.
    % PROPOSAL: Automatic test: tf+inversion, compare ~spectras? Same result for multiple resolutions (even/odd). Pure delay/advance.
    % TODO-NEED-INFO: How robust is automatic test code in the event of better algorithms? (Hann windows, de-trending and stuff...)
    % PROPOSAL: Use specially-written "equals"-like function that permits (greater) discrepancies at the edge(s).
    %   PROPOSAL: Use weighting function.
    % PROPOSAL: Not only use delay TFs. It is possible that this makes the test insensitive to the function's association of
    %           frequencies with DFT components.
    %
    % TODO-NEED-INFO: Should it not be possible to get more accurate agreement between expectation and result (without de-trending)?
    %       Gets up to y2diffMax = 0.0889.
    %       Compare nDelay=3; y1 = rand(1,100); N=numel(y1); yDft1 = fft(y1); yDft2 = yDft1 .* exp(1i*2*pi*(0:N-1)* -nDelay/N); y2 = ifft(yDft2); y2Pred = circshift(y1, [1,nDelay]); transpose([y2; y2Pred]), max(abs(y2-y2Pred))
    %       which produces maxDiff < 1e-15.
    %   PROPOSAL: Only set up examples using no "time t" and only array indices.
    %       PRO: Increase likelyhood of exact predictions.
    %       CON: Harder to change only the resolution, but keep the curve.
    %   PROPOSAL: Use circshift, instead of mod to cycle/delay functions/arrays.
    %   PROPOSAL: Use random input arrays.
    %
    % PROPOSAL: Test scaling X_0 and X_k, k>0 differently.
    %
    % PROPOSAL: Proper functions (with longer descriptive names) for creating functions, instead of function pointers.
    
    % EPSILON = 1e-6;
    EPSILON = 1e-4;
    
    
    
    [inputCa, outputCa] = define_tests();
%     inputCa  = inputCa(1);
%     outputCa = outputCa(1);

    %======================
    % Run & plot all tests
    %======================
    close all
    for iTest = 1:numel(inputCa)
        
        y2_exp = outputCa{iTest};
        y2_res = bicas.utils.apply_TF_freq( inputCa{iTest}{:} );
        
        n = (1:length(y2_exp))';
        y1 = inputCa{iTest}{2};
        figure
        plot(n, y1,     '-');   hold on
        plot(n, y2_exp, '-k', 'LineWidth', 2.0)
        plot(n, y2_res, '*')
        legend('y1', 'y2\_expected', 'y2\_result')
        xlabel('Array index (not t)')
        
        if ~EJ_library.utils.approx_equals(y2_res, y2_exp, EPSILON, 'NaN equal to itself')
            %error('BICAS:TEST', 'TEST FAILED')
            warning('TEST FAILED')
            
            [y2diffMax, iY2DiffMax] = max(abs(y2_res - y2_exp))
            %         keyboard
            %         close all
        else
            disp('TEST OK')
        end
    end
    
end



function [inputCa, outputCa] = define_tests()
    inputCa  = {};
    outputCa = {};
    
    % Function for creating time samples vector.
    tVec = @(N,dt) (0 : dt : ((N-1)*dt) )';
    
    %===========================================================================
    % Function for creating another function (time domain) which is a delayed
    % version of a specified function AND treats it as cyclic. delay>0 pushes it
    % in the t+ direction.
    %
    % f  : Function pointer
    % N  : Number of samples in time series.
    % dt :
    %===========================================================================
    delayedFunc = @(f,t,delay,N,dt) (f(mod(t-delay, N*dt)));
    
    %===========================================================================
    % TF that delays function, i.e. it is moved in the t+ direction (time
    % domain), i.e. one time delay for all frequencies
    %===========================================================================
    delayTfZ    = @(omega, delay)   (exp(1i*omega*(-delay)));
    
    %===========================================================================
    % TF for (almost) constant Z. Implements Z(omega=0)=z0 and Z(omega>0)=z1.
    %
    % NOTE: Must work for omega=vector.
    % NOTE: z0 should be real.
    %===========================================================================
    constantTfZ = @(omega, z0, z1) ( (omega==0)*z0 + (omega~=0)*z1 );
    
    if 1
        % Signal: Quadratic function
        % TF    : Constant delay.
        N  = 100;
        %     N  = 10;
        dt = 0.1;
        delay = 8*dt;
        
        %tfOmega = linspace(0, 4/dt, 1e6)';
        tf     = @(omega) delayTfZ(omega, delay);
        
        t  = tVec(N, dt);
        f = @(t) (0.5*(t-3).^2 + 5);
        y1 = f(t);
        y2 = delayedFunc(f, t, delay, N, dt);
        
        inputCa{end+1} = {dt, y1, tf};
        outputCa{end+1} = y2;
    end
    
    if 1
        % Signal: Constant function
        % TF    : Constant Z != 1
        N  = 100;
        dt = 0.1;
        
        tf    = @(omega) constantTfZ(omega, 2, 2);
        
        t  = tVec(N, dt);
        f1 = @(t) (1 * ones(size(t)) );
        f2 = @(t) (2 * ones(size(t)) );
        y1 = f1(t);
        y2 = f2(t);
        
        inputCa{end+1} = {dt, y1, tf};
        outputCa{end+1} = y2;
    end
    
    if 1
        % Signal: Approximately one non-zero DFT component
        % TF:     Constant Z(!), except tfZ(omega=0). Different time delays on different frequencies, which produces a chosen time delay for this specific signal.
        for N = 100:101
            %     for N = 10
            delay  = 0.60;
            dt     = 2*pi/N;
            omega0 = 1;    % Fits time interval perfectly. Perfectly periodic.
            
            tf      = @(omega) constantTfZ(omega, 1, exp(1i*omega0*(-delay)));
            t  = tVec(N, dt);
            f  = @(t) (3+cos(omega0 * (t-pi/5)));
            
            y1 = f(t);
            y2 = delayedFunc(f, t, delay, N, dt);
            
            inputCa{end+1} = {dt, y1, tf};
            outputCa{end+1} = y2;
        end
    end
    
    if 1
        % Signal: Single non-zero DFT component
        % TF:     One delay for all frequencies.
        for N = 100:101
            %     for N = 10
            % TF for a delay.
            dt    = 0.1;
            delay = 3*dt;
            
            tf    = @(omega) delayTfZ(omega, delay);
            
            omega0 = 2*pi * 5/(N*dt);    % Exact DFT frequency.  ==> Good match
            %omega0 = 2*pi * 1/3;          % Arbitrary frequency.  ==> Edge effects, generally
            if (N*dt) / (1/omega0) < 10   % Assert minimum number of oscillation periods in function (in radians).
                error('BICAS:TEST', 'Bad test config.?')
            end
            
            t  = tVec(N, dt);
            f = @(t) cos(omega0 * (t-pi/5));
            
            y1 = f(t);
            y2 = delayedFunc(f,t,delay,N,dt);
            
            inputCa{end+1} = {dt, y1, tf};
            outputCa{end+1} = y2;
        end
    end
    
    if 1
        % Signal: Arbitrary function + delay
        %     for N = 100
        for N = 100:101
            % TF for a delay.
            dt    = 0.01;
            delay = 13*dt;
            
            tf    = @(omega) delayTfZ(omega, delay);
            
            t = tVec(N, dt);
            
            f = @(t) exp(t);
            y1 = f(t);
            y2 = delayedFunc(f,t,delay,N,dt);
            
            inputCa{end+1} = {dt, y1, tf};
            outputCa{end+1} = y2;
        end
    end
    
    if 1
        % Arbitrary function + delay
        for N = 100:101
            % TF for a delay.
            dt    = 1 / (N-1);
            delay = 10*dt;
            
            tfOmega = linspace(0, 4*N/dt, 1e6)';
            tf     = @(omega) delayTfZ(omega, delay);
            
            t = tVec(N,dt);
            
            f = @(t) ((t-5).^3 - 20*t + 25);
            y1 = f(t);
            y2 = delayedFunc(f,t,delay,N,dt);
            %         y1(3) = NaN;    % TEST
            
            inputCa{end+1} = {dt, y1, tf};
            outputCa{end+1} = y2;
        end
    end
end
