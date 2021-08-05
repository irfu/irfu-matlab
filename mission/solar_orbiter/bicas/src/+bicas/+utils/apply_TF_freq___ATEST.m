%
% Semi-automatic test code for function "bicas.utils.apply_TF_freq".
%
%
% TS = Time Series. Data in time domain.
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
    % TODO-NI: How robust is automatic test code in the event of better algorithms? (Hann windows, de-trending and stuff...)
    % PROPOSAL: Use specially-written "equals"-like function that permits (greater) discrepancies at the edge(s).
    %   PROPOSAL: Use weighting function.
    % PROPOSAL: Not only use delay TFs. It is possible that this makes the test insensitive to the function's association of
    %           frequencies with DFT components.
    %
    % TODO-NI: Should it not be possible to get more accurate agreement between expectation and result (without de-trending)?
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
    %   PRO: Can use functions everywhere, also in local functions.
    %
    % PROPOSAL: Use EJ_library.atest
    %   PRO: Better isolated test setups.
    
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
        y2_act = bicas.utils.apply_TF_freq( inputCa{iTest}{:} );
        
        n = (1:length(y2_exp))';
        y1 = inputCa{iTest}{2};
        
        %==============================================
        % Plot
        % -- Input           into apply_TF_time()
        % -- Actual output   from apply_TF_time()
        % -- Expected output from apply_TF_time()
        % NOTE: All are in time-domain.
        %==============================================
        figure
        plot(n, y1,     '-');   hold on
        plot(n, y2_exp, '-k', 'LineWidth', 2.0)
        plot(n, y2_act, '*')
        legend('y1', 'y2\_expected', 'y2\_result')
        xlabel('Array index (not t)')

        
        
        if ~EJ_library.utils.approx_equals(y2_act, y2_exp, EPSILON, 'NaN equal to itself')
            %error('BICAS:TEST', 'TEST FAILED')
            warning('TEST FAILED')
            
            [y2diffMax, iY2DiffMax] = max(abs(y2_act - y2_exp))
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
        % N  = 10;
        dt    = 0.1;
        delay = 8*dt;
        
        %tfOmega = linspace(0, 4/dt, 1e6)';
        tf = @(omega) tf_delay(omega, delay);
        
        t  = time_vector_N_dt(N, dt);
        f  = @(t) (0.5*(t-3).^2 + 5);
        y1 = f(t);
        y2 = ts_delay_func(f, t, delay, N, dt);
        
        inputCa{end+1} = {dt, y1, tf};
        outputCa{end+1} = y2;
    end
    
    if 1
        % Signal: Constant function
        % TF    : Constant Z != 1
        N  = 100;
        dt = 0.1;
        
        tf    = @(omega) constantTfZ(omega, 2, 2);
        
        t  = time_vector_N_dt(N, dt);
        f1 = @(t) (1 * ones(size(t)) );
        f2 = @(t) (2 * ones(size(t)) );
        y1 = f1(t);
        y2 = f2(t);
        
        inputCa{end+1} = {dt, y1, tf};
        outputCa{end+1} = y2;
    end
    
    if 1
        % Signal: Approximately one non-zero DFT component
        % TF:     Constant Z(!), except tfZ(omega=0). Different time delays on
        % different frequencies, which produces a chosen time delay for this
        % specific signal.
        for N = 100:101
            %     for N = 10
            delay  = 0.60;
            dt     = 2*pi/N;
            omega0 = 1;    % Fits time interval perfectly. Perfectly periodic.
            
            tf = @(omega) constantTfZ(omega, 1, exp(1i*omega0*(-delay)));
            t  = time_vector_N_dt(N, dt);
            
            f  = @(t) (3+cos(omega0 * (t-pi/5)));
            
            y1 = f(t);
            y2 = ts_delay_func(f, t, delay, N, dt);

            
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
            
            tf = @(omega) tf_delay(omega, delay);
            
            omega0 = 2*pi * 5/(N*dt);    % Exact DFT frequency.  ==> Good match
            %omega0 = 2*pi * 1/3;          % Arbitrary frequency.  ==> Edge effects, generally
            if (N*dt) / (1/omega0) < 10   % Assert minimum number of oscillation periods in function (in radians).
                error('BICAS:TEST', 'Bad test config.?')
            end
            
            t  = time_vector_N_dt(N, dt);
            f = @(t) cos(omega0 * (t-pi/5));
            
            y1 = f(t);
            y2 = ts_delay_func(f, t, delay, N, dt);

            
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
            
            tf = @(omega) tf_delay(omega, delay);
            
            t  = time_vector_N_dt(N, dt);
            
            f = @(t) exp(t);
            y1 = f(t);
            y2 = ts_delay_func(f, t, delay, N, dt);
            
            
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
            
            %tfOmega = linspace(0, 4*N/dt, 1e6)';
            %tf     = @(omega) delayTfZ(omega, delay);
            tf = @(omega) tf_delay(omega, delay);
            
            t  = time_vector_N_dt(N, dt);
            
            f = @(t) ((t-5).^3 - 20*t + 25);
            y1 = f(t);
            y2 = ts_delay_func(f, t, delay, N, dt);
            
            inputCa{end+1} = {dt, y1, tf};
            outputCa{end+1} = y2;
        end
    end
    
    
    if 1
        %inputCa  = {};
        %outputCa = {};
        [inputCa, outputCa] = test_Only_keep_constant(inputCa, outputCa);
    end
    
end



function [inputCa, outputCa] = test_Only_keep_constant(inputCa, outputCa)
    % TF for keeping only omega=0 component.
    
    N = 100;
    dt = 1 / (N-1);
    Z0 = -2;
    
    tf      = @(omega) tf_scale(omega, Z0);
    
    t = time_vector_N_dt(N, dt);
    %t_middle = mean(t);
    
    f_const    = @(t) (10 * ones(size(t)));   % Constant part.
    f_nonconst = @(t) (5 * sin(10*t.^3));        % Non-constant part.
    % IMPLEMENTATION NOTE: Must remove constant part (mean) of f2.
    y1 = f_const(t) + detrend(f_nonconst(t), 0);
    y2 = f_const(t) * Z0;
    
    inputCa{end+1} = {dt, y1, tf};
    outputCa{end+1} = y2;
end



% Create TF that delays function, i.e. it is moved in the t+ direction (time
% domain), i.e. the same time delay (in time units) for all frequencies.
function Z = tf_delay(omega, delay)
%delayTfZ    = @(omega, delay)   (exp(1i*omega*(-delay)));

    Z = exp(1i*omega*(-delay));
end



% Create TF that is zero except for omega=0.
%
% Should represent keeping only the mean and multiplying it by Z0 (real).
%
function Z = tf_scale(omega, Z0)
    assert(isreal(Z0))
    
    Z = (omega == 0 ) * Z0;
end



% Function for creating TIME SERIES which is a delayed version of a specified
% function (not time series) AND treats the time interval as cyclic. 
%
% f     : Function pointer
% t     : Array of timestamps.
% delay : delay>0 pushes function/data/shapes in the t+ direction.
% N     : Number of samples in time series.
% dt    :
%
function y2 = ts_delay_func(funcHandle,t, delay,N,dt)
    % PROPOSAL: Replace by function that has argument for y data instead.
    %   ~PROBLEM: Must determine delay in units of dt. => Rounding.
    %       ==> Can not be used for approximate results (more manual atests).
    y2 = funcHandle(mod(t-delay, N*dt));    
    
    
%     delayedFunc = @(f,t,delay,N,dt) (f(mod(t-delay, N*dt)));
end
% function y2 = ts_delay(y1, delay, dt)
%     N = numel(y1);
%     y2 = y1( mod(1:N-delay, N*dt) );
% end



% Function for creating timestamps vector.
function [t, t2] = time_vector_N_dt(N,dt)
    % tVec = @(N,dt) (0 : dt : ((N-1)*dt) )';
    
    t2 = (N-1)*dt;
    t  = [0 : dt : t2]';
end
