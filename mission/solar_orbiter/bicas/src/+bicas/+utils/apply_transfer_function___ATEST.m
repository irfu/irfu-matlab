% Semi-automated test code for function "apply_transfer_function".
%
%
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-16
%
function apply_transfer_function___ATEST
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-13
%
% NOTE: Want to test even & odd N, dt<>1.
%
% UNFINISHED.

% BOGIQ:
% ------
% PROPOSAL: Function for generating delay/advance TF, single frequency change.
% PROPOSAL: Inner function that generates TF for delay. Function pointer to function?
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



% EPSILON = 1e-6;
EPSILON = 1e-3;

input  = {};
output = {};

tVec = @(N,dt) (0 : dt : ((N-1)*dt) )';
delayedFunc = @(f,t,delay,N,dt) (f(mod(t-delay, (N-dt)*dt)));    % Delays function (in time domain), AND treats it as cyclic. delay>0 pushes it in t+ direction.
delayTfZ    = @(omega, delay)   (exp(1i*omega*(-delay)));        % TF that delays function (in time domain), i.e. one time delay for all frequencies.



if 1
    % Signal: Quadratic function
    % TF    : Constant delay.
    N  = 100;
    dt = 0.1;
    delay = 12*dt;
    
    tfOmega = linspace(1e-9, 4/dt, 1e6)';
    tfZ     = delayTfZ(tfOmega, delay);
    
    t  = tVec(N, dt);
%     f = @(t) (0*t + 2);
%     f = @(t) (1*t + 0);
    f = @(t) (0.5*t.^2 - 1*t + 2);
    y1 = f(t);
    y2 = delayedFunc(f, t, delay, N, dt);

    input{end+1} = {dt, y1, tfOmega, tfZ, 'enableDetrending', 0};
    output{end+1} = y2;
end

if 1
    % Signal: Approximately one non-zero DFT component
    % TF:     Constant Z(!) Different time delays on different frequencies, which produces a chosen time delay for this specific signal.
    for N = 100:101
        delay  = 0.60;
        dt     = 2*pi/N;
        omega0 = 1;    % Fits time interval perfectly. Perfectly periodic.
                
        tfOmega = [1e-8, 100]';   % Deliberately low-resolved.
        tfZ     = ones(size(tfOmega)) * exp(1i*omega0*(-delay));
        t  = tVec(N, dt);
        f  = @(t) (3+cos(omega0 * (t-pi/5)));
        
        y1 = f(t);
        y2 = delayedFunc(f, t, delay, N, dt);
        
        input{end+1} = {dt, y1, tfOmega, tfZ, 'enableDetrending', 0};
        output{end+1} = y2;
    end
end

if 1
    % Signal: Single non-zero DFT component
    % TF:     One delay for all frequencies.
    for N = 100:101
        % TF for a delay.
        dt    = 0.1;
        delay = 3*dt;
        
        tfOmega = linspace(1e-9, 4/dt, 1e6)';
        tfZ     = delayTfZ(tfOmega, delay);
        
        omega0 = 2*pi * 12/(N*dt);   % Exact DFT frequency.  ==> Good match
        %omega0 = 2*pi * 1/3;         % Arbitrary frequency.  ==> Edge effects, generally
        if (N*dt) / (1/omega0) < 10   % Assert minimum number of oscillation periods in function (in radians).
            error('Bad test config.?')
        end
        
        t  = tVec(N, dt);
        f = @(t) cos(omega0 * (t-pi/5));
        
        y1 = f(t);
        y2 = delayedFunc(f,t,delay,N,dt);
        
        input{end+1} = {dt, y1, tfOmega, tfZ, 'enableDetrending', 0};
        output{end+1} = y2;
    end
end

if 1
    % Signal: Arbitrary function + delay
    for N = 100
%     for N = 100:101
        % TF for a delay.
        dt    = 0.01;
        delay = 13*dt;
        
        tfOmega = linspace(1e-8, 4/dt, 1e6)';
        tfZ = delayTfZ(tfOmega, delay);
        
        t = tVec(N, dt);
        
        f = @(t) exp(t);
        y1 = f(t);
        y2 = delayedFunc(f,t,delay,N,dt);
        
%         input{end+1} = {dt, y1, tfOmega, tfZ, 'enableDetrending', 1};    % De-trend
        input{end+1} = {dt, y1, tfOmega, tfZ, 'enableDetrending', 0};    % De-trend
        output{end+1} = y2;
    end
end

if 1
    % Arbitrary function + delay
    for N = 100:101
        % TF for a delay.
        dt = 1 / (N-1);
        delay = 10*dt;
        
        tfOmega = linspace(1e-8, 4*N/dt, 1e6)';
        tfZ     = delayTfZ(tfOmega, delay);
        
        t = tVec(N,dt);
        
        f = @(t) ((t-5).^3 - 20*t + 25);
        y1 = f(t);
        y2 = delayedFunc(f,t,delay,N,dt);
%         y1(3) = NaN;    % TEST
        
        input{end+1} = {dt, y1, tfOmega, tfZ, 'enableDetrending', 0};
%         input{end+1} = {dt, y1, tfOmega, tfZ, 'enableDetrending', 1};
        output{end+1} = y2;
    end
end



close all
for iTest = 1:length(input)
    
    y2_exp = output{iTest};
    y2_res = bicas.utils.apply_transfer_function( input{iTest}{:} );
    
    n = (1:length(y2_exp))';
    y1 = input{iTest}{2};
    figure
%     plot(n, [y1, y2_exp])
    plot(n, [y1, y2_exp, y2_res])
    legend('y1', 'y2\_expected', 'y2\_result')
    xlabel('array index (not t)')
    
    if ~EJ_library.utils.approx_equals(y2_res, y2_exp, EPSILON, 'NaN equal to itself')
        %error('TEST FAILED')
        warning('TEST FAILED')
        
        y2diffMax = max(abs(y2_res - y2_exp))
%         keyboard
%         close all
    else
        disp('TEST OK')
    end
end

end
