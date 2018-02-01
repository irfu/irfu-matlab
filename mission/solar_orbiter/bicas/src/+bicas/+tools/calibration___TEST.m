% Semi-automated test code for "calibration".
% 
% Author: Erik P G Johansson, IRF-U, Uppsala, Sweden
% First created 2017-02-16
%
function calibration_TEST
    apply_transfer_function_in_freq_TEST
end



function apply_transfer_function_in_freq_TEST
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
% QUESTION: How robust is automatic test code in the event of better algorithms? (Hann windows, de-trending and stuff...)
% PROPOSAL: Use specially-written "equals"-like function that permits (greater) discrepancies at the edge(s).

EPSILON = 1e-6;

input = {};
output = {};

% Unit TF.
%tfOmega = [1e-8, 100]';
%tfZ     = [1,      1]';
%y1 = ones(10,1);
%input{end+1} = {0.123456, y1, tfOmega, tfZ, 0};
%output{end+1} = y1;

if 1
    % Single non-zero DFT component
    for N = 10:11
        tfOmega = [1e-8, 100]';
        tfZ     = [-1i,  -1i].';
        dt = 0.3;
        t = (0 : dt : ((N-1)*dt) )';
        y1 = cos(2*pi * (1/(N*dt) * t       ));
        y2 = cos(2*pi * (1/(N*dt) * t - 0.25));
        
        input{end+1} = {dt, y1, tfOmega, tfZ, 'EnableDetrending', 0};
        output{end+1} = y2;
    end
end

if 1
    % Single non-zero DFT component + delay
    for N = 100:101
        % TF for a delay.
        delay = 0.25;
        tfOmega = linspace(1e-8, 200, 1e6)';
        tfZ     = exp(-1i* tfOmega.*delay);
        
        dt = 0.1;
        omegaY = 2*pi * 12/(N*dt);   % Exact DFT frequency.  ==> Good match
        %omegaY = 2*pi * 1/3;         % Arbitrary frequency.  ==> Edge effects, generally
        if (N*dt) / (1/omegaY) < 10   % Assert minimum number of oscillation periods in function (in radians).
            error('Bad test config.?')
        end
        t = (0 : dt : ((N-1)*dt) )';
        
        f = @(t) cos(omegaY * t);
        y1 = f(t);
        y2 = f(t - delay);
        
        input{end+1} = {dt, y1, tfOmega, tfZ, 'EnableDetrending', 0};
        output{end+1} = y2;
    end
end

if 1
    % Arbitrary function + delay
    for N = 100:101
        % TF for a delay.
        delay = 0.25;
        tfOmega = linspace(1e-8, 200, 1e6)';
        tfZ     = exp(-1i* tfOmega.*delay);
        
        dt = 0.1;
        %omegaY = 2*pi * 12/(N*dt);   % Exact DFT frequency.  ==> Good match
        %omegaY = 2*pi * 1/3;         % Arbitrary frequency.  ==> Edge effects, generally
        t = (0 : dt : ((N-1)*dt) )';
        
        f = @(t) exp(2 * t / (N*dt));
        y1 = f(t);
        y2 = f(t - delay);
        %y1(3) = Inf;    % TEST
        
        input{end+1} = {dt, y1, tfOmega, tfZ, 'EnableDetrending', 1};    % De-trend
        output{end+1} = y2;
    end
end

if 1
    % Arbitrary function + delay
    for N = 100:101
        % TF for a delay.
        delay = 0.25;
        tfOmega = linspace(1e-8, 200, 1e6)';
        tfZ     = exp(-1i* tfOmega.*delay);
        
        %omegaY = 2*pi * 12/(N*dt);   % Exact DFT frequency.  ==> Good match
        %omegaY = 2*pi * 1/3;         % Arbitrary frequency.  ==> Edge effects, generally
        t = linspace(-5, 5, N)';
        dt = (max(t) - min(t)) / (N-1);
        
        f = @(t) (t.^3 - 20*t);
        y1 = f(t);
        y2 = f(t - delay);
        %y1(3) = Inf;    % TEST
        
        input{end+1} = {dt, y1, tfOmega, tfZ, 'EnableDetrending', 1};    % De-trend
        output{end+1} = y2;
    end
end



for iTest = 1:length(input)
    
    y_exp = output{iTest};
    y_res = bicas.calibration.apply_transfer_function_in_freq( input{iTest}{:} );
    
    if ~bicas.utils.equals_tolerance(y_res, y_exp, EPSILON)
        %error('TEST FAILED')
        warning('TEST FAILED')
        
        %            yDeviationMax = max(abs(y_res - y_exp))
        
        close all
        n = (1:length(y_exp))';
        y1 = input{iTest}{2};
        plot(n, [y1, y_exp, y_res])
        legend('y1', 'y\_expected', 'y\_result')
        keyboard
    end
    disp('TEST OK')
end

end
