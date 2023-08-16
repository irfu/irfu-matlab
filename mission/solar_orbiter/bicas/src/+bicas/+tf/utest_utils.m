%
% Collection of utility functions used by automatic test code for the bicas.tf
% package.
%
%
% Author: Erik P G Johansson, Uppsala, Sweden
% First created 2021-09-23
%
classdef utest_utils
    % PROPOSAL: Shorter name.
    
    %#######################
    %#######################
    % PUBLIC STATIC METHODS
    %#######################
    %#######################
    methods(Static)
        
        
        
        % Return TF (function handle) that delays function, i.e. it is moved in
        % the t+ direction (time domain), i.e. the same time delay (in time
        % units) for all frequencies.
        function tf = get_TF_delay(delaySec)
            assert(isscalar(delaySec))
        
            tf = @(omegaRps) (exp(1i*omegaRps*(-delaySec)));            
        end
        
        
        
        % Return TF (function handle) for TF with (almost) constant Z, where
        % Z(omega=0)=z0 and Z(omega~=0)=z1.
        %
        % z1 == 0: 
        %   Keeping only the mean and multiplying it by z0 (real).
        %
        % z0 == z1:
        %   Constant TF. Multiply signal by z0.
        %
        % z0 == 0:
        %   Removing the mean and multiply the remainder by z1???
        %
        %
        % NOTE: z0 should be real to maintain that real input signal ==> real
        % output signal.
        %
        function tf = get_TF_constant(z0, z1)
            assert(isscalar(z0) & isreal(z0))
            assert(isscalar(z1))
            
            tf = @(omegaRps) ((omegaRps==0)*z0 + (omegaRps~=0)*z1);
        end
        
        
        
    end    % methods(Static)
    
    
    
end
