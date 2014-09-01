function [spinnedDataProbePair] = mms_add_fake_spin(dataProbePair, phase )
%MMS_ADD_FAKE_SPIN add a sinusodial variation to input and return.
% as the MMS values are not at all periodic in the test data provided use
% this function to add a fake sun periodic value to the signal. This signal
% can then be identified and removed just as the proper code should do in
% production of real data.

narginchk(2,2);

A = 0.10; % Scale amplitude of Sinusodial wave.
d = 0.00; % Shift phase of Sinusodial wave.

% Simply add A*cos(phase + d) to signal.
spinnedDataProbePair = dataProbePair + A * cosd(phase + d);

end

