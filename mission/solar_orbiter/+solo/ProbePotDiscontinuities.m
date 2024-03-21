function discontTimes=ProbePotDiscontinuities()
%SOLO.ProbePotDiscontinuities() contains a list of times where there are
%discontinuities in the probe-to-spacecraft potentials which may affect
%calibration.
% Input: N/A
%
% Output: discontTimes - EpochTT object containing times of solar orbiter
% probe-to-spacecraft potential discontinuities.
%
% Example: discontTimes=solo.ProbePotDiscontinuities()
%
% NOTE: This function is used by BICAS for producing calibration-files
% used in official datasets.

discontTimes=EpochTT(['2020-12-21T01:39:23.246487Z';... %Potential jump in V3 related to solar panel currents
  '2021-01-10T23:33:05.376840Z';... %Potential jump in V2 related to solar panel currents
  '2021-01-17T00:19:26.999006Z';... %Potential jump in all probes related to solar panel currents
  '2021-03-01T00:45:20.392809Z';... %Potential jump in V? related to solar panel currents
  '2022-02-22T18:37:29.476500Z']);  %Potential jump in all probes, mainly V3.

end
