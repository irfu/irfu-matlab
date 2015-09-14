function MMS_CONST = mms_constants
%MMS_CONST  initialize MMS constants
%
%  MMS_CONST = mms_constants()


MMS_CONST.MinFileVer = 0; % min version of l1b files accepted

MMS_CONST.MMSids = 1:4;

% Spin rate max and min, nominally 3.05 rpm +/-0.05 rpm.
% After June 18, 2015
MMS_CONST.Spinrate.max = 3.1; % Rev per Minute.
MMS_CONST.Spinrate.min = 3.0; % Rev per Minute.
% Spin rate nominally 3.1 rpm +/-0.1 rpm.
% After May 1, 2015
MMS_CONST.Spinrate.max_comm = 3.2; % Rev per Minute.
% Spin rate max 7.3 rpm +/-0.1 rpm during thin wire deployments
% Until May 1, 2015
MMS_CONST.Spinrate.max_deploy = 7.4; % Rev per Minute.

% Spin phase is zero when BSC x-axis points sunward, and
% increases monotonically. As the spacecraft spins, the order in
% which the probes and BCS axes point sunward is the following:
% at   0 degrees, BSC +X is sunward
% at  60 degrees, probe 4 is sunward
% at  90 degrees, BSC -Y is sunward
% at 150 degrees, probe 2 is sunward
% at 180 degrees, BSC -X is sunward
% at 240 degrees, probe 3 is sunward
% at 270 degrees, BSC +Y is sunward
% at 330 degrees, probe 1 is sunward
%
% Angles when phase=0 (X BSC direction)
MMS_CONST.Phaseshift.e12 =  2*pi*150/360; % probe 2 sunward
MMS_CONST.Phaseshift.e34 =  2*pi* 60/360; % probe 4 sunward

% Nominal Amplitude Correction factor multiplied to DCE data.
MMS_CONST.NominalAmpCorr.e12 = 1.1;
MMS_CONST.NominalAmpCorr.e34 = 1.1;
MMS_CONST.NominalAmpCorr.e56 = 1.5;

% Telemetry mode
MMS_CONST.TmModes = {'slow','fast','brst', 'comm'};
MMS_CONST.TmMode.slow = 1; % Number must corrspond to position in the list
MMS_CONST.TmMode.fast = 2;
MMS_CONST.TmMode.brst = 3;
MMS_CONST.TmMode.comm = 4; % Commissioning data.

% Sample rates
MMS_CONST.Samplerate.slow = 8; % Samples per second (dce & dcv), TM mode slow
MMS_CONST.Samplerate.fast = 32; % TM mode fast
MMS_CONST.Samplerate.comm_8 = 8; % Commissioning "Slow"
MMS_CONST.Samplerate.comm_32 = 32; % Commissioning "I&T" phase
MMS_CONST.Samplerate.comm_64 = 64; % Commissioning "Turn ON" phase
MMS_CONST.Samplerate.comm_128 = 128; % Commissioning "Boom deployment" phase
MMS_CONST.Samplerate.brst = 8192; % Or 1024? TM mode burst

% SDC process names
MMS_CONST.SDCProcs = {'sitl','ql','scpot','l2pre','l2a'};
MMS_CONST.SDCProc.sitl  = 1; % Number must corrspond to position in the list
MMS_CONST.SDCProc.ql    = 2;
MMS_CONST.SDCProc.scpot = 3;
MMS_CONST.SDCProc.l2pre = 4;
MMS_CONST.SDCProc.l2a   = 5;

% Limits used in processing
MMS_CONST.Limit.LOW_DENSITY_SATURATION = -100; % Probe stuck and below limit.
MMS_CONST.Limit.DIFF_PROBE_TO_SCPOT_MEDIAN = 1.5; % Probe not used for probe2scpot if moving average is off by this from the mean of all probes moving average, in V.
MMS_CONST.Limit.DCE_DCV_DISCREPANCY = 0.28; % Max discrepancy DCE{12,34}=(DCV{1,3}-DCV{2,4})/NominalLength, for data with all probes.


% Bitmask values; 2^(bit_number - 1):
MMS_CONST.Bitmask.SIGNAL_OFF               =  uint16(1);       % Bit 1
MMS_CONST.Bitmask.BAD_BIAS                 =  uint16(2);       % Bit 2
MMS_CONST.Bitmask.PROBE_SATURATION         =  uint16(4);       % Bit 3
MMS_CONST.Bitmask.LOW_DENSITY_SATURATION   =  uint16(8);       % Bit 4
MMS_CONST.Bitmask.SWEEP_DATA               =  uint16(16);      % Bit 5
MMS_CONST.Bitmask.ADP_SHADOW               =  uint16(32);      % Bit 6

MMS_CONST.Error = -Inf; % Indicates error in computation

% Version numbering, start with X, Y, Z = 0, 0, 0. When releasing new
% software (with significant changes) update values here and subsequent
% output files created will have these numbers.
% Major new Software version, X
% New Calibration version, Y
MMS_CONST.Version = struct(...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.sitl},  struct('X', 0, 'Y', 3), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.ql},    struct('X', 0, 'Y', 3), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.scpot}, struct('X', 0, 'Y', 2), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.l2pre}, struct('X', 0, 'Y', 2), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.l2a},   struct('X', 0, 'Y', 1));
% Version Notes Y, for us. Not written to CDF files.
% Y - 0, initial release
% Y - 1, ADP shadow removal, use of irf_filt in scpot calculation.
% Y - 2 (scpot), apply_nom_amp_corr to DCE after all DCV was calculated.
MMS_CONST.Version.Z = 0; % File revision, increased by 1 for each re-run.
% Version.MODS - MODS cdf GlobalAttribute should contain a description of
% all significant changes to the data set, essentially capturing a log of
% high-level release notes. Can have as many entries as necessary and
% should be updated if the "X" value of the version number changes.
% Each cell corresponds to one version, append like: mods=[mods; {'new text'}];
MMS_CONST.Version.MODS = {'V.0. Initial release.'};

end
