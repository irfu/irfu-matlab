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
% (Ref: e-mail forwarded from PAL on 2015-03-26T15:53, originally from I.Dors, mentions a "29.85 deg", i.e. -30 for p1)
%
% Angles when phase=0 (X BSC direction)
MMS_CONST.Phaseshift.e12 =  2*pi*150/360; % probe 2 sunward
MMS_CONST.Phaseshift.e34 =  2*pi* 60/360; % probe 4 sunward
MMS_CONST.Phaseshift.p1  = 2*pi* 330/360; % probe 1 sunward
MMS_CONST.Phaseshift.p2  = 2*pi* 150/360; % probe 2 sunward
MMS_CONST.Phaseshift.p3  = 2*pi* 240/360; % probe 3 sunward
MMS_CONST.Phaseshift.p4  = 2*pi*  60/360; % probe 4 sunward
MMS_CONST.Phaseshift.dfg = 2*pi* 315/360; % DFG sunward
MMS_CONST.Phaseshift.afg = 2*pi* 135/360; % AFG sunward

% Nominal Amplitude Correction factor multiplied to DCE data.
MMS_CONST.NominalAmpCorr.e12 = 1.25;
MMS_CONST.NominalAmpCorr.e34 = 1.25;
MMS_CONST.NominalAmpCorr.e56 = 1.0;

% Inner Magnetosphere Amplitude Correction factor multiplied to DCE data.
MMS_CONST.InnerMSPAmpCorr.e12 = 1.0;
MMS_CONST.InnerMSPAmpCorr.e34 = 1.0;
MMS_CONST.InnerMSPAmpCorr.e56 = 1.0;

MMS_CONST.InnerMSPradius = 5*6.3712e+03; % 5 R_E in km.

% Telemetry mode
MMS_CONST.TmModes = {'slow','fast','brst', 'comm'};
MMS_CONST.TmMode.slow = 1; % Number must corrspond to position in the list
MMS_CONST.TmMode.fast = 2;
MMS_CONST.TmMode.brst = 3;
MMS_CONST.TmMode.comm = 4; % Commissioning data.

% Sample rates
MMS_CONST.Samplerate.slow = {8; 16; 32}; % Samples per second (dce & dcv), TM mode slow.
MMS_CONST.Samplerate.fast = 32; % TM mode fast
MMS_CONST.Samplerate.comm_8 = 8; % Commissioning "Slow"
MMS_CONST.Samplerate.comm_32 = 32; % Commissioning "I&T" phase
MMS_CONST.Samplerate.comm_64 = 64; % Commissioning "Turn ON" phase
MMS_CONST.Samplerate.comm_128 = 128; % Commissioning "Boom deployment" phase
MMS_CONST.Samplerate.brst = {8192; 1024; 16384}; % 8192 in tail, 1024 in sub-solar TM mode burst, special 16384 Hz at PI discretion.


% SDC process names
MMS_CONST.SDCProcs = {'ql','scpot','l2pre','l2a','l1ace','l2ace'};
MMS_CONST.SDCProc.ql    = 1; % Number must corrspond to position in the list
MMS_CONST.SDCProc.scpot = 2;
MMS_CONST.SDCProc.l2pre = 3;
MMS_CONST.SDCProc.l2a   = 4;
MMS_CONST.SDCProc.l1ace = 5; % ACE data (type of brst)
MMS_CONST.SDCProc.l2ace = 6;

% Limits used in processing
MMS_CONST.Limit.LOW_DENSITY_SATURATION = -100; % Probe stuck and below limit.
MMS_CONST.Limit.DCE_DCV_DISCREPANCY = 0.28; % Max discrepancy DCE{12,34}=(DCV{1,3}-DCV{2,4})/NominalLength, for data with all probes.
MMS_CONST.Limit.SPINFIT_INTERV = int64(20*10^9); % Perform spinfits covering this interval, in [ns].
MMS_CONST.Limit.MERGE_FREQ = 600; % Merging frequency (Hz) of measured and reconstructed electric fields

% Bitmask values; 2^(bit_number - 1):
MMS_CONST.Bitmask.SIGNAL_OFF             = mms_sdp_typecast('bitmask',1);  % Bit 1
MMS_CONST.Bitmask.BAD_BIAS               = mms_sdp_typecast('bitmask',2);  % Bit 2
MMS_CONST.Bitmask.PROBE_SATURATION       = mms_sdp_typecast('bitmask',4);  % Bit 3
MMS_CONST.Bitmask.LOW_DENSITY_SATURATION = mms_sdp_typecast('bitmask',8);  % Bit 4
MMS_CONST.Bitmask.SWEEP_DATA             = mms_sdp_typecast('bitmask',16); % Bit 5
MMS_CONST.Bitmask.ADP_SHADOW             = mms_sdp_typecast('bitmask',32); % Bit 6
MMS_CONST.Bitmask.ASPOC_RUNNING          = mms_sdp_typecast('bitmask',64); % Bit 7
% MMS_CONST.Bitmask.EDI_CORRECTION = mms_sdp_typecast('bitmask', 128); % Bit 8.
MMS_CONST.Bitmask.ASYMM_CONF             = mms_sdp_typecast('bitmask',256); % Bit 9
MMS_CONST.Bitmask.MANEUVERS              = mms_sdp_typecast('bitmask',512); % Bit 10
MMS_CONST.Bitmask.SW_WAKE_REMOVED        = mms_sdp_typecast('bitmask',1024); % Bit 11
MMS_CONST.Bitmask.ECLIPSE                = mms_sdp_typecast('bitmask',2048); % Bit 12

MMS_CONST.Error = -Inf; % Indicates error in computation

% Version numbering, start with X, Y, Z = 0, 0, 0. When releasing new
% software (with significant changes) update values here and subsequent
% output files created will have these numbers.
% Major new Software version, X
% New Calibration version, Y
MMS_CONST.Version = struct(...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.ql},    struct('X', 1, 'Y', 11), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.scpot}, struct('X', 2, 'Y', 9), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.l2pre}, struct('X', 2, 'Y', 2), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.l2a},   struct('X', 3, 'Y', 3), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.l1ace}, struct('X', 0, 'Y', 0), ...
  MMS_CONST.SDCProcs{MMS_CONST.SDCProc.l2ace}, struct('X', 0, 'Y', 0));
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
MMS_CONST.Version.MODS = [MMS_CONST.Version.MODS; {'V.1. QL (v1.0.z), SCPOT (v1.0.z), L2A (v0.1.z) now uses ASPOC srvy l2 and DEFATT, if these are available. Brst QL uses intermediate L2A file from Fast mode for delta offsets. Bitmask changed to uint16 and Quality to uint8.'}];
MMS_CONST.Version.MODS = [MMS_CONST.Version.MODS; {'V.2. SCPOT (v2.0.z), L2A (v1.0.z) now uses variable names in accordance with new recommended standard for FIELDS, All products change shortening factor to 1.25 on SDP, offsets applied indicated by GlobalAttribute Calibration_file.'}];
MMS_CONST.Version.MODS = [MMS_CONST.Version.MODS; {'V.2. L2a (v2.0.z), QL (v1.6.z) now try to remove solar wind wake which previously left a clear sinusodial signal in the data.'}];
MMS_CONST.Version.MODS = [MMS_CONST.Version.MODS; {'V.3. L2a (v3.0.z) Slow Mode probe Gain set to 1.0 when orbital radius less than 5 RE (1.25 otherwise), L2pre (v2.0.z) DSL offsets removed from field is now included in the file as the Slow mode is dependent on scpot product (Fast/Brst is simply based on offset in Calibration_file).'}];
end
