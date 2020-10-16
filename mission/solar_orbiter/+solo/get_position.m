function soloPosition = get_position(tint);
% A function to get the position of Solar Orbiter for a given time interval
% Output is a TSeries object with position in the data variable "pos",
% if any positons found for the time interval "tint".
%
% Example
% tint = irf.tint('2020-07-31T00:00:00Z/2020-08-01T15:00:00Z'); % time interval in TT2000 UTC
% soloPosition = solo.get_position(tint);
% irf_plot(soloPosition, '.');

t_step = 60*60; %3600 sec

% Find the most recent metakernel file and load it
sharedPath = '/share/SPICE/'; % SPICE kernels for different missions are found in this folder on IRFU servers
dirs = dir([sharedPath,'Solar-Orbiter/kernels/mk/*pred-mk_v*.tm']);
if size(dirs, 1) > 1
  % Multiple kernels could be found if executing this script at the same time as syncing new kernel files
  error('Found multiple metakernels, please check your folder.');
end
kernelFile = [dirs.folder, filesep, dirs.name];
cspice_furnsh(kernelFile);

% Compute et (SPICE ephemeries time, make use of input in TT2000)
et = tint.start.tts:t_step:tint.stop.tts;

% The position of Solar Orbiter in units of km
pos = cspice_spkpos('solo', et, 'ECLIPJ2000', 'LT+S', 'Sun');
% More frames ('ECLIPJ2000') described in/share/SPICE/Solar-Orbiter/kernels/fk/solo_ANC_soc-sci-fk_V06.tf

% Convert to utc and then to TT2000 in TSeries object
utc_tmp = cspice_et2utc(et, 'ISOC', 0);

% Note pos' since it is returned as 3xN but TSeries expexcts Nx3 (where N is number of records).
soloPosition = irf.ts_vec_xyz(EpochTT(utc_tmp), pos');
soloPosition.units = 'km'; % Add some metadata information (read when plotting in irf_plot)

end
