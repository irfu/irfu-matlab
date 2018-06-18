function out = thor_orbit(orbitFile,dt)
% returns THOR orbit in orbitfile with time step dt in GSE.
% If dt is not given, the default resolution of output is 10min.
% The output is a TSeries vector.
%
% R = THOR_ORBIT(orbitFile,dt)
%
% Orbit is read from SPICE directory specified by 
% datastore('spice','dir',fullSPICEDirectoryPath)
% The code assumes that under SPICE is directory is subdirectory THOR with
% containing the orbit files.
%
% R = THOR_ORBIT
%
% Returns the default pre-MCR orbit (new1a.bsp) at 10min resolution.

%% Check input
if nargin == 0 && nargout == 0
	help thor_read_orbit_spis;
	return;
elseif nargin == 0 && nargout > 0
	orbitFile = 'new1a.bsp';
end
if nargin <= 1
	dt = 10*60;
end

%% Define defaults
out = datastore('spice');
if isfield(out,'dir')
	spiceDirectory = datastore('spice','dir');
else
	disp('SPICE directory underfined!');
	disp('Please define, for example:');
	disp('  > datastore(''spice'',''dir'',''/Users/andris/calc/SPICE'');');
	return;
end
orbitFileFullPath  = [spiceDirectory '/THOR/' orbitFile];
idTHOR = cspice_spkobj(orbitFileFullPath,100); % get the THOR id number used by ESOC (-666 in alt1a.bsp and -999 in new1a.bsp

%% Read in orbits
% Provide a name for the satellite ID in the bsp-files
% Need to do this since spice hasn't got official support for Solar Orbiter
% yet.
cspice_boddef('THOR',idTHOR)
cspice_furnsh([spiceDirectory '/generic_kernels/spk/planets/de421.bsp']); % planet orbits
cspice_furnsh([spiceDirectory '/generic_kernels/lsk/naif0011.tls']);      % leap seconds
cspice_furnsh([spiceDirectory '/metakernel_frames.txt']); % defines GSE,..
cspice_furnsh(orbitFileFullPath);
	
%% Calculate interval
%ids = cspice_spkobj(orbitFileFullPath,1000);
cover=cspice_spkcov(orbitFileFullPath,idTHOR,1000);
tstart=cover(1);
tend=cover(2);
t=tstart:dt:tend;

% some convenient timestamps
%t0=cspice_str2et('2017-01-03T00:00:00.000');
%t1=cspice_str2et('2027-01-01T00:00:00.000');


%% Calculate orbit
%	rTHOR = cspice_spkpos('THOR',t, 'j2000', 'none', 'Earth');
rTHOR = cspice_spkpos('THOR',t,'GSE', 'none', 'Earth');

RTHOR =	irf.ts_vec_xyz(t,rTHOR');
RTHOR.units = 'km';
RTHOR.name = 'THOR orbit';
RTHOR.coordinateSystem = 'GSE';

%% output

out = RTHOR;

