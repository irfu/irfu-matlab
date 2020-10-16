%A function to get the position of Solar Orbiter for a given time interval
% Usage: solo = get_solo_pos(t_start, t_end);  plot(solo.t, solo.pos,'.')
% t_start and t_end in some time format accepted by datestr.m

function solopos = get_position(t_start, t_end);

t_step=60*60; %3600 sec

%Find the most recent metakernel file and load it
dirs=dir('/share/SPICE/Solar-Orbiter/kernels/mk/*pred-mk_v*.tm');
kernelFile = ['/share/SPICE/Solar-Orbiter/kernels/mk/' dirs.name];
cspice_furnsh(kernelFile);

%Set start and stop time of wanted interval and convert to et (SPICE ephemeries time)
int1_cspice= datestr(t_start,'YYYY mm dd HH MM SS'); % e.g. '2020 02 11 00 00 00';
int2_cspice= datestr(t_end,'YYYY mm dd HH MM SS'); % e.g .'2020 12 05 00 00 00';
et=cspice_str2et(int1_cspice):t_step:cspice_str2et(int2_cspice);

%The position of Solar Orbiter in units of km
solopos.pos=cspice_spkpos('solo',et,'ECLIPJ2000', 'LT+S','Sun'); %More frames ('ECLIPJ2000') described in/share/SPICE/Solar-Orbiter/kernels/fk/solo_ANC_soc-sci-fk_V06.tf

%Convert to utc and then matlab-time
utc_tmp=cspice_et2utc(et,'ISOC',0);
solopos.t=irf_time(utc_tmp,'utc>date');


