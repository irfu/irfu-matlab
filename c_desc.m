function varargout = c_desc(vs,v_info)
%C_DESC provide a description of a Cluster variable
% C_DESC(V_S) prints out a descriptoon of variable VS
%
% DESC = C_DESC(V_S) returns a descriptoon of variable VS as structure DESC.
%
% Input:
%	V_S - string defining a variable
%
% Output:
%	DESC - structure containing a description of variable VS. It has the 
%	following fields:
%   cl_id		%Clister ID
%   inst		%Instrument
%   frame		%Reference frame
%   sig			%Signal
%   sen			%Sensor
%   cs			%Coord
%   units		%Units
%   si_conv		%SI conversion (CEF)
%   size		%Data dimention (scalar=1)
%   name		%Name of a CEF variable
%   labels		%Label of a CEF variable
%   field_name	%Description of a CEF variable
%   com			%Comment
%   file		%Matlab file name (mXXX.mat)
%   quant		%Quantity name to use with getData
%   lev         %Data level: 0 - ClusterDB, 1 - ClusterProc, 2 - manual
%
% Examples:
% c_desc('diE2')
%	desc = c_desc('diE2');
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)

error(nargchk(1,2,nargin))
if ~isstr(vs), error('VS must be string'), end
if nargin<2, v_info = []; end

com_Ez = 'Ez is not reliable when magnetic field is close to the spin plane';

vvs = 'XXXXXXXXXX';
vvs(1:length(vs)) = vs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P & Ps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if regexp(vs,'^P(s)?[1-4]$')==1
	v.data = 1;
	if vs(2)=='s', v.cl_id = vs(3);
	else, v.cl_id = vs(2);
	end
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	if ~isempty(v_info) & isfield(v_info,'probe'), v.sen = num2str(v_info.probe);
	else, v.sen = '1234';
	end 
	v.cs = {'na'};
	v.rep = {'scalar'};
 	v.units =  {'V'};
	v.si_conv = {''};
	v.size = 1;
	v.name = {'Spacecraft_potential'};
	v.labels = {'-Sc pot'};
	if vs(2)=='s'
		v.quant = 'ps';
		v.field_name = {'Spacecraft potential (spin resolution)'};
	else
		v.quant = 'p';
		v.field_name = {'Spacecraft potential'};
	end
	v.ent = {'Instrument'};
	v.prop = {'Probe_Potential'};
	v.fluc = {'Waveform'};
	v.com = ['this signal is averaged from probes ' v.sen];
	v.file = 'mP';
	v.lev = 1;
elseif regexp(vs,'^P(s)?[1-4]_info$')==1
	v.data = 0;
	if vs(2)=='s', v.cl_id = vs(3);
	else, v.cl_id = vs(2);
	end
	v.inst = 'EFW';
	v.com = 'Spacecraft potential INFO';
	v.file = 'mP';
	if vs(2)=='s', v.quant = 'ps';
	else, v.quant = 'p';
	end
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P - individual probes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^P10Hz[1-4]p[1-4]$')==1
	v.data = 1;
	v.cl_id = vs(6);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = vs(8);
	v.cs = {'na'};
	v.rep = {'scalar'};
 	v.units =  {'V'};
	v.si_conv = {''};
	v.size = 1;
	v.name = {['P' v.sen]};
	v.labels = v.name;
	v.field_name = {['Probe ' v.sen ' to spacecraft potential']};
	v.ent = {'Instrument'};
	v.prop = {'Probe_Potential'};
	v.fluc = {'Waveform'};
	v.com = '';
	v.file = 'mPR';
	v.quant = 'p';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P - individual probes from internal burst
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^P(32|4)kHz[1-4]p[1-4]$')==1
	v.data = 1;
	v.cl_id = vs(end-2);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = vs(end);
	v.cs = {'na'};
	v.rep = {'scalar'};
 	v.units =  {'V'};
	v.si_conv = {''};
	v.size = 1;
	v.name = {['P' v.sen]};
	v.labels = v.name;
	v.field_name = {['Probe ' v.sen ' to spacecraft potential']};
	v.com = '';
	v.file = 'mEFWburst';
	v.quant = 'pburst';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw E p12 and p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^wE[1-4]p(12|32|34)$')
	v.data = 1;
	v.cl_id = vs(3);
	v.inst = 'EFW';
	v.frame = 'na';
	v.sig = 'E';
	v.sen = vs(5:6);
	v.cs = {'SC'};
	v.rep = {'scalar'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V m^-1'};
	v.size = [1];
	v.name = {['P' v.sen]};
	v.labels = v.name;
	v.field_name = {['Electric field component measured between the probes '...
		v.sen(1) ' and ' v.sen(2)]};
	v.ent = {'Electric_Field'};
	v.prop = {'Component'};
	v.fluc = {'Waveform'};
	v.com = '';
	v.file = 'mER';
	v.quant = 'e';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spin fits E p12 and p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?diEs[1-4]p(12|32|34)')==1
	v.data = 1;
	if vvs(1)=='i', 
		vvs = vvs(2:end);
		v.frame = 'inertial'; 
		v.file = 'mEIDSI';
		v.quant = 'idies';
	else
		v.frame = 'sc';
		v.file = 'mEDSI';
		v.quant = 'dies';
	end
	v.cl_id = vvs(5);
	v.inst = 'EFW';
	v.sig = 'E';
	v.sen = vvs(7:8);
	v.cs = {'ISR2','na'};
	v.rep = {'xy',''};
 	v.units =  {'mV/m','unitless'};
	v.si_conv = {'1.0e-3>V m^-1',''};
	v.size = [2 1];
	v.name = {'E_Vec_xy_ISR2', 'E_sigma'};
	v.labels = {'E','St dev'};
	v.label_1 = {'"Ex", "Ey"',''};
    v.col_labels = {{'x','y','z'},''};
    v.rep_1 = {'"x", "y"',''};
	v.field_name = {'Electric field (spin resolution)','Standard deviation'};
	v.ent = {'Electric_Field','Electric_Field'};
	v.prop = {'Vector','Vector'};
	v.fluc = {'Waveform','Fluctuation_Level'};
	v.com = 'Ez=0 by definition (not measured).';
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full resolution E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?di(b)?E[1-4]p1234$')==1
	v.data = 1;
    switch vvs(1:findstr(vvs,'E')-1) % characters before 'E'
        case 'di'
            v.frame = 'sc';
            v.file = 'mEDSIf';
            v.quant = 'die';
        case 'idi',
            v.frame = 'inertial';
            v.file = 'mEIDSI';
            v.quant = 'idie';
        case 'dib' % internal burst mode E field
            v.frame = 'sc';
            v.file = 'mEFWburst';
            v.quant = 'dibe';
        case 'idib' % internal burst mode E field
            v.frame = 'inertial';
            v.file = 'mEFWburst';
            v.quant = 'idibe';
    end
	v.cl_id = vvs(findstr(vvs,'E')+1); % next character after 'E'
	v.inst = 'EFW';
	v.sig = 'E';
	% TODO: we need to check here for p32-p34 - v.sen = '3234';
	v.sen = '1234';
	v.cs = {'ISR2'};
	v.rep = {'xy'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V m^-1'};
	v.size = [2];
	v.name = {'E_Vec_xy_ISR2'};
	v.labels = {'E'};
	v.label_1 = {'"Ex", "Ey"'};
	v.rep_1 = {'"x", "y"'};
	v.col_labels = {{'x','y','z'}};
	v.field_name = {'Electric field'};
	v.ent = {'Electric_Field'};
	v.prop = {'Vector'};
	v.fluc = {'Waveform'};
	v.com = '';
	v.lev = 1;
elseif regexp(vs,'^diE[1-4]p1234_info$')==1
	v.data = 0;
	v.cl_id = vs(4);
	v.inst = 'EFW';
	v.com = 'E full res INFO';
	v.file = 'mEDSIf';
	v.quant = 'die';
	v.lev = 1;	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^IBIAS[1-4]p[1-4]$')==1
	v.data = 1;
	v.cl_id = vs(6);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'BIAS';
	v.sen = vs(end);
	v.cs = {'na'};
	v.rep = {'scalar'};
 	v.units =  {'nA'};
	v.si_conv = {'1.0e-9>A'};
	v.size = 1;
	v.name = {['I-bias-p' v.sen]};
	v.labels = v.name;
	v.field_name = {['Probe ' v.sen ' bias current']};
	v.com = '';
	v.file = 'mFDM';
	v.quant = 'ibias';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^DSC[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(4);
	v.inst = 'EFW';
	v.com = 'DSC';
	v.file = 'mFDM';
	v.quant = 'dsc';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FDM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^FDM[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(4);
	v.inst = 'EFW';
	v.com = 'FDM';
	v.file = 'mFDM';
	v.quant = 'fdm';
	v.lev = 0;	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Whisper pulses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^WHIP[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(5);
	v.inst = 'EFW';
	v.com = 'Whisper pulses';
	v.file = 'mFDM';
	v.quant = 'whip';
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sweep + dump
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^SWEEP[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(6);
	v.inst = 'EFW';
	v.com = 'Sweep';
	v.file = 'mFDM';
	v.quant = 'sweep';
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Burst dump
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^BDUMP[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(6);
	v.inst = 'EFW';
	v.com = 'Burst dump';
	v.file = 'mFDM';
	v.quant = 'bdump';
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad bias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^BADBIAS[1-4]p[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(8);
	v.sen = vs(end);
	v.inst = 'EFW';
	v.com = 'Bad bias current';
	v.file = 'mFDM';
	v.quant = 'badbias';
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC offsets corse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^Da[1-4]p(12|32|34)$')==1
	v.data = 0;
	v.cl_id = vs(3);
	v.sen = vvs(5:6);
	v.inst = 'EFW';
	v.com = 'ADC offset';
	v.file = 'mEDSIf';
	v.quant = 'die';
	v.lev = 1;	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADC offsets from spinfits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^Dadc[1-4]p(12|32|34)$')==1
	v.data = 0;
	v.cl_id = vs(5);
	v.sen = vvs(7:8);
	v.inst = 'EFW';
	v.com = 'ADC offset';
	v.file = 'mEDSI';
	v.quant = 'dies';
	v.cs = {'na'};
	v.rep = {'scalar'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V m^-1'};
	v.size = [1];
	v.name = {['Dadc-p' v.sen]};
	v.labels = v.name;
	v.field_name = {'ADC offset'};
	v.ent = {'Instrument'};
	v.prop = {'ADC_Offset'};
	v.fluc = {'Waveform'};
	v.com = '';
	v.lev = 1;	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Delta offsets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^D[1-4]p12p34$')==1
	v.data = 0;
	v.cl_id = vs(2);
	v.inst = 'EFW';
	v.com = 'Delta offset';
	v.file = 'mEDSI';
	v.quant = 'dies';
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DSI offsets
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^Ddsi[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(5);
	v.inst = 'EFW';
	v.com = 'DSI offsets';
	v.file = 'mEDSI';
	v.quant = '';
	v.lev = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amplitude correction factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^Damp[1-4]$')==1
	v.data = 0;
	v.cl_id = vs(5);
	v.inst = 'EFW';
	v.com = 'Amplitude correction factor for E';
	v.file = 'mEDSI';
	v.quant = '';
	v.lev = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full/spin resolution E with assumption E.B = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?diE(s)?[1-4]$')
	v.data = 1;
	if vvs(1)=='i', 
		vvs = vvs(2:end);
		v.frame = 'inertial'; 
		v.file = 'mEdBI';
		if vvs(4)=='s', v.quant = 'iedbs'; else, v.quant = 'iedb'; end
	else
		v.frame = 'sc';
		v.file = 'mEdB';
		if vvs(4)=='s', v.quant = 'edbs'; else, v.quant = 'edb'; end
	end
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.sig = 'E';
	if vvs(4)=='s', v.sen = 's'; else, v.sen = ''; end
	v.cs = {'ISR2', 'na', 'na'};
 	v.units =  {'mV/m','deg','unitless'};
	v.si_conv = {'1.0e-3>V m^-1','1>degree','1>unitless'};
	v.size = [3 1 1];
	v.name = {'E','Theta','E_sigma'};
	v.labels = v.name;
	v.label_1 = {'"Ex", "Ey", "Ez"','',''};
	v.col_labels = {{'x','y','z'},'',''};
	v.field_name = {'Electric field','Elevation of B above the sc spin plane','Standard deviation'};
	v.com = com_Ez;
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full/spin resolution E in GSE coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?E(s)?[1-4]$')
	v.data = 1;
	if vvs(1)=='i', 
		vvs = vvs(2:end);
		v.frame = 'inertial'; 
		v.file = 'mEdBI';
		if vvs(4)=='s', v.quant = 'iedbs'; else, v.quant = 'iedb'; end
	else
		v.frame = 'sc';
		v.file = 'mEdB';
		if vvs(4)=='s', v.quant = 'edbs'; else, v.quant = 'edb'; end
	end
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.sig = 'E';
	if vs(2)=='s', v.sen = 's'; else, v.sen = ''; end
	v.cs = {'GSE', 'na','na'};
 	v.units =  {'mV/m','deg','unitless'};
	v.si_conv = {'1.0e-3>V m^-1','1>degree','1>unitless'};
	v.size = [3 1];
	v.name = {'E','Theta','E_sigma'};
	v.labels = v.name;
	v.label_1 = {'"Ex", "Ey", "Ez"','',''};
	v.col_labels = {{'x','y','z'},'',''};
	v.field_name = {'Electric field','Elevation of B above the sc spin plane','Standard deviation'};
	v.com = com_Ez;
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?VExB(s)?[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'V=ExB';
	if vvs(5)=='s' | vvs(7)=='s', v.sen = 's';
	else, v.sen = ''; 
	end
	if strcmp(vvs(1:2),'di')
		v.cs = {'ISR2', 'na'};
	else
		v.cs = {'GSE', 'na'};
	end
 	v.units =  {'km/s','deg'};
	v.si_conv = {'1.0e3>m/s','1>degree'};
	v.size = [3 1];
	v.name = {'V', 'Theta'};
	v.labels = v.name;
	v.label_1 = {'"Vx", "Vy", "Vz"',''};
	v.col_labels = {{'x','y','z'},''};
	v.field_name = {'Convection velocity','Elevation of B above the sc spin plane'};
	v.com = com_Ez;
	v.file = 'mEdB';
	if vvs(5)=='s' | vvs(7)=='s', v.quant = 'vedbs';
	else, v.quant = 'vedb';
	end
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full resolution satellite potential and derived density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^NVps[1-4]$')==1
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = 'all';
	v.cs = {'na','na'};
	v.units =  {'cc','V'};
	v.si_conv = {'1.0e-6>1/m^6',''};
	v.size = [1 1];
	v.name = {'Plasma density', 'Spacecraft potential'};
	v.labels = {'Nscp', '-Sc pot'};
	v.field_name = {'Plasma density','Spacecraft potential'};
	v.com = 'density NVps is derived from Vps based on empirical fit. It is NOT a true density';
	v.file = 'mP';
	v.quant = 'p';	
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase and phase_2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^A(two)?[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'Ephemeris';
	v.frame = 'sc';
	v.sig = 'Phase';
	v.sen = '';
	v.cs = {'na'};
	v.units =  {'deg'};
	v.si_conv = {'1>degree'};
	v.size = [1];
	if strcmp(vs(2),'t')
		v.name = {'Atwo'};
		v.labels = {'Phase_2'};
	else
		v.name = {'A'};
		v.labels = {'Phase'};
	end
	v.field_name = {'Spacecraft phase'};
	v.com = '';
	v.file = 'mA';
	v.quant = 'a';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spin axis orientation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^SAX[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'Ephemeris';
	v.frame = 'sc';
	v.sig = 'Attitude';
	v.sen = '';
	v.cs = {'GSE'};
	v.units =  {''};
	v.si_conv = {''};
	v.size = [1];
	v.name = {'SAX'};
	v.labels = {'Spin axis'};
	v.field_name = {'Spacecraft spin axis'};
	v.com = '';
	v.file = 'mEPH';
	v.quant = 'sax';	
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spacecraft velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?V[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'Ephemeris';
	v.frame = 'sc';
	v.sig = 'Velocity';
	v.sen = '';
	if strcmp(vvs(1:2),'di'), v.cs = {'ISR2'};
	else, v.cs = {'GSE'};
	end
	v.units =  {'km/s'};
	v.si_conv = {'1e3>m/s'};
	v.size = [3];
	v.name = {'V'};
	v.labels = v.name;
	v.label_1 = {'"Vx", "Vy", "Vz"'};
	v.col_labels = {{'x','y','z'},''};
	v.field_name = {'Spacecraft velocity'};
	v.com = '';
	v.file = 'mR';
	v.quant = 'v';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spacecraft position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?R[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'Ephemeris';
	v.frame = 'sc';
	v.sig = 'Position';
	v.sen = '';
	if strcmp(vvs(1:2),'di'), v.cs = {'ISR2'};
	else, v.cs = {'GSE'};
	end
	v.units =  {'km'};
	v.si_conv = {'1e3>m'};
	v.size = [3];
	v.name = {'R'};
	v.labels = v.name;
	v.label_1 = {'"Rx", "Ry", "Rz"'};
	v.col_labels = {{'x','y','z'},''};
	v.field_name = {'Spacecraft position'};
	v.com = '';
	v.file = 'mR';
	v.quant = 'r';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIS N PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^NC(h|p)[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'CIS';
	v.frame = 'sc';
	v.cs = {'na'};
	if vvs(3)=='h'
		v.sig = 'N';
		v.sen = 'HIA';
		v.field_name = {'Ion density'};
	else
		v.sig = 'Np';
		v.sen = 'COD';
		v.field_name = {'Proton density'};
	end
 	v.units =  {'cc'};
	v.si_conv = {'1.0e6>1/m^6'};
	v.size = [1];
	v.name = {'N'};
	v.labels = v.name;
	v.com = '';
	v.file = 'mCIS';
	v.quant = 'ncis';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIS V PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?VC(h|p)[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'CIS';
	v.frame = 'sc';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'ISR2'};
	else
		v.cs = {'GSE'};
	end
	if vvs(3)=='h'
		v.sig = 'V';
		v.sen = 'HIA';
		v.field_name = {'Ion flow velocity'};
	else
		v.sig = 'Vp';
		v.sen = 'COD';
		v.field_name = {'Proton flow velocity'};
	end
 	v.units =  {'km/s'};
	v.si_conv = {'1.0e3>m/s'};
	v.size = [3];
	v.name = {'V'};
	v.labels = v.name;
	v.label_1 = {'"Vx", "Vy", "Vz"'};
	v.col_labels = {{'x','y','z'},''};
	v.com = '';
	v.file = 'mCIS';
	v.quant = 'vcis';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIS VxB PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?VCE(h|p)[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'CIS';
	v.frame = 'sc';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'ISR2'};
	else
		v.cs = {'GSE'};
	end
	if vvs(3)=='h'
		v.sig = 'V';
		v.sen = 'HIA';
		v.field_name = {'Ion VxB'};
	else
		v.sig = 'Vp';
		v.sen = 'COD';
		v.field_name = {'Proton VxB'};
	end
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V m^-1'};
	v.size = [3];
	v.name = {'E'};
	v.labels = v.name;
	v.label_1 = {'"Ex", "Ey", "Ez"'};
	v.col_labels = {{'x','y','z'},''};
	v.com = '';
	v.file = 'mCIS';
	v.quant = 'vce';
	v.lev = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDI E PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?(di)?EDI[1-4]$')
	v.data = 1;
	if vvs(1)=='i', 
		vvs = vvs(2:end);
		v.frame = 'inertial'; 
		v.lev = 0;
	else
		v.frame = 'sc';
		v.lev = 1;
	end

	v.cl_id = vs(end);
	v.inst = 'EDI';
	v.sig = 'E';
	v.sen = '';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'ISR2'};
	else
		v.cs = {'GSE'};
	end
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V m^-1'};
	v.size = [3];
	v.name = {'E'};
	v.labels = v.name;
	v.label_1 = {'"Ex", "Ey", "Ez"'};
	v.col_labels = {{'x','y','z'},''};
	v.field_name = {'Perpendicular Electric field'};
	v.com = '';
	v.file = 'mEDI';
	v.quant = 'edi';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FGM B PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?BPP[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'FGM';
	v.frame = 'sc';
	v.sig = 'B';
	v.sen = '';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'ISR2'};
	else
		v.cs = {'GSE'};
	end
 	v.units =  {'nT'};
	v.si_conv = {'1.0e-9>T'};
	v.size = [3];
	v.name = {'B'};
	v.labels = v.name;
	v.label_1 = {'"Bx", "By", "Bz"'};
	v.col_labels = {{'x','y','z'},''};
	v.field_name = {'Magnetic field PP'};
	v.com = '';
	v.file = 'mBPP';
	v.quant = 'b';
	v.lev = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FGM B full resolution and resampled
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?B(r|rs)?[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'FGM';
	v.frame = 'sc';
	v.sig = 'B';
	v.sen = '';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'ISR2'};
	else
		v.cs = {'GSE'};
	end
 	v.units =  {'nT'};
	v.si_conv = {'1.0e-9>T'};
	v.size = [3];
	v.name = {'B'};
	v.labels = v.name;
	v.label_1 = {'"Bx", "By", "Bz"'};
	v.col_labels = {{'x','y','z'}};
	if regexp(vs,'^(di)?B(r|rs)[1-4]$')
		v.file = 'mBr';
		if regexp(vs,'^(di)?Brs[1-4]$')
			v.field_name = {'Magnetic field resampled to E (spin resolution)'};
			v.com = 'Resampled to E (spin resolution)';
			v.quant = 'brs';
		else
			v.field_name = {'Magnetic field resampled to E (full resolution)'};
			v.com = 'Resampled to E (full resolution)';
			v.quant = 'br';
		end
		v.lev = 1;
	else
		v.field_name = {'Magnetic field'};
		v.com = '';
		v.file = 'mB';
		v.quant = 'bfgm';
		v.lev = 0;
	end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDITIONAL HELP IN PLOTTING, NOT PARTICULAR TOCLUSTER 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(vs,'B') | strcmp(vs,'j') | strcmp(vs,'jz') | strcmp(vs,'jxB')
	v.data = 1;
	v.cl_id = '';
	v.inst = NaN;
	v.frame = NaN;
    if strcmp(vs,'B')
        v.units =  {'nT'};
        v.name = {'B'};
        v.labels = v.name;
        v.label_1 = {'"Bx", "By", "Bz"'};
    elseif strcmp(vs,'j')
        v.units =  {'A/m^2'};
        v.name = {'J'};
        v.labels = v.name;
        v.label_1 = {'"Jx", "Jy", "Jz"'};
    elseif strcmp(vs,'jz')
        v.units =  {'A/m^2'};
        v.name = {'J_{||}'};
        v.labels = v.name;
        v.label_1 = {'"Jx", "Jy", "Jz"'};
    elseif strcmp(vs,'jxB')
        v.units =  {'T A'}; 
        v.name = {'JxB'};
        v.labels = v.name;
        v.label_1 = {'"JxBx", "JxBy", "JxBz"'};
    end
else
	error('Variable name not recognized')
end

if nargin>2, have_options = 1; args = varargin;
% construct the output
else 
	if nargout==1, varargout = {v};
	else
		%just print out everthing
		disp(['SC#         : ' v.cl_id ]);
		disp(['Instrument  : ' v.inst ]);
		if v.data
			disp(['Ref Frame   : ' v.frame ]);
			disp(['Signal      : ' v.sig ]);
			disp(['Sensor      : ' v.sen ]);
		end
		disp(['Comment     : ' v.com ]);
		disp(['File        : ' v.file ]);
		if v.lev<2
			disp(['getData q   : ' v.quant ]);
			if v.lev, disp('getData cl  : ClusterProc');
			else, disp('getData cl  : ClusterDB');
			end
		else
			disp('getData     : manual processing');
		end
		
		if v.data
			disp(['Size        : ' num2str(v.size) ]);
			for j=1:length(v.size)
				disp(['Var Name    : ' v.name{j} ]);
				disp(['  Labels    : ' v.labels{j} ]);
				disp(['  Field     : ' v.field_name{j} ]);
				disp(['  Coord Sys : ' v.cs{j} ]);
				if isfield(v,'rep'), disp(['  Represent : ' v.rep{j} ]); end
				disp(['  Units     : ' v.units{j} ]);
				disp(['  SI conv   : ' v.si_conv{j} ]);
			end
		end
	end
end

function r = is14(s)

r = 0;
if s=='1' | s=='2' | s=='3' | s=='4', r = 1; end
