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
% P 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if regexp(vs,'^P[1-4]$')==1
	v.data = 1;
	v.cl_id = vs(2);
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
	v.labels = {'P'};
	v.field_name = {'Spacecraft potential'};
	v.ent = {'Spacecraft'};
	v.prop = {'Potential'};
	v.fluc = {'Waveform'};
	v.com = ['this signal is averaged from probes ' v.sen];
	v.file = 'mP';
	v.quant = 'p';
elseif regexp(vs,'^P[1-4]_info$')==1
	v.data = 0;
	v.cl_id = vs(2);
	v.inst = 'EFW';
	v.com = 'Spacecraft potential INFO';
	v.file = 'mP';
	v.quant = 'p';
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
	v.field_name = {['Probe #' sen ' potential']};
	v.com = '';
	v.file = 'mP';
	v.quant = 'p';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw E p12 and p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^wE[1-4]p(12|32|34)$')
	v.data = 1;
	v.cl_id = vs(3);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'E';
	v.sen = vs(4:6);
	v.cs = {'na'};
	v.rep = {'scalar'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V/m'};
	v.size = [1];
	v.name = {['P' v.sen]};
	v.labels = v.name;
	v.label_1 = {};
	v.field_name = {'probe potential difference'};
	v.com = '';
	v.file = 'mER';
	v.quant = 'e';

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
	v.sen = vvs(6:8);
	v.cs = {'DSI'};
	v.rep = {'xy'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V/m'};
	v.size = [2];
	v.name = {['Es' v.sen]};
	v.labels = v.name;
	v.label_1 = {'"x", "y"'};
	v.field_name = {'Electric field'};
	v.com = 'Ez=0 by definition (not measured).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full resolution E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?di(b)?E[1-4]p1234$')==1
	v.data = 1;
    switch vvs(1:findstr(vvs,'E')-1) % characters before 'E'
        case 'di'
            v.frame = 'sc';
            v.file = 'mEDSI';
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
	v.sen = '1234';
	v.cs = {'DSI'};
	v.rep = {'xy'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V/m'};
	v.size = [2];
	v.name = {['E' v.sen]};
	v.labels = v.name;
	v.label_1 = {'"x", "y"'};
	v.field_name = {'Electric field'};
	v.com = '';
	
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
	v.cs = {'vector>DSI_xyz','scalar>na'};
 	v.units =  {'mV/m','deg'};
	v.si_conv = {'1.0e-3>V/m','1>degree'};
	v.size = [3 1];
	v.name = {'E', 'Theta'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"',''};
	v.field_name = {'Electric field','Elevation of B above the sc spin plane'};
	v.com = com_Ez;

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
	v.cs = {'vector>gse_xyz','scalar>na'};
 	v.units =  {'mV/m','deg'};
	v.si_conv = {'1.0e-3>V/m','1>degree'};
	v.size = [3 1];
	v.name = {'E', 'Theta'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"',''};
	v.field_name = {'Electric field','Elevation of B above the sc spin plane'};
	v.com = com_Ez;
	
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
		v.cs = {'vector>DSI_xyz','scalar>na'};
	else
		v.cs = {'vector>gse_xyz','scalar>na'};
	end
 	v.units =  {'km/s','deg'};
	v.si_conv = {'1.0e3>m/s','1>degree'};
	v.size = [3 1];
	v.name = {'V', 'Theta'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"',''};
	v.field_name = {'Convection velocity','Elevation of B above the sc spin plane'};
	v.com = com_Ez;
	v.file = 'mEdB';
	if vvs(5)=='s' | vvs(7)=='s', v.quant = 'vedbs';
	else, v.quant = 'vedb';
	end

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
	v.cs = {'scalar>na','scalar>na'};
	v.units =  {'cc','V'};
	v.si_conv = {'1.0e-6>1/m^6',''};
	v.size = [1 1];
	v.name = {['V' v.sen], 'Theta'};
	v.labels = {'V', 'Theta'};
	v.label_1 = {'"x", "y", "z"',''};
	v.field_name = {'Plasma density','Spacecraft potential'};
	v.com = 'density NVps is derived from Vps based on empirical fit. It is NOT a true density';
	v.file = 'mP';
	v.quant = 'p';	

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^A[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'Ephemeris';
	v.frame = 'sc';
	v.sig = 'Phase';
	v.sen = '';
	v.cs = {'scalar>na'};
	v.units =  {'deg'};
	v.si_conv = {'1>degree'};
	v.size = [1];
	v.name = {'A'};
	v.labels = {'Phase'};
	v.label_1 = {''};
	v.field_name = {'Spacecraft phase'};
	v.com = '';
	v.file = 'mA';
	v.quant = 'a';
	
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
	v.cs = {'vector>gse'};
	v.units =  {''};
	v.si_conv = {''};
	v.size = [1];
	v.name = {'SAX'};
	v.labels = {'Spin axis'};
	v.label_1 = {''};
	v.field_name = {'Spacecraft spin axis'};
	v.com = '';
	v.file = 'mEPH';
	v.quant = 'sax';	
	
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
	if strcmp(vvs(1:2),'di'), v.cs = {'vector>DSI_xyz'};
	else, v.cs = {'vector>gse_xyz'};
	end
	v.units =  {'km/s'};
	v.si_conv = {'1e3>m/s'};
	v.size = [3];
	v.name = {'V'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"'};
	v.field_name = {'Spacecraft velocity'};
	v.com = '';
	v.file = 'mR';
	v.quant = 'v';
	
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
	if strcmp(vvs(1:2),'di'), v.cs = {'vector>DSI_xyz'};
	else, v.cs = {'vector>gse_xyz'};
	end
	v.units =  {'km'};
	v.si_conv = {'1e3>m'};
	v.size = [3];
	v.name = {'R'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"'};
	v.field_name = {'Spacecraft position'};
	v.com = '';
	v.file = 'mR';
	v.quant = 'r';
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIS N PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^NC(h|p)[1-4]$')
	v.data = 1;
	v.cl_id = vs(end);
	v.inst = 'CIS';
	v.frame = 'sc';
	v.cs = {'scalar>na'};
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
	v.label_1 = {};
	v.com = '';
	v.file = 'mCIS';
	v.quant = 'ncis';
	
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
		v.cs = {'vector>DSI_xyz','scalar>na'};
	else
		v.cs = {'vector>gse_xyz','scalar>na'};
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
	v.label_1 = {'"x", "y", "z"'};
	v.com = '';
	v.file = 'mCIS';
	v.quant = 'vcis';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDI E PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?(di)?EDI[1-4]$')
	v.data = 1;
	if vvs(1)=='i', 
		vvs = vvs(2:end);
		v.frame = 'inertial'; 
	else, v.frame = 'sc';
	end

	v.cl_id = vs(end);
	v.inst = 'EDI';
	v.sig = 'E';
	v.sen = '';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'vector>DSI_xyz','scalar>na'};
	else
		v.cs = {'vector>gse_xyz','scalar>na'};
	end
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V/m'};
	v.size = [3];
	v.name = {'E'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"'};
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
		v.cs = {'vector>DSI_xyz','scalar>na'};
	else
		v.cs = {'vector>gse_xyz','scalar>na'};
	end
 	v.units =  {'nT'};
	v.si_conv = {'1.0e-12>T'};
	v.size = [3];
	v.name = {'B'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"'};
	v.field_name = {'Magnetic field PP'};
	v.com = '';
	v.file = 'mBPP';
	v.quant = 'b';
	
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
		v.cs = {'vector>DSI_xyz','scalar>na'};
	else
		v.cs = {'vector>gse_xyz','scalar>na'};
	end
 	v.units =  {'nT'};
	v.si_conv = {'1.0e-12>T'};
	v.size = [3];
	v.name = {'B'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"'};
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
	else
		v.field_name = {'Magnetic field'};
		v.com = '';
		v.file = 'mB';
		v.quant = 'bfgm';
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
		disp(['getData q   : ' v.quant ]);
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
