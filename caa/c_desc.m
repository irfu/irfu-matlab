function varargout = c_desc(vs,varargin)
%C_DESC provide a description of a Cluster variable
% C_DESC(VS) prints out a descriptoon of variable VS
%
% DESC = C_DESC(VS) returns a descriptoon of variable VS as structure DESC.
%
% Input:
%	VS - variable string
%
% Output:
%	DESC - structure containing a description of variable VS. It has the 
%	following fields:
%	cl_id		%Clister ID
%	inst		%Instrument
%	frame		%Reference frame
%	sig			%Signal
%	sen			%Sensor
%	cs			%Coord
% 	units		%Units
%	si_conv		%SI conversion (CEF)
%	size		%Data dimention (scalar=1)
%	name		%Name of a CEF variable
%	labels		%Label of a CEF variable
%	field_name	%Description of a CEF variable
%	com			%Comment
%	file		%Matlab file name (mXXX.mat)
%	quant		%Quantity name to use with getData
%
% Example:
%	desc = c_desc('diE2');
%
% $Id$

% Copyright 2004 Yuri Khotyaintsev (yuri@irfu.se)
%
error(nargchk(1,10,nargin))

if ~isstr(vs), error('VS must be string'), end

com_Ez = 'Ez is not reliable when magnetic field is close to the spin plane';

vvs = 'XXXXXXXXXX';
vvs(1:length(vs)) = vs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if regexp(vs,'^P[1-4]$')==1
%if strcmp(vs,'P1')|strcmp(vs,'P2')|strcmp(vs,'P3')|strcmp(vs,'P4')
	v.cl_id = vs(2);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = 'all';
	v.cs = {'scalar>na'};
 	v.units =  {'V'};
	v.si_conv = {''};
	v.size = 1;
	v.name = {'Spacecraft_potential'};
	v.labels = {['P' v.sen]};
	v.field_name = {'Averaged probe potential from all probes'};
	v.com = 'this signal is averaged from all probes available at the time';
	v.file = 'mP';
	v.quant = 'p';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P - individual probes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^P10Hz[1-4]p[1-4]$')==1
%elseif strcmp(vvs(1:4),'P10Hz')==1 & is14(vvs(5)) & vvs(6)=='p' & is14(vvs(7))
	v.cl_id = vs(6);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = vs(8);
	v.cs = {'scalar>na'};
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
elseif regexp(vs,'^wE[1-4]p(12|34)$')==1
%elseif strcmp(vvs(1:2),'wE') & is14(vvs(3)) & vvs(4)=='p' & ...
%(strcmp(vvs(5:6),'12') | strcmp(vvs(5:6),'34')) 
	v.cl_id = vs(3);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'E';
	v.sen = vs(4:6);
	v.cs = {'scalar>na'};
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
elseif regexp(vs,'^diEs[1-4]p(12|34)')==1
%elseif strcmp(vvs(1:4),'diEs') & is14(vvs(5)) & vvs(6)=='p' & ...
%(strcmp(vvs(7:8),'12') | strcmp(vvs(7:8),'34'))
	v.cl_id = vs(5);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'E';
	v.sen = vs(6:8);
	v.cs = {'vector>dsi_xy'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V/m'};
	v.size = [2];
	v.name = {['Es' v.sen]};
	v.labels = v.name;
	v.label_1 = {'"x", "y"'};
	v.field_name = {'Electric field'};
	v.com = '';
	v.file = 'mEDSI';
	v.quant = 'dies';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full resolution E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^diE[1-4]p1234$')==1
%elseif (length(vs)==9 & strcmp(vs(1:3),'diE') & strcmp(vvs(5:9),'p1234'))
	v.cl_id = vs(4);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'E';
	v.sen = 'all';
	v.cs = {'vector>dsi_xy'};
 	v.units =  {'mV/m'};
	v.si_conv = {'1.0e-3>V/m'};
	v.size = [2];
	v.name = {['E' v.sen]};
	v.labels = v.name;
	v.label_1 = {'"x", "y"'};
	v.field_name = {'Electric field'};
	v.com = '';
	v.file = 'mEDSI';
	v.quant = 'die';
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full resolution E with assumption E.B = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(diE[1-4]|diEs[1-4])$')==1
%elseif (length(vs)==4 & strcmp(vvs(1:3),'diE') & is14(vvs(4))) | ...
%(length(vs)==5 & strcmp(vvs(1:4),'diEs') & is14(vvs(5)))
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'E';
	if vs(4)=='s', v.sen = 's'; else, v.sen = ''; end
	v.cs = {'vector>dsi_xyz','scalar>na'};
 	v.units =  {'mV/m','deg'};
	v.si_conv = {'1.0e-3>V/m','1>degree'};
	v.size = [3 1];
	v.name = {'E', 'Theta'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"',''};
	v.field_name = {'Electric field','Elevation of B above the sc spin plane'};
	v.com = com_Ez;
	v.file = 'mEdB';
	if vs(4)=='s', v.quant = 'edbs'; else, v.quant = 'edb'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full resolution E in GSE coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(E[1-4]|Es[1-4])$')==1
%elseif (length(vs)==2 & vvs(1)=='E' & is14(vvs(2))) | ...
%(length(vs)==3 & strcmp(vvs(1:2),'Es') & is14(vvs(3)))
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
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
	v.file = 'mEdB';
	if vs(2)=='s', v.quant = 'edbs'; else, v.quant = 'edb'; end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(diVExBs[1-4]|VExBs[1-4])$')==1
%elseif (((length(vs)==8 & vvs(7)=='s')|length(vs)==7) & ...
%strcmp(vvs(1:6),'diVExB') & is14(vs(end))) | ...
%(((length(vs)==6 & vvs(5)=='s')|length(vs)==5) & strcmp(vvs(1:4),'VExBs') & is14(vs(end)))
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'V=ExB';
	if vvs(5)=='s' | vvs(7)=='s', v.sen = 's';
	else, v.sen = ''; 
	end
	if strcmp(vvs(1:2),'di')
		v.cs = {'vector>dsi_xyz','scalar>na'};
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
%elseif (length(vs)==5 & strcmp(vvs(1:4),'NVps') & is14(vvs(5)))
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
%elseif length(vs)==2 & vs(1)=='A' & is14(vvs(2))
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
elseif regexp(vs,'^V[1-4]$')
	v.cl_id = vs(end);
	v.inst = 'Ephemeris';
	v.frame = 'sc';
	v.sig = 'Attitude';
	v.sen = '';
	v.cs = {'vector>gse'};
	v.units =  {'km/s'};
	v.si_conv = {'1e3>m/s'};
	v.size = [3];
	v.name = {'V'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"'};
	v.field_name = {'Spacecraft velocity'};
	v.com = '';
	v.file = 'mR';
	v.quant = 'Vsc';	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIS N PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^NC(h|p)[1-4]$')

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
elseif regexp(vs,'^diVC(h|p)[1-4]$') | regexp(vs,'^VC(h|p)[1-4]$')
	v.cl_id = vs(end);
	v.inst = 'CIS';
	v.frame = 'sc';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'vector>dsi_xyz','scalar>na'};
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
elseif regexp(vs,'^diEDI[1-4]|EDI[1-4]$')
	v.cl_id = vs(end);
	v.inst = 'EDI';
	v.frame = 'inertial';
	v.sig = 'E';
	v.sen = '';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'vector>dsi_xyz','scalar>na'};
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
elseif regexp(vs,'^diBPP[1-4]|BPP[1-4]$')
	v.cl_id = vs(end);
	v.inst = 'FGM';
	v.frame = 'sc';
	v.sig = 'B';
	v.sen = '';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'vector>dsi_xyz','scalar>na'};
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
elseif regexp(vs,'^diB[1-4]|B[1-4]')
	v.cl_id = vs(end);
	v.inst = 'FGM';
	v.frame = 'sc';
	v.sig = 'B';
	v.sen = '';
	if strcmp(vvs(1:2),'di')
		vvs = vvs(3:end);
		v.cs = {'vector>dsi_xyz','scalar>na'};
	else
		v.cs = {'vector>gse_xyz','scalar>na'};
	end
 	v.units =  {'nT'};
	v.si_conv = {'1.0e-12>T'};
	v.size = [3];
	v.name = {'B'};
	v.labels = v.name;
	v.label_1 = {'"x", "y", "z"'};
	v.field_name = {'Magnetic field'};
	v.com = '';
	v.file = 'mB';
	v.quant = 'bfgm';
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
		disp(['Ref Frame   : ' v.frame ]);
		disp(['Signal      : ' v.sig ]);
		disp(['Sensor      : ' v.sen ]);
		disp(['Comment     : ' v.com ]);
		disp(['Size        : ' num2str(v.size) ]);
		for j=1:length(v.size)
			disp(['Var Name    : ' v.name{j} ]);
			disp(['  Labels    : ' v.labels{j} ]);
			disp(['  Field     : ' v.field_name{j} ]);
			disp(['  Coord Sys : ' v.cs{j} ]);
			disp(['  Units     : ' v.units{j} ]);
			disp(['  SI conv   : ' v.si_conv{j} ]);
		end
	end
end

function r = is14(s)

r = 0;
if s=='1' | s=='2' | s=='3' | s=='4', r = 1; end
