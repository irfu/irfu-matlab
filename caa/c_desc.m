function varargout = c_desc(vs,varargin)
%C_DESC provide a description of a Cluster variable
%
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

error(nargchk(1,10,nargin))

if ~isstr(vs), error('VS must be string'), end

com_Ez = 'Ez is not reliable when magnetic field is close to the spin plane';

vvs = 'XXXXXXXXXX';
vvs(1:length(vs)) = vs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if regexp(vs,'^P[1-4]$')==1
if strcmp(vs,'P1')|strcmp(vs,'P2')|strcmp(vs,'P3')|strcmp(vs,'P4')
	v.cl_id = vs(2);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = '1234';
	v.cs = {'scalar>na'};
 	v.units =  {'V'};
	v.si_conv = {''};
	v.size = 1;
	v.name = {['P_' v.sen]};
	v.labels = {['P' v.sen]};
	v.field_name = {'Averaged probe potential from all probes'};
	v.com = 'this signal is averaged from all probes available at the time';
	v.file = 'mP';
	v.quant = 'p';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P - individual probes 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elseif regexp(vs,'^P10Hz[1-4]p[1-4]$')==1
elseif strcmp(vvs(1:4),'P10Hz')==1 & is14(vvs(5)) & vvs(6)=='p' & is14(vvs(7))
	v.cl_id = vs(6);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = vs(8);
	v.cs = {'scalar>na'};
 	v.units =  {'V'};
	v.si_conv = {''};
	v.size = 1;
	v.name = {['P_' v.sen]};
	v.labels = {['P' v.sen]};
	v.field_name = {['Probe #' sen ' potential']};
	v.com = '';
	v.file = 'mP';
	v.quant = 'p';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% raw E p12 and p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elseif regexp(vs,'^wE[1-4]p(12|34)')==1
elseif strcmp(vvs(1:2),'wE') & is14(vvs(3)) & vvs(4)=='p' & ...
(strcmp(vvs(5:6),'12') | strcmp(vvs(5:6),'34')) 
	cl_id = vs(3);
	inst = 'EFW';
	sig = 'E';
	sen = vs(4:6);
	var_units =  {'mV/m'};
	if CEF
		var_size = 1;
		var_name = {['E_' sen]};
		field_name = {['Electric field']};
		frame = {'component>sc_xy'};
		switch sen
		case 'p12'
			comp_desc = {'x>WEC Z axis (12)'};
		case 'p34'
			comp_desc = {'y>WEC Y axis (34)'};
		end
		si_conv = {'1.0e-3>V/m'};
		var_labels = {'E'};
	else
		frame = 'sc';
		var_labels = {['E' sen]};
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% spin fits E p12 and p34
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elseif regexp(vs,'^diEs[1-4]p(12|34)')==1
elseif strcmp(vvs(1:4),'diEs') & is14(vvs(5)) & vvs(6)=='p' & ...
(strcmp(vvs(7:8),'12') | strcmp(vvs(7:8),'34'))
	cl_id = vs(5);
	inst = 'EFW';
	sig = 'E';
	sen = vs(6:8);
	if CEF
		var_units =  {'mV/m'};
		if CAA
			% for CAA we save only Ex and Ey
			var_size = 2;
			var = var(1:3);
			frame = {'vector>dsi_xy'};
			var_label_1 = {'"x", "y"'};
		else
			var_size = 3;
			frame = {'vector>dsi_xyz'};
			var_label_1 = {'"x", "y", "z"'};
		end
		var_name = {['Es_' sen]};
		field_name = {['Electric field (spin fit ' sen ')']};
		si_conv = {'1.0e-3>V/m'};
		var_labels = {'E'};
	else
		sen = ['spin fits ' sen];
		frame = 'DSI,  Ez==0 : not measured';
		var_labels = {'Ex','Ey','Ez'};
		var_units =  {'mV/m','mV/m','mV/m'};
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full resolution E
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elseif regexp(vs,'^diE[1-4]p1234')==1
elseif (length(vs)==9 & strcmp(vs(1:3),'diE') & strcmp(vvs(5:9),'p1234'))
	cl_id = vs(4);
	inst = 'EFW';
	sig = 'E';
	sen = vs(5:9);
	if CEF
		var_units =  {'mV/m'};
		if CAA
			% for CAA we save only Ex and Ey
			var_size = 2;
			var = var(:,1:3);
			frame = {'vector>dsi_xy'};
			var_label_1 = {'"x", "y"'};
		else
			var_size = 3;
			frame = {'vector>dsi_xyz'};
			var_label_1 = {'"x", "y", "z"'};
		end
		var_name = {['E' sen]};
		field_name = {'Electric field'};
		si_conv = {'1.0e-3>V/m'};
		var_labels = {'E'};
	else
		frame = 'DSI,  Ez==0 : not measured';
		var_labels = {'Ex','Ey','Ez'};
		var_units =  {'mV/m','mV/m','mV/m'};
	end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% despun full resolution E with assumption E.B = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elseif regexp(vs,'^(diE[1-4]|diEs[1-4])$')==1
elseif (length(vs)==4 & strcmp(vvs(1:3),'diE') & is14(vvs(4))) | ...
(length(vs)==5 & strcmp(vvs(1:4),'diEs') & is14(vvs(5)))
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'E';
	if vs(4)=='s', v.sen = 's'; else, v.sen = ''; end
	v.cs = {'vector>dsi_xyz','scalar>na'};
 	v.units =  {'mV/m','deg'};
	v.si_conv = {'1.0e-3>V/m','1>degree'};
	v.size = [3 1];
	v.name = {['E' v.sen], 'Theta'};
	v.labels = {['P' v.sen], 'Theta'};
	v.label_1 = {'"x", "y", "z"',''};
	v.field_name = {'Electric field','Elevation of B above the sc spin plane'};
	v.com = com_Ez;
	v.file = 'mEdB';
	if vs(4)=='s', v.quant = 'edbs'; else, v.quant = 'edb'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% full resolution E in GSE coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elseif regexp(vs,'^(E[1-4]|Es[1-4])')==1
elseif (length(vs)==2 & vvs(1)=='E' & is14(vvs(2))) | ...
(length(vs)==3 & strcmp(vvs(1:2),'Es') & is14(vvs(3)))
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'E';
	if vs(2)=='s', v.sen = 's'; else, v.sen = ''; end
	v.cs = {'vector>gse_xyz','scalar>na'};
 	v.units =  {'mV/m','deg'};
	v.si_conv = {'1.0e-3>V/m','1>degree'};
	v.size = [3 1];
	v.name = {['E' v.sen], 'Theta'};
	v.labels = {['P' v.sen], 'Theta'};
	v.label_1 = {'"x", "y", "z"',''};
	v.field_name = {'Electric field','Elevation of B above the sc spin plane'};
	v.com = com_Ez;
	v.file = 'mEdB';
	if vs(2)=='s', v.quant = 'edbs'; else, v.quant = 'edb'; end
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%elseif regexp(vs,'^(diVExBs[1-4]|VExBs[1-4])')==1
elseif (((length(vs)==8 & vvs(7)=='s')|length(vs)==7) & ...
strcmp(vvs(1:6),'diVExB') & is14(vs(end))) | ...
(((length(vs)==6 & vvs(5)=='s')|length(vs)==5) & strcmp(vvs(1:4),'VExBs') & is14(vs(end)))
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
	v.name = {['V' v.sen], 'Theta'};
	v.labels = {'V', 'Theta'};
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
%elseif regexp(vs,'^NVps[1-4]')==1
elseif (length(vs)==5 & strcmp(vvs(1:4),'NVps') & is14(vvs(5)))
	v.cl_id = vs(end);
	v.inst = 'EFW';
	v.frame = 'sc';
	v.sig = 'P';
	v.sen = 'p1234';
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
%elseif regexp(vs,'^A[1-4]')
elseif length(vs)==2 & vs(1)=='A' & is14(vvs(2))
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
% CIS V PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^NC(h|p)[1-4]')

	if CAA
		c_log('fcal', ['Variable ' vs ' is not intended for the CAA'])
		CAA = 0;
	end
	if CEF
		c_log('fcal', ['CEF export is not (yet) supported for ' vs])
		CEF = 0;
	end
	cl_id = vs(4);
	inst = 'CIS PP';
	if vs(3)=='h'
		sig = 'N';
		sen = 'HIA';
	else
		sig = 'Np';
		sen = 'CODIF';
	end
	frame = '';
	var_labels = {'N'};
	var_units =  {'cc'};
	com = 'This data is CSDS PP';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIS N PP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^VC(h|p)[1-4]')

	if CAA
		c_log('fcal', ['Variable ' vs ' is not intended for the CAA'])
		CAA = 0;
	end
	if CEF
		c_log('fcal', ['CEF export is not (yet) supported for ' vs])
		CEF = 0;
	end
	cl_id = vs(4);
	inst = 'CIS PP';
	if vs(3)=='h'
		sig = 'V';
		sen = 'HIA';
	else
		sig = 'Vp';
		sen = 'CODIF';
	end
	frame = 'GSE';
	var_labels = {'Vx','Vy','Vz'};
	var_units =  {'km/s','km/s','km/s'};
	com = 'This data is CSDS PP';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dump without headers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(vs,'dump')
	if isempty(file_name), filename = inputname(1); end
	fid = fopen([file_name '.dat'],'w');
	for j=1:size(var,1)
		if var(j,1)>5e8 % assume first column time in isdat epoch
      d=sprintf('%4.0f %2.0f %2.0f %2.0f %2.0f %7.4f ',fromepoch(var(j,1)));
	  		fprintf(fid,[d num2str(var(j,2:end)) '\n']);
		else, fprintf(fid,[num2str(var(j,1:end)) '\n']);
		end
	end
	fclose(fid);
	return
else
	error('Wariable name not recognized')
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
