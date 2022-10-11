function c_export_ascii(var,varargin)
%C_EXPORT_ASCII export variable into ascii/CEF file with header
%
% c_export_ascii(var,[option,value])
%
% Options:
%  var_name - use another name (standart) for identification of a
%             custom named variable;
%             If var_name is set 'dump', the variable will be saved without
%             any headers.
%  comm     - add extra comment;
%  f_name   - specify file name. Dafault is {var_name_s}.dat
%  mode     - 'plain' (default), 'cef' (Cluster Exchange Format) or 'caa'
%
% Information written into the header of file is based on the variable name
% If comment is given it is added to the header
%
% Example:
%	c_export_ascii(diE3p1234)
%	c_export_ascii(Etemp,'var_name','diE3p1234')
%	c_export_ascii(Etemp,'var_name','diE3p1234',...
%      'comm','electric field is high pass filtered at 3Hz')
%
%   c_export_ascii(dE,'var_name','dump') - no headers
%
% $Id$

% ----------------------------------------------------------------------------
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

if nargin<1, help c_export_ascii; return, end
if nargin>2, have_options = 1; args = varargin;
else, have_options = 0;
end

CEF = 0;
CAA = 0;
comment = '';
vs = inputname(1);
file_name = '';

while have_options
  l = 2;
  if length(args)>1
    switch(args{1})
      case 'var_name'
        if ischar(args{2}), vs = args{2};
        else, irf_log('fcal','wrong ArgType : var_name must be string')
        end
      case 'comm'
        if ischar(args{2}), comment = args{2};
        else, irf_log('fcal','wrong ArgType : comm must be string')
        end
      case 'f_name'
        if ischar(args{2}), file_name = args{2};
        else, irf_log('fcal','wrong ArgType : f_name must be string')
        end
      case 'mode'
        if ischar(args{2})
          switch(args{2})
            case 'plain'
            case 'cef'
              CEF = 1;
              DATA_VERSION = '01';
            case 'caa'
              CAA = 1;
              CEF = 1;
              DATA_VERSION = '01';
            otherwise
              irf_log('fcal',...
                'wrong ArgType : mode must be one of : ''plain'', ''cef'' or ''caa''')
              irf_log('fcal','Using default mode : ''plain''');
          end
        else, irf_log('fcal',...
            'wrong ArgType : mode must be one of : ''plain'', ''cef'' or ''caa''')
        end
      otherwise
        irf_log('fcal',['Option ''' args{1} '''not recognized'])
    end
    if length(args) > l, args = args(l+1:end);
    else, break
    end
  else
    error('caa:wrongArgType','use c_export_ascii(..,''option'',''value'')')
  end
end

com_Ez = 'Ez is not reliable when magnetic field is close to the spin plane\n%% The last column shows the angle of B with respect to the spin plane (B,spin)';
com = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if regexp(vs,'^P[1-4]$')==1
  
  cl_id = vs(2);
  inst = 'EFW';
  sig = 'P';
  sen = '1234';
  var_units =  {'V'};
  if CEF
    var_size = 1;
    var_name = {['P_' sen]};
    field_name = {'Averaged probe potential from all probes'};
    frame = {'scalar>na'};
    si_conv = {''};
    var_labels = {'P'};
  else
    frame = 'sc';
    var_labels = {['P' sen]};
    com = 'this signal is averaged from all probes available at the time';
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % P - individual probes
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^P10Hz[1-4]p[1-4]$')==1
  cl_id = vs(6);
  inst = 'EFW';
  sig = 'P';
  sen = vs(8);
  var_units =  {'V'};
  if CEF
    var_size = 1;
    var_name = {['P_' sen]};
    field_name = {['Probe #' sen ' potential']};
    frame = {'scalar>na'};
    si_conv = {''};
    var_labels = {'P'};
  else
    frame = 'sc';
    var_labels = {['P' sen]};
  end
elseif regexp(vs,'^P32kHz[1-4]p[1-4]$')==1
  cl_id = vs(7);
  inst = 'EFW';
  sig = 'P';
  sen = vs(8);
  var_units =  {'V'};
  frame = 'sc';
  var_labels = {['P' sen]};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % raw E p12 and p34
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^wE[1-4]p(12|32|34|24)')==1
  
  cl_id = vs(3);
  inst = 'EFW';
  sig = 'E';
  sen = vs(4:6);
  var_units =  {'mV/m'};
  if CEF
    var_size = 1;
    var_name = {['E_' sen]};
    field_name = {'Electric field'};
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
elseif regexp(vs,'^(i)?diEs[1-4]p(12|32|34)')==1
  
  
  cl_id = vs(5);
  inst = 'EFW';
  sig = 'E';
  sen = vs(6:8);
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
elseif regexp(vs,'^di(b)?E[1-4]p1234')==1
  
  desc = c_desc(vs);
  cl_id = desc.cl_id;
  inst = desc.inst;
  sig = desc.sig;
  sen = desc.sen;
  
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
elseif regexp(vs,'^(i)?(diE[1-4]|diEs[1-4])$')==1
  vs_old = vs;
  if vs(1)=='i'
    frame = 'DSI (inertial frame),  Ez is computed from E.B=0';
    vs = vs(2:end);
  else
    frame = 'DSI,  Ez is computed from E.B=0';
  end
  
  cl_id = vs(end);
  inst = 'EFW';
  sig = 'E';
  if CEF
    if vs(4)=='s', sen = 's'; else, sen = ''; end
    var_units =  {'mV/m'};
    var_name = {['E' sen]};
    field_name = {'Electric field'};
    si_conv = {'1.0e-3>V/m'};
    var_labels = {'E'};
    var_label_1 = {'"x", "y", "z"'};
    frame = {'vector>dsi_xyz'};
    if CAA
      % for CAA we save only E
      var_size = 3;
      var = var(:,1:4);
    else
      var_size = [3 1];
      var_name = [var_name {'B_field_elevation'}];
      frame = [frame {'scalar>na'}];
      field_name = [field_name{:} {'Elevation of B above the sc spin plane'}];
      var_labels = [var_labels{:} {'Theta'}];
      var_label_1 = [var_label_1{:} {''}];
      si_conv = [si_conv {'1>degree'}];
      var_units =  [var_units {'deg'}];
    end
  else
    if vs(4)=='s', sen = 'spin fit';
    else, sen = 'p1234';
    end
    var_labels = {'Ex','Ey','Ez','(B,spin)','sdev'};
    var_units =  {'mV/m','mV/m','mV/m','deg','unitless'};
    com = com_Ez;
  end
  vs = vs_old;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % B from FGM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?B(r)?[1-4]$')==1
  
  if CAA
    irf_log('fcal', ['Variable ' vs ' is not intended for the CAA'])
    CAA = 0;
  end
  
  cl_id = vs(end);
  inst = 'FGM';
  sig = 'B';
  sen = '';
  if strcmp(vs(1:2),'di'), frame = 'DSI,  approximately the same as GSE';
  else, frame = 'GSE';
  end
  var_labels = {'Bx','By','Bz'};
  var_units =  {'nT','nT','nT'};
  if vs(2)=='r' || (length(vs)>3 && vs(4)=='r')
    com = 'B data is interpolated to E';
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % full resolution E in GSE coordinates
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?(E[1-4]|Es[1-4])$')==1
  
  vs_old = vs;
  if vs(1)=='i'
    frame = 'GSE (inertial frame)';
    vs = vs(2:end);
  else
    frame = 'GSE';
  end
  
  if CAA
    irf_log('fcal', ['Variable ' vs ' is not intended for the CAA'])
    CAA = 0;
  end
  cl_id = vs(2);
  inst = 'EFW';
  sig = 'E';
  if CEF
    if vs(2)=='s', sen = 's'; else, sen = ''; end
    var_units =  {'mV/m'};
    var_name = {['E' sen]};
    field_name = {'Electric field'};
    si_conv = {'1.0e-3>V/m'};
    var_labels = {'E'};
    var_label_1 = {'"x", "y", "z"'};
    frame = {'vector>gse_xyz'};
    if CAA
      % for CAA we save only E
      var_size = 3;
      var = var(:,1:4);
    else
      var_size = [3 1];
      var_name = {var_name{:}, 'B_field_elevation'};
      frame = {frame{:}, 'scalar>na'};
      field_name = {field_name{:}, 'Elevation of B above the sc spin plane'};
      var_labels = {var_labels{:}, 'Theta'};
      var_label_1 = {var_label_1{:}, ''};
      si_conv = {si_conv{:}, '1>degree'};
      var_units =  {var_units{:}, 'deg'};
    end
  else
    if vs(2)=='s', sen = 'spin fit';
    else, sen = 'p1234';
    end
    var_labels = {'Ex','Ey','Ez','(B,spin)'};
    var_units =  {'mV/m','mV/m','mV/m','deg'};
    com = com_Ez;
  end
  vs = vs_old;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % ExB
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(diVExB(s)?[1-4]|VExB(s)?[1-4])$')==1
  
  cl_id = vs(end);
  inst = 'EFW';
  sig = 'V=ExB';
  sen = 'spin fit';
  if strcmp(vs(1:2),'di'), frame = 'DSI,  approximately the same as GSE';
  else
    frame = 'GSE';
    if CAA
      irf_log('fcal',...
        ['Variable ' vs ' is not intended for the CAA'])
      CAA = 0;
    end
  end
  if CEF
    var_units =  {'km/s'};
    var_name = {'EsxB'};
    field_name = {'Convection velocity'};
    si_conv = {'1.0e3>m/s'};
    var_labels = {'V'};
    var_label_1 = {'"x", "y", "z"'};
    frame = {'vector>dsi_xyz'};
    if CAA
      % for CAA we save only V
      var_size = 3;
      var = var(:,1:4);
    else
      var_size = [3 1];
      var_name = {var_name{:}, 'B_field_elevation'};
      frame = {frame{:}, 'scalar>na'};
      field_name = {field_name{:}, 'Elevation of B above the sc spin plane'};
      var_labels = {var_labels{:}, 'Theta'};
      var_label_1 = {var_label_1{:}, ''};
      si_conv = {si_conv{:}, '1>degree'};
      var_units =  {var_units{:}, 'deg'};
    end
  else
    var_labels = {'Vx','Vy','Vz','(B,spin)'};
    var_units =  {'km/s','km/s','km/s','deg'};
    com = com_Ez;
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % full/4 sec resolution satellite potential
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^P(s)?[1-4]$')==1
  
  desc = c_desc(vs);
  cl_id = desc.cl_id;
  inst = desc.inst;
  sig = desc.sig;
  sen = desc.sen;
  frame = desc.frame;
  var_labels = desc.labels;
  var_units =  desc.units;
  com = desc.com;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % full resolution satellite potential and derived density
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(b)?NVps[1-4]$')==1
  
  desc = c_desc(vs);
  cl_id = desc.cl_id;
  inst = desc.inst;
  sig = desc.sig;
  sen = desc.sen;
  frame = desc.frame;
  var_labels = desc.labels;
  var_units =  desc.units;
  com = desc.com;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % phase
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^A(two)?[1-4]$')
  
  if CAA
    irf_log('fcal', ['Variable ' vs ' is not intended for the CAA'])
    CAA = 0;
  end
  if CEF
    irf_log('fcal', ['CEF export is not (yet) supported for ' vs])
    CEF = 0;
  end
  cl_id = vs(end);
  inst = 'ephemeris';
  sig = 'phase';
  sen = '';
  frame = 'SC';
  var_labels = {'phase'};
  var_units =  {'deg'};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % CIS V PP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^NC(h|p)[1-4]$')
  
  if CAA
    irf_log('fcal', ['Variable ' vs ' is not intended for the CAA'])
    CAA = 0;
  end
  if CEF
    irf_log('fcal', ['CEF export is not (yet) supported for ' vs])
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
    irf_log('fcal', ['Variable ' vs ' is not intended for the CAA'])
    CAA = 0;
  end
  if CEF
    irf_log('fcal', ['CEF export is not (yet) supported for ' vs])
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
  % R
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^R?[1-4]$')==1
  
  cl_id = vs(end);
  inst = 'EPHEMERIS';
  sig = 'R';
  sen = '';
  frame = 'GSE';
  var_labels = {'Rx','Ry','Rz'};
  var_units =  {'km','km','km'};
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % dump without headers
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(vs,'dump')
  if isempty(file_name), file_name = inputname(1); end
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
  irf_log('fcal','Variable name not recognized, will do nothing.')
  return
end

sz = size(var);
n_data = sz(2) - 1; % number of data columns - time

if CEF
  % determine the sampling frequency
  dt = var(2,1) - var(1,1);
  s_f = [.25 1 5 22.5 25 67.5 450 9000 18000];
  s_fr = '';
  ee = .05;
  for j=1:length(s_f)
    if (1/dt > s_f(j)*(1-ee)) && (1/dt < s_f(j)*(1+ee))
      s_fr = s_f(j);
      break
    end
  end
  if isempty(s_fr), s_fr = 1/dt; end
  resolution = sprintf('%.2f Hz',s_fr);
  delta = sprintf('%.2f',.5/s_fr);
  irf_log('save',[var_name{1} ' resolution: ' resolution ' Delta: ' delta])
  mask = '';
  for j=1:n_data, mask = [mask ', %8.3f']; end
  time_s = epoch2iso(var(:,1));
else
  t0=toepoch(fromepoch(var(1,1)).*[1 1 1 0 0 0]);
  t0_s=datestr(datenum(fromepoch(t0)),0);
  var(:,1) = var(:,1) - t0;
  var_s =    'time        ';
  var_unit = '[s]         ';
  mask = '%11.6f ';
  for j=1:n_data
    var_length_nine=strvcat(var_labels{j},'         ');
    var_unit_length_nine=strvcat(['[' var_units{j} ']'],'         ');
    var_s = [var_s var_length_nine(1,:) ];
    var_unit = [var_unit var_unit_length_nine(1,:)];
    mask = [mask '%8.3f '];
  end
end

% construct filename
if CEF, ext_s = '.cef'; else, ext_s = '.dat'; end
if isempty(file_name)
  if CAA % we have special names for CAA
    file_name = ...
      ['C' num2str(cl_id) '_CAA_EFW_' var_name{1} '_' irf_fname(var(1,1)) '_V' DATA_VERSION];
  else, file_name = vs;
  end
end
fid = fopen([file_name ext_s],'w');
if CEF
  fprintf(fid,'!-------------------- CEF ASCII FILE --------------------|\n');
  nnow = now;
  fprintf(fid,['! created on ' datestr(nnow) '\n']);
  fprintf(fid,'!--------------------------------------------------------|\n');
  fprintf(fid,['FILE_NAME = "' file_name ext_s '"\n']);
  fprintf(fid,'FILE_TYPE_VERSION = "CEF-2.0"\n!\n');
  fprintf(fid,'START_META = Project\n');
  fprintf(fid,'ENTRY = "PROJ>CLUSTER II"\n');
  fprintf(fid,'END_META = Project\n!\n');
  fprintf(fid,'START_META = Discipline\n');
  fprintf(fid,'ENTRY = "SPACE PHYSICS> MAGNETOSPHERIC PHYSICS"\n');
  fprintf(fid,'END_META = Discipline\n!\n');
  fprintf(fid,'START_META = Data_type\n');
  fprintf(fid,['ENTRY = "RES>' resolution '"\n']);
  fprintf(fid,'END_META = Data_type\n!\n');
  fprintf(fid,'START_META = Descriptor\n');
  fprintf(fid,['ENTRY = "INS>' inst '"\n']);
  fprintf(fid,'END_META = Descriptor\n!\n');
  fprintf(fid,'START_META = Data_version\n');
  fprintf(fid,['ENTRY = "' DATA_VERSION '"\n']);
  fprintf(fid,'END_META = Data_version\n!\n');
  fprintf(fid,'START_META = Generation_date\n');
  ss = sprintf('%6f',nnow-fix(nnow));
  ss = [datestr(nnow,'yyyy-mm-ddTHH:MM:SS') ss(2:8) 'Z'];
  fprintf(fid,['ENTRY = "' ss '"\n']);
  fprintf(fid,'END_META = Generation_date\n!\n');
  fprintf(fid,'START_META = Caveats\n');
  fprintf(fid,['ENTRY = "' com '"\n']);
  fprintf(fid,'END_META = Caveats\n!\n');
  fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
  fprintf(fid,'!                                        Variables    !\n');
  fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
  fprintf(fid,'START_VARIABLE    = time_tags\n');
  fprintf(fid,'  VALUE_TYPE      = ISO_TIME\n');
  fprintf(fid,['  DELTA_PLUS      = ' delta '\n']);
  fprintf(fid,['  DELTA_MINUS     = ' delta '\n']);
  fprintf(fid,'  LABLAXIS        = "UT"\n');
  fprintf(fid,'  FIELDNAM        = "Universal Time"\n');
  fprintf(fid,'END_VARIABLE      = time_tags\n\n');
  
  for j=1:length(var_size)
    fprintf(fid,['START_VARIABLE    = ' var_name{j} '\n']);
    if var_size(j) > 1
      fprintf(fid,'  SIZES           = %d\n',var_size(j));
    end
    fprintf(fid,'  VALUE_TYPE      = FLOAT\n');
    fprintf(fid,['  FIELDNAM        = "' field_name{j} '"\n']);
    fprintf(fid,['  FRAME           = "' frame{j} '"\n']);
    if exist('comp_desc','var')
      if ~isempty(comp_desc{j})
        fprintf(fid,['  COMPONENT_DESC  = "' comp_desc{j} '"\n']);
      end
    end
    if ~isempty(si_conv{j})
      fprintf(fid,['  SI_CONVERSION   = "' si_conv{j} '"\n']);
    end
    fprintf(fid,['  UNITS           = "' var_units{j} '"\n']);
    %FILLVAL         = -1.0E-10
    fprintf(fid,['  LABLAXIS        = "' var_labels{j} '"\n']);
    if var_size(j) > 1
      fprintf(fid,['  LABEL_1         = ' var_label_1{j} '\n']);
    end
    fprintf(fid,'  DEPEND_0        = time_tags\n');
    fprintf(fid,['END_VARIABLE      = ' var_name{j} '\n\n']);
  end
  fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
  fprintf(fid,'!                                        Data         !\n');
  fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
  fprintf(fid,'DATA_UNTIL = EOF\n');
  for j=1:length(var(:,1))
    fprintf(fid,time_s(j,:));
    fprintf(fid,[mask '\n'],var(j,2:end)');
  end
else
  fprintf(fid,['%% this file was created on ' char(datetime("now","Format","dd-MMM-uuuu HH:mm:ss")) ' \n%%\n']);
  fprintf(fid,['%% SC:        Cluster ' cl_id ' \n']);
  fprintf(fid,['%% Intrument: ' inst ' \n']);
  fprintf(fid,['%% Signal:    ' sig ' \n']);
  fprintf(fid,['%% Sensor:    ' sen ' \n']);
  fprintf(fid,['%% Coord Sys: ' frame ' \n%%\n']);
  fprintf(fid,['%% comment:   ' com ' \n%%\n']);
  if comment
    comm=tokenize(comment,'\n');
    for j=1:size(comm,2), fprintf(fid,['%%   ' comm{1,j} ' \n']); end
  end
  fprintf(fid,['%%\n%% Time from: ' t0_s ' \n%%\n']);
  fprintf(fid,['%% ' var_s ' \n']);
  fprintf(fid,['%% ' var_unit ' \n']);
  fprintf(fid,[mask '\n'],var');
end

fclose(fid);
