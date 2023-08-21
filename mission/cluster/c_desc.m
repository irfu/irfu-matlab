function varargout = c_desc(vs,v_info)
%C_DESC  Provide a description of a Cluster variable
%
% C_DESC(V_S [,V_S_INFO])
%		prints out a description of variable VS
%
% DESC = C_DESC(V_S [,V_S_INFO])
%		returns a description of variable VS as structure DESC.
%		returns empty cell string if variable name unknown
%
% Input:
%	V_S      - string defining a variable
%   V_S_INFO - infor structure for the variable
%
% Output:
%	DESC - structure containing a description of variable VS. It has the
%	following fields:
%   cl_id		%Cluster ID
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
%   ptype       %Parameter type of a CEF variable
%   valtype     %Value type of a CEF variable
%   sigdig      %Significant digits of a CEF variable
%   com			%Comment
%   file		%Matlab file name (mXXX.mat)
%   quant		%Quantity name to use with getData
%   lev         %Data level: 0 - ClusterDB, 1 - ClusterProc, 2 - manual
%
% Examples:
% c_desc('diE2')
%	desc = c_desc('diE2');

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(1,2)
if ~ischar(vs), error('VS must be a string'), end
if nargin<2, v_info = []; end

com_Ez = 'Ez is not reliable when magnetic field is close to the spin plane';

vvs = 'XXXXXXXXXX';
vvs(1:length(vs)) = vs;

v.file_old = ''; % compatibility mode for old files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% P Spacecraft potential level 2/3 caa_export only no variable in .mat
%% files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if regexp(vs,'^P[1-4]$')==1
  v.data = 1;
  v.cl_id = vs(2);
  v.lev = 2;
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'P';
  v.sen = '1234';
  v.cs = {'na','na','na','na','na'};
  v.rep = {'scalar','scalar','scalar','scalar','scalar'};
  v.units =  {'V','na','na','na','na'};
  v.si_conv = {'1>V','','','',''};
  v.size = [ 1 1 1 1 1 ];
  v.name = {'Spacecraft_Potential','Probe','ASPOC_Status','P_Bitmask','P_Quality'};
  v.quant = 'p';
  v.labels = {'-Sc pot','Probe','ASPOC Active','Bitmask','Quality'};
  v.field_name = {'Spacecraft potential','Probe','ASPOC','Bitmask','Quality'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT','INT','INT','INT','INT'};
  v.sigdig = [ 6 4 1 5 1 ];
  v.ent = {'Instrument'};
  v.prop = {'Probe_Potential'};
  v.fluc = {'Waveform'};
  v.com = '';
  v.file = 'mP';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% P internal burst level 2 caa_export only no variable in .mat
  %% files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^PB[1-4]$')==1
  v.data = 1;
  v.cl_id = vs(2);
  v.lev = 2;
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'P';
  v.sen = '1234';
  v.cs = {'na','na','na','na','na'};
  v.rep = {'scalar','scalar','scalar','scalar','scalar'};
  v.units =  {'V','na','na','na','na'};
  v.si_conv = {'1>V','','','',''};
  v.size = [ 1 1 1 1 1 ];
  v.name = {'Spacecraft_Potential','Probe','ASPOC_Status','P_Bitmask','P_Quality'};
  v.quant = 'p';
  v.labels = {'-Sc pot','Probe','ASPOC Active','Bitmask','Quality'};
  v.field_name = {'Spacecraft potential','Probe','ASPOC','Bitmask','Quality'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT','INT','INT','INT','INT'};
  v.sigdig = [ 6 4 1 5 1 ];
  v.ent = {'Instrument'};
  v.prop = {'Probe_Potential'};
  v.fluc = {'Waveform'};
  v.com = '';
  v.file = 'cef export only';  % caa_export_new() cef export only
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% P & Ps
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(b)?P(s)?[1-4]$')==1
  v.data = 1;
  if vs(2)=='s' || vs(1)=='b' , v.cl_id = vs(3);
  else, v.cl_id = vs(2);
  end
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'P';
  if ~isempty(v_info) && isfield(v_info,'probe'), v.sen = num2str(v_info.probe);
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
    v.field_name = {'Spacecraft potential (4 sec resolution)'};
  elseif vs(1)=='b'
    v.quant = 'pburst';
    v.field_name = {'Spacecraft potential (internal burst)'};
  else
    v.quant = 'p';
    v.field_name = {'Spacecraft potential'};
  end
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.ent = {'Instrument'};
  v.prop = {'Probe_Potential'};
  v.fluc = {'Waveform'};
  v.com = ['this signal is averaged from probes ' v.sen];
  if vs(1)=='b', v.file = 'mEFWburst';
  else, v.file = 'mP';
  end
  v.lev = 1;
elseif regexp(vs,'^(b)?P(s)?[1-4]_info$')==1
  v.data = 0;
  if vs(2)=='s' || vs(1)=='b', v.cl_id = vs(3);
  else, v.cl_id = vs(2);
  end
  v.inst = 'EFW';
  v.com = 'Spacecraft potential INFO';
  if vs(1)=='b', v.file = 'mEFWburst';
  else, v.file = 'mP';
  end
  if vs(2)=='s'
    v.quant = 'ps';
  elseif vs(1)=='b'
    v.quant = 'pburst';
  else
    v.quant = 'p';
  end
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% P - individual probes
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
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.ent = {'Instrument'};
  v.prop = {'Probe_Potential'};
  v.fluc = {'Waveform'};
  v.com = '';
  v.file = 'mPR';
  v.quant = 'p';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% P - individual probes from internal burst
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif any(regexp(vs,'^P(32|4)kHz[1-4]p[1-4]$')==1) || ...
    any(regexp(vs,'^P180Hz[1-4]p[1-4]$')==1)
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
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mEFWburstR';
  v.file_old = 'mEFWburst';
  v.quant = 'pburst';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Raw and corrected E p12 and p34
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif any(regexp(vs,'^w(b|c)?E[1-4]p(12|32|34|42)$'))
  v.data = 1;
  v.inst = 'EFW';
  v.frame = 'na';
  v.sig = 'E';
  v.sen = vs(end-1:end);
  v.cs = {'SC'};
  v.rep = {'scalar'};
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.size = 1;
  v.name = {['P' v.sen]};
  v.labels = v.name;
  v.field_name = {['Electric field component measured between the probes '...
    v.sen(1) ' and ' v.sen(2)]};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.ent = {'Electric_Field'};
  v.prop = {'Component'};
  v.fluc = {'Waveform'};
  if vs(2)=='E'
    v.cl_id = vs(3);
    v.file = 'mER';
    v.com = '';
    v.lev = 0;
    v.quant = 'e';
  elseif vs(2)=='c'
    v.cl_id = vs(4);
    v.file = 'mERC';
    v.com = 'This data is not original raw data. It has been cleaned.';
    v.lev = 1;
    v.quant = 'ec';
  else
    v.cl_id = vs(4);
    v.file = 'mEFWburstR';
    v.com = 'This data is from EFW internal burst.';
    v.lev = 0;
    v.quant = 'e';
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Raw and corrected E p12 and p34
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif any(regexp(vs,'^wE8kHz[1-4]p(12|34)$'))
  v.data = 1;
  v.inst = 'EFW';
  v.frame = 'na';
  v.sig = 'E';
  v.sen = vs(end-1:end);
  v.cs = {'SC'};
  v.rep = {'scalar'};
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.size = 1;
  v.name = {['P' v.sen]};
  v.labels = v.name;
  v.field_name = {['Electric field component measured between the probes '...
    v.sen(1) ' and ' v.sen(2)]};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.ent = {'Electric_Field'};
  v.prop = {'Component'};
  v.fluc = {'Waveform'};
  v.cl_id = vs(7);
  v.file = 'mEFWburstR';
  v.com = 'This data is from EFW internal burst.';
  v.lev = 0;
  v.quant = 'e';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Wake description
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^WAKE[1-4]p(12|32|34|42)$')==1
  v.data = 1;
  v.cl_id = vs(5);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'WAKE';
  v.sen = vs(end-1:end);
  v.cs = {'ISR2','na','na'};
  v.rep = {'scalar','scalar','scalar'};
  v.units =  {'deg','mV/m','deg'};
  v.si_conv = {'1>degree','1.0e-3>V m^-1','1>degree'};
  v.size = [1 1 1];
  v.name = {['Wake-p' v.sen ' location'], ['Wake-p' v.sen ' amp'],...
    ['Wake-p' v.sen ' h-width']};
  v.labels = v.name;
  v.field_name = {'Wake location','Wake amplitude','Wake half-width'};
  v.ptype = {'Data','Data','Data'};
  v.valtype = {'FLOAT','FLOAT','FLOAT'};
  v.sigdig = [6 6 6];
  v.com = '';
  v.file = 'mERC';
  v.quant = 'ec';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Plasmaspheric/Lobe Wake
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(NONSIN|PS|LO)WAKE[1-4]p(12|32|34|42)$')==1
  v.data = 1;
  v.cl_id = vs(7);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'WAKE';
  v.sen = vs(end-1:end);
  v.cs = {'na'};
  v.rep = {'scalar'};
  v.units =  {'sec'};
  v.si_conv = {'1>sec'};
  v.size = 1;
  if vs(1)=='P', reg = 'Plasmaspheric';
  elseif vs(1)=='N', reg = 'Non-sinusoidal';
  else, reg = 'Lobe';
  end
  v.name = {[reg ' wake-p' v.sen]};
  v.labels = {['Wake-p' v.sen ' stop']};
  v.field_name = {'Wake stop'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mEFW';
  v.quant = 'wake';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Plasmaspheric/Lobe Wake
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^HBSATDSC[1-4]p(12|32|34)$')==1
  v.data = 1;
  v.cl_id = vs(9);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'HBSATDSC';
  v.sen = vs(end-1:end);
  v.cs = {'na','na'};
  v.rep = {'scalar','scalar'};
  v.units =  {'mV/m','deg'};
  v.si_conv = {'1.0e-3>V m^-1','1>degrees'};
  v.size = [1 1];
  v.name = {'HBSat_Max','HBSat_Width'};
  v.field_name = {'HBSat Max','HBSat Width'};
  v.labels = v.field_name;
  v.ptype = {'Data','Data'};
  v.valtype = {'FLOAT','FLOAT'};
  v.sigdig = 3;
  v.com = '';
  v.file = 'mEFW';
  v.quant = 'hbiassa';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% RSPEC
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^RSPEC[1-4]p(12|32|34|42)$')
  v.data = 1;
  v.cl_id = vs(6);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'E-RSPEC';
  v.sen = vs(8:9);
  v.cs = {'ISR2', 'ISR2','ISR2', 'ISR2','ISR2'};
  v.units =  {'mV/m','mV/m','mV/m','mV/m','mV/m'};
  v.si_conv = {'1.0e-3>V m^-1','1.0e-3>V m^-1','1.0e-3>V m^-1',...
    '1.0e-3>V m^-1','1.0e-3>V m^-1'};
  v.size = [2 2 2 2 2];
  v.name = {'ER_1omega','ER_2omega','ER_3omega','ER_4omega','ER_5omega'};
  v.labels = v.name;
  v.label_1 = {'"E1w_x", "E1w_y"','"E2w_x", "E2w_y"','"E3w_x", "E3w_y"',...
    '"E4w_x", "E4w_y"','"E5w_x", "E5w_y"'};
  v.col_labels = {{'x','y'},{'x','y'},{'x','y'},{'x','y'},{'x','y'}};
  v.field_name = {'Raw electric field : 1 omega',...
    'Raw electric field : 2 omega','Raw electric field : 3 omega',...
    'Raw electric field : 4 omega','Raw electric field : 5 omega'};
  v.ptype = {'Data','Data','Data','Data','Data'};
  v.valtype = {'FLOAT','FLOAT','FLOAT','FLOAT','FLOAT'};
  v.sigdig = [6 6 6 6 6];
  v.file = 'mEFW';
  v.quant = 'rawspec';
  v.com = '';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spin fits E p12 and p34
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?diEs[1-4]p(12|32|34)$')==1
  v.data = 1;
  if vvs(1)=='i'
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
  %	v.cs = {'ISR2','na','na','na'};
  %	v.rep = {'xy','','',''};
  %	v.units =  {'mV/m','mV/m','unitless','unitless'};
  %	v.si_conv = {'1.0e-3>V m^-1','1.0e-3>V m^-1','',''};
  %	v.size = [3 1 1 1];
  %	v.tensor_order = [1 0 0 0];
  %	v.name = {'E_Vec_xy_ISR2', 'E_sigma','E_bitmask','E_quality'};
  %	v.labels = {'E','St dev','Bitmask','Quality'};
  %	v.label_1 = {'"Ex", "Ey"','','',''};
  %	v.col_labels = {{'x','y','z'},'','',''};
  %	v.rep_1 = {'"x", "y"','','',''};
  %	v.field_name = {'Electric field (4 sec resolution)',...
  %		'Electric field standard deviation',...
  %		'Electric field measurement quality bitmask',...
  %		'Electric field measurement quality flag (9=best)'};
  %	v.ptype = {'Data','Data','Support_Data','Support_Data'};
  %	v.valtype = {'FLOAT','FLOAT','INT','INT'};
  %	v.sigdig = [6 6 5 1];
  %	v.ent = {'Electric_Field','Electric_Field','Electric_Field','Electric_Field'};
  %	v.prop = {'Vector','Vector','Status','Status'};
  %	v.fluc = {'Waveform','Fluctuation_Level','',''};
  v.cs = {'ISR2','na'};
  v.rep = {'xy',''};
  v.units =  {'mV/m','mV/m'};
  v.si_conv = {'1.0e-3>V m^-1','1.0e-3>V m^-1'};
  v.size = [3 1];
  v.tensor_order = [1 0];
  v.name = {'E_Vec_xy_ISR2', 'E_sigma'};
  v.labels = {'E','St dev'};
  v.label_1 = {'"Ex", "Ey"',''};
  v.col_labels = {{'x','y','z'},''};
  v.rep_1 = {'"x", "y"',''};
  v.field_name = {'Electric field (4 sec resolution)',...
    'Electric field standard deviation'};
  v.ptype = {'Data','Data'};
  v.valtype = {'FLOAT','FLOAT'};
  v.sigdig = [6 6];
  v.ent = {'Electric_Field','Electric_Field'};
  v.prop = {'Vector','Vector'};
  v.fluc = {'Waveform','Fluctuation_Level'};
  v.com = 'Ez=0 by definition (not measured).';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spin fits E p12 and p34
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?diELXs[1-4]p(12|32|34|42)$')==1
  v.data = 1;
  if vvs(1)=='i'
    vvs = vvs(2:end);
    v.frame = 'inertial';
    v.file = 'mEIDSI';
    v.quant = 'idies';
  else
    v.frame = 'sc';
    v.file = 'mEDSI';
    v.quant = 'dies';
  end
  v.cl_id = vvs(7);
  v.inst = 'EFW';
  v.sig = 'E';
  v.sen = vvs(9:10);
  v.cs = {'ISR2','na'};
  v.rep = {'xy',''};
  v.units =  {'mV/m','mV/m'};
  v.si_conv = {'1.0e-3>V m^-1','1.0e-3>V m^-1'};
  v.size = [3 1];
  v.tensor_order = [1 0];
  v.name = {'E_Vec_xy_ISR2', 'E_sigma'};
  v.labels = {'E','St dev'};
  v.label_1 = {'"Ex", "Ey"',''};
  v.col_labels = {{'x','y','z'},''};
  v.rep_1 = {'"x", "y"',''};
  v.field_name = {'Electric field (4 sec resolution)',...
    'Electric field standard deviation'};
  v.ptype = {'Data','Data'};
  v.valtype = {'FLOAT','FLOAT'};
  v.sigdig = [6 6];
  v.ent = {'Electric_Field','Electric_Field'};
  v.prop = {'Vector','Vector'};
  v.fluc = {'Waveform','Fluctuation_Level'};
  v.com = 'Ez=0 by definition (not measured).';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spin fits E p12/32 and p34 level 3 caa_export only no variable in .mat
  %% files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^SFIT[1-4]$')==1
  v.data = 1;
  v.frame = 'sc';
  v.file = 'cef export only';  % caa_export_new() cef export only
  v.quant = 'dies';
  v.cl_id = vvs(5);
  v.inst = 'EFW';
  v.sig = 'E';
  v.sen = '';
  v.cs = {'ISR2','na','ISR2','na'};
  v.units =  {'mV/m','mV/m','mV/m','mV/m'};
  v.si_conv = {'1.0e-3>V m^-1','1.0e-3>V m^-1','1.0e-3>V m^-1',...
    '1.0e-3>V m^-1'};
  v.size = [2 1 2 1];
  v.name = {'E_Vec_xy_ISR2','E_sigma','E_Vec_xy_ISR2','E_sigma'};
  v.labels = {'E','E','St dev','E','E','St dev'};
  v.field_name = {'Electric field (4 sec resolution)',...
    'Electric field standard deviation',...
    'Electric field (4 sec resolution)',...
    'Electric field standard deviation'};
  v.ptype = {'Data','Data','Data','Data'};
  v.valtype = {'FLOAT','FLOAT','FLOAT','FLOAT'};
  v.sigdig = [6 6 6 6];
  v.ent = {'Electric_Field','Electric_Field',...
    'Electric_Field','Electric_Field'};
  v.prop = {'Vector','Vector','Vector','Vector'};
  v.fluc = {'Waveform','Fluctuation_Level','Waveform',...
    'Fluctuation_Level'};
  v.com = ''; % Set in caa_export_new()
  v.lev = 2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% IB Internal Burst P1-4 BX-Z P12/P34 level 1 caa_export only no variable in .mat
  %% files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^IB[1-4]$')==1
  v.data = 1;
  v.frame = 'sc';
  v.file = 'cef export only';  % caa_export_new() cef export only
  v.quant = 'pburst';
  v.cl_id = vvs(3);
  v.inst = 'EFW';
  v.sig = 'P';
  v.sen = ''; % No fixed data order
  v.cs = {'na','na','na','na','na','na','na','na'};
  v.units = {'na','na','na','na','na','na','na','na'};
  v.si_conv = {'na','na','na','na','na','na','na','na'};
  v.size = [1 1 1 1 1 1 1 1];
  v.name = {'na','na','na','na','na','na','na','na'};
  v.labels = {'na','na','na','na','na','na','na','na'};
  v.field_name = {'na','na','na','na','na','na','na','na'};
  v.ptype = {'Data','Data','Data','Data','Data','Data','Data','Data'};
  v.valtype = {'INT','INT','INT','INT','INT','INT','INT','INT'};
  v.sigdig = [6 6 6 6 6 6 6 6];
  v.ent = {'Instrument'};

  v.com = ''; % Set in caa_export_new()
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spin fits of 2 omega for p32
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^w2W[1-4]p32$')==1
  v.data = 1;

  v.frame = 'sc';
  v.file = 'mEDSI';
  v.quant = 'dies';
  v.cl_id = vvs(4);
  v.inst = 'EFW';
  v.sig = 'E';
  v.sen = 32;
  v.cs = 'ISR2';
  v.rep = 'xy';
  v.units =  'mV/m';
  v.si_conv = '1.0e-3>V m^-1';
  v.size = 2;
  v.name = 'E_2omega_ISR2';
  v.labels = 'E2w';
  v.label_1 = '"Ex", "Ey"';
  v.col_labels = {'x','y'};
  v.rep_1 = '"x", "y"';
  v.field_name = 'Electric field 2 omega';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% E internal burst level 2 caa_export only no variable in .mat
  %% files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^EB[1-4]$')==1
  v.data = 1;
  v.cl_id = vs(2);
  v.lev = 2;
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'E';
  if ~isempty(v_info) && isfield(v_info,'probe'), v.sen = num2str(v_info.probe);
  else, v.sen = '1234';
  end
  v.cs = {'na','na','na','na'};
  v.rep = {'scalar','scalar','scalar','scalar'};
  v.units =  {'mV/m','mV/m','na','na'};
  v.si_conv = {'1.0e-3>V m^-1','1.0e-3>V m^-1','',''};
  v.size = [ 1 1 1 1 ];
  v.name = {'Electric field','Electric field','P_Bitmask','P_Quality'};
  v.quant = 'dibe';
  v.labels = {'Ex','Ey','Bitmask','Quality'};
  v.field_name = {'Electric field','Electric field','Bitmask','Quality'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT','FLOAT','INT','INT'};
  v.sigdig = [ 6 6 5 1 ];
  v.ent = {'Electric_Field'};
  v.prop = {'Vector'};
  v.fluc = {'Waveform'};
  v.com = '';
  v.file = 'cef export only';  % caa_export_new() cef export only
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Despun full resolution E
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?di(b)?E(F|LX)?[1-4]p1234$')==1
  v.data = 1;
  switch vvs(1:strfind(vvs,'E')-1) % characters before 'E'
    case 'di'
      v.frame = 'sc';
      v.file = 'mEDSIf';
      v.quant = 'die';
    case 'idi'
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
  if vvs(strfind(vvs,'E')+1)=='F'
    v.cl_id = vvs(strfind(vvs,'E')+2);
    v.name = {'EF_Vec_xy_ISR2'};
    v.labels = {'EF'};
    v.label_1 = {'"EFx", "EFy"'};
    v.field_name = {'Electric field (high-pass filtered)'};
    v.quant = 'dief';
  elseif vvs(strfind(vvs,'E')+1)=='L'
    v.cl_id = vvs(strfind(vvs,'E')+3);
    v.name = {'ELX_Vec_xy_ISR2'};
    v.labels = {'ELX'};
    v.label_1 = {'"ELXx", "ELXy"'};
    v.field_name = {'Electric field (LX)'};
    v.quant = 'dielx';
  else
    v.cl_id = vvs(strfind(vvs,'E')+1); % next character after 'E'
    v.name = {'E_Vec_xy_ISR2'};
    v.labels = {'E'};
    v.label_1 = {'"Ex", "Ey"'};
    v.field_name = {'Electric field'};
  end
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.inst = 'EFW';
  v.sig = 'E';
  if ~isempty(v_info) && isfield(v_info,'probe'), v.sen = num2str(v_info.probe);
  else, v.sen = '1234';
  end
  v.cs = {'ISR2'};
  v.rep = {'xy'};
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.size = 3;
  v.tensor_order = 1;
  v.rep_1 = {'"x", "y"'};
  v.col_labels = {{'x','y','z'}};
  v.ent = {'Electric_Field'};
  v.prop = {'Vector'};
  v.fluc = {'Waveform'};
  v.com = '';
  v.lev = 1;
elseif regexp(vs,'^di(b)?E(F|LX)?[1-4]p1234_info$')==1
  v.data = 0;
  if vs(4)=='F'
    v.cl_id = vs(5);
    v.com = 'E filtered INFO';
    v.quant = 'dief';
    v.file = 'mEDSIf';
  elseif vs(4)=='L'
    v.cl_id = vs(6);
    v.com = 'E LX filtered INFO';
    v.quant = 'dielx';
    v.file = 'mEDSIf';
  elseif vs(3)=='b'
    v.cl_id = vs(3);
    v.com = 'E internal burst INFO';
    v.quant = 'dieburst';
    v.file = 'mEFWburst';
  else
    v.cl_id = vs(4);
    v.com = 'E full res INFO';
    v.quant = 'die';
    v.file = 'mEDSIf';
  end
  v.inst = 'EFW';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% E spectrum
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^diESPEC[1-4]p1234$')==1
  v.data = 0;
  v.cl_id = vs(8);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'E';
  if ~isempty(v_info) && isfield(v_info,'probe'), v.sen = num2str(v_info.probe);
  else, v.sen = '1234';
  end
  v.size = 2;
  v.labels = {'E'};
  v.col_labels = {{'x','y'}};
  v.units =  {'(mV/m)^2/Hz'};
  v.com = 'E-field spectrum';
  v.file = 'mEDSI';
  v.lev = 1;
  v.quant = 'diespec';
elseif regexp(vs,'^diESPEC?[1-4]p1234_info$')==1
  v.data = 0;
  v.cl_id = vs(8);
  v.com = 'E spectrum (DSI) INFO';
  v.quant = 'diespec';
  v.file = 'mEDSI';
  v.inst = 'EFW';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% E spectrum LX
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^diELXSPEC[1-4]p1234$')==1
  v.data = 0;
  v.cl_id = vs(10);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'E';
  if ~isempty(v_info) && isfield(v_info,'probe'), v.sen = num2str(v_info.probe);
  else, v.sen = '1234';
  end
  v.size = 2;
  v.labels = {'E'};
  v.col_labels = {{'x','y'}};
  v.units =  {'(mV/m)^2/Hz'};
  v.com = 'E-field spectrum (LX)';
  v.file = 'mEDSI';
  v.lev = 1;
  v.quant = 'dielxspec';
elseif regexp(vs,'^diELXSPEC?[1-4]p1234_info$')==1
  v.data = 0;
  v.cl_id = vs(10);
  v.com = 'E spectrum LX (DSI) INFO';
  v.quant = 'dielxspec';
  v.file = 'mEDSI';
  v.inst = 'EFW';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% I bias
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
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mEFWR';
  v.file_old = 'mFDM';
  v.quant = 'ibias';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% EFW time
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^EFWT[1-4]$')==1
  v.data = 1;
  v.cl_id = vs(5);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'DSC';
  v.sen = '';
  v.cs = {'na'};
  v.rep = {'scalar'};
  v.units =  {'s'};
  v.si_conv = {'s>s'};
  v.size = 1;
  v.name = {'EFW clock'};
  v.labels = v.name;
  v.field_name = {'EFW clock'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mEFWR';
  v.file_old = 'mFDM';
  v.quant = 'efwt';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% DSC
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^DSC[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(4);
  v.inst = 'EFW';
  v.com = 'DSC';
  v.file = 'mEFWR';
  v.file_old = 'mFDM';
  v.quant = 'dsc';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% HK
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^HK[1-4]$')==1
  v.data = 1;
  v.cl_id = vs(3);
  v.inst = 'EFW';
  v.frame = 'na';
  v.sig = 'HK';
  v.com = 'Housekeeping';
  v.sen = '';
  v.size = [1 1 1 1 1 1 1 1 1 1 1 1];
  v.file = 'mEFW';
  v.name = {'BIAS1','BIAS2','BIAS3','BIAS4','PUCK1','PUCK2','PUCK3','PUCK4','GUARD1','GUARD2','GUARD3','GUARD4'};
  v.labels = v.name;
  v.field_name = {'BIAS current probe1','BIAS current probe2','BIAS current probe3','BIAS current probe4',...
    'PUCK voltage probe1','PUCK voltage probe2','PUCK voltage probe3','PUCK voltage probe4',...
    'GUARD voltage probe1','GUARD voltage probe2','GUARD voltage probe3','GUARD voltage probe4'};
  v.cs = {'na'};
  v.units =  {'nA','nA','nA','nA','V','V','V','V','V','V','V','V'};
  v.si_conv = {'1.0e-9>A','1.0e-9>A','1.0e-9>A','1.0e-9>A','V>V','V>V','V>V','V>V','V>V','V>V','V>V','V>V'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT','FLOAT','FLOAT','FLOAT','FLOAT','FLOAT','FLOAT','FLOAT','FLOAT','FLOAT','FLOAT','FLOAT'};
  v.sigdig = [6 6 6 6 6 6 6 6 6 6 6 6];
  v.quant = 'dsc';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% TMMode
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^mTMode[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(7);
  v.inst = 'EFW';
  v.com = 'MT Mode';
  v.file = 'mEFWR';
  v.file_old = 'mTMode';
  v.quant = 'tmode';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% FDM
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^FDM[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(4);
  v.inst = 'EFW';
  v.file = 'mEFWR';
  v.file_old = 'mFDM';
  v.quant = 'fdm';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Whisper pulses
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^WHIP[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(5);
  v.inst = 'EFW';
  v.com = 'Whisper pulses';
  v.file = 'mEFW';
  v.file_old = 'mFDM';
  v.quant = 'whip';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spike internal burst
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^SPIKE[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(6);
  v.inst = 'EFW';
  v.com = 'Spike';
  v.file = 'mEFWburstR';
  v.file_old = '';
  v.quant = 'spike';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Sweep + dump
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^SWEEP[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(6);
  v.inst = 'EFW';
  v.com = 'Sweep';
  v.file = 'mEFW';
  v.file_old = 'mFDM';
  v.quant = 'sweep';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Burst dump
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^BDUMP[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(6);
  v.inst = 'EFW';
  v.com = 'Burst dump';
  v.file = 'mEFW';
  v.file_old = 'mFDM';
  v.quant = 'bdump';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Bad bias
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^BADBIAS[1-4]p[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(8);
  v.sen = vs(end);
  v.inst = 'EFW';
  v.com = 'Bad bias current';
  v.file = 'mEFW';
  v.file_old = 'mFDM';
  v.quant = 'badbias';
  v.lev = 1;
elseif regexp(vs,'^BADBIASRESET[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(13);
  v.sen = '';
  v.inst = 'EFW';
  v.com = 'Bad bias current due to EFW reset';
  v.file = 'mEFW';
  v.file_old = 'mFDM';
  v.quant = 'badbias';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Probe low density saturation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^PROBELD[1-4]p[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(8);
  v.sen = vs(end);
  v.inst = 'EFW';
  v.com = 'Low density saturation';
  v.file = 'mEFW';
  v.file_old = 'mFDM';
  v.quant = 'probesa';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Probe saturation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^PROBESA[1-4]p[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(8);
  v.sen = vs(end);
  v.inst = 'EFW';
  v.com = 'Probe saturation';
  v.file = 'mEFW';
  v.file_old = 'mFDM';
  v.quant = 'probesa';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Saturation due to probe shadow for SAA=90 deg
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^SAASA(SE|DI)[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(end);
  v.sen = '';
  v.inst = 'EFW';
  if vs(6)=='S', sTmp = 'Single ended'; else, sTmp = 'Differential'; end
  v.com = ['Saturation due to high SAA (' sTmp ' signals)'];
  v.file = 'mEFW';
  v.quant = 'probesa';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Probe high bias saturation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^HBIASSA[1-4]p([1-4]|12|32|34|42)$')==1
  v.data = 0;
  v.cl_id = vs(8);
  v.sen = vs(end);
  v.inst = 'EFW';
  v.com = 'Lhigh bias saturation';
  v.file = 'mEFW';
  v.quant = 'hbiassa';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Probe bad DAC interval
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^BADDAC[1-4]p(12|32|34)$')==1
  v.data = 0;
  v.cl_id = vs(7);
  v.sen = vs(9:end);
  v.inst = 'EFW';
  v.com = 'Bad settings on bias DAC';
  v.file = 'mEFW';
  v.quant = 'baddac';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Nonstandard operations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^NSOPS[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(6);
  v.inst = 'EFW';
  v.com = 'Nonstandard operations';
  v.file = 'mEFW';
  v.quant = 'nsops';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ADC offsets coarse
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
  %% ADC offsets from spinfits
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^Dadc[1-4]p(12|32|34)$')==1
  v.data = 0;
  v.cl_id = vs(5);
  v.sen = vvs(7:8);
  v.inst = 'EFW';
  v.com = 'raw signal DC offset';
  v.file = 'mEDSI';
  v.quant = 'dies';
  v.cs = {'na'};
  v.rep = {'scalar'};
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.size = 1;
  v.name = {'dER_Mag'};
  v.labels = v.name;
  v.field_name = {'raw signal DC offset'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.ent = {'Electric_Field'};
  v.prop = {'Magnitude'};
  v.fluc = {'Waveform'};
  v.com = '';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ADC offsets from spinfits
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^DadcLX[1-4]p(12|32|34|42)$')==1
  v.data = 0;
  v.cl_id = vs(7);
  v.sen = vvs(9:10);
  v.inst = 'EFW';
  v.com = 'raw signal DC offset';
  v.file = 'mEDSI';
  v.quant = 'dies';
  v.cs = {'na'};
  v.rep = {'scalar'};
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.size = 1;
  v.name = {'dER_Mag'};
  v.labels = v.name;
  v.field_name = {'raw signal DC offset'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.ent = {'Electric_Field'};
  v.prop = {'Magnitude'};
  v.fluc = {'Waveform'};
  v.com = '';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Delta offsets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^D[1-4]p12p34$')==1
  v.data = 0;
  v.cl_id = vs(2);
  v.inst = 'EFW';
  v.com = 'Delta offset';
  v.file = 'mEDSI';
  v.name = {'Delta offset'};
  v.labels = v.name;
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.quant = 'dies';
  v.size = 2;
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% DSI offsets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^Ddsi[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(5);
  v.inst = 'EFW';
  v.com = 'DSI offsets';
  v.file = 'mEDSI';
  v.quant = '';
  v.name = {'dE'};
  v.labels = v.name;
  v.units =  {'mV/m'};
  v.size = 1;
  v.lev = 2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% X-TRA DSI offsets
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^DdsiX[1-4]$')==1
  v.data = 0;
  v.cl_id = vs(6);
  v.inst = 'EFW';
  v.com = 'X-TRA DSI offsets';
  v.size = 1;
  v.labels = {'dE'};
  v.units = {'mV/m'};
  v.file = 'mXTRA';
  v.name = {'dE'};
  v.labels = v.name;
  v.units =  {'mV/m'};
  v.size = 1;
  v.quant = 'wake';
  v.lev = 2;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Amplitude correction factor
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
  %% Despun full/4 sec resolution E with assumption E.B = 0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?diE(s)?[1-4]$')
  v.data = 1;
  if vvs(1)=='i'
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
  v.ptype = {'Data','Data','Data'};
  v.valtype = {'FLOAT','FLOAT','FLOAT'};
  v.sigdig = [6 6 6];
  v.com = com_Ez;
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Full/4 sec resolution E in GSE coordinates
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?E(s)?[1-4]$')
  v.data = 1;
  if vvs(1)=='i'
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
  v.ptype = {'Data','Data','Data'};
  v.valtype = {'FLOAT','FLOAT','FLOAT'};
  v.sigdig = [6 6 6];
  v.com = com_Ez;
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ExB
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?VExB(s)?[1-4]$')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'EFW';
  v.frame = 'sc';
  v.sig = 'V=ExB';
  if vvs(5)=='s' || vvs(7)=='s', v.sen = 's';
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
  v.ptype = {'Data','Data'};
  v.valtype = {'FLOAT','FLOAT'};
  v.sigdig = [6 6];
  v.com = com_Ez;
  v.file = 'mEdB';
  if vvs(5)=='s' || vvs(7)=='s', v.quant = 'vedbs';
  else, v.quant = 'vedb';
  end
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Full resolution spacecraft potential and derived density
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(b)?NVps[1-4]$')==1
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
  if vs(1)=='b'
    v.field_name = {'Plasma density (internal burst)',...
      'Spacecraft potential (internal burst)'};
    v.file = 'mEFWburst';
    v.quant = 'pburst';
  else
    v.field_name = {'Plasma density', 'Spacecraft potential'};
    v.file = 'mP';
    v.quant = 'p';
  end
  v.ptype = {'Data','Data'};
  v.valtype = {'FLOAT','FLOAT'};
  v.sigdig = [6 6];
  v.labels = {'Nscp', '-Sc pot'};
  v.name = {'Plasma density','Spacecraft potential'};
  v.com = 'density NVps is derived from Vps based on empirical fit. It is NOT a true density';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% phase and phase_2
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
  v.size = 1;
  if strcmp(vs(2),'t')
    v.name = {'Atwo'};
    v.labels = {'Phase_2'};
  else
    v.name = {'A'};
    v.labels = {'Phase'};
  end
  v.field_name = {'Spacecraft phase'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mA';
  v.quant = 'a';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spin axis orientation
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
  v.size = 1;
  v.name = {'SAX'};
  v.labels = {'Spin axis'};
  v.field_name = {'Spacecraft spin axis'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mEPH';
  v.quant = 'sax';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spacecraft velocity
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
  v.size = 3;
  v.name = {'V'};
  v.labels = v.name;
  v.label_1 = {'"Vx", "Vy", "Vz"'};
  v.col_labels = {{'x','y','z'},''};
  v.field_name = {'Spacecraft velocity'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mR';
  v.quant = 'v';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Spacecraft position
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
  v.size = 3;
  v.name = {'R'};
  v.labels = v.name;
  v.label_1 = {'"Rx", "Ry", "Rz"'};
  v.col_labels = {{'x','y','z'},''};
  v.field_name = {'Spacecraft position'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mR';
  v.quant = 'r';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CIS N PP
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
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.units =  {'cc'};
  v.si_conv = {'1.0e6>1/m^6'};
  v.size = 1;
  v.name = {'N'};
  v.labels = v.name;
  v.com = '';
  v.file = 'mCISR';
  v.file_old = 'mCIS';
  v.quant = 'ncis';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CIS T PP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^T(perp|par)?C(h|p)[1-4]$')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'CIS';
  v.frame = 'na';
  v.cs = {'na'};
  if strcmp(vs(2:4),'per'), comp = 'perp'; cf = 'Perpendicular'; lcomp='\perp';
  else, comp = 'par'; cf = 'Parallel'; lcomp='{||}';
  end
  if vvs(strfind(vvs,'C')+1)=='h' % characters after 'C'
    v.sig = ['T_' comp];
    v.sen = 'HIA';
    v.field_name = {[cf ' ion temperature']};
  else
    v.sig = ['Tp_' comp];
    v.sen = 'COD';
    v.field_name = {[cf ' proton temperature']};
  end
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.units =  {'MK'};
  v.si_conv = {'1.0e6>K'};
  v.size = 1;
  v.name = {['T_' lcomp]};
  v.labels = v.name;
  v.com = '';
  v.file = 'mCISR';
  v.file_old = 'mCIS';
  v.quant = 'tcis';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CIS V PP
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
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.units =  {'km/s'};
  v.si_conv = {'1.0e3>m/s'};
  v.size = 3;
  v.name = {'V'};
  v.labels = v.name;
  v.label_1 = {'"Vx", "Vy", "Vz"'};
  v.col_labels = {{'x','y','z'},''};
  v.com = '';
  v.file = 'mCISR';
  v.file_old = 'mCIS';
  v.quant = 'vcis';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% CIS VxB PP
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
  if vvs(4)=='h'
    v.sig = 'V';
    v.sen = 'HIA';
    v.field_name = {'Ion VxB'};
  else
    v.sig = 'Vp';
    v.sen = 'COD';
    v.field_name = {'Proton VxB'};
  end
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.size = 3;
  v.name = {'E'};
  v.labels = v.name;
  v.label_1 = {'"Ex", "Ey", "Ez"'};
  v.col_labels = {{'x','y','z'},''};
  v.com = '';
  v.file = 'mCIS';
  v.file_old = 'mCIS';
  v.quant = 'vce';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% EDI E PP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(i)?(di)?EDI[1-4]$')
  v.data = 1;
  if vvs(1)=='i'
    vvs = vvs(2:end);
    v.frame = 'inertial';
    v.lev = 0;
    v.file = 'mEDIR';
  else
    v.frame = 'sc';
    v.lev = 1;
    v.file = 'mEDI';
  end

  v.cl_id = vs(end);
  v.inst = 'EDI';
  v.sig = 'E';
  v.sen = '';
  if strcmp(vvs(1:2),'di')
    v.cs = {'ISR2'};
  else
    v.cs = {'GSE'};
  end
  v.units =  {'mV/m'};
  v.si_conv = {'1.0e-3>V m^-1'};
  v.size = 3;
  v.name = {'E'};
  v.labels = v.name;
  v.label_1 = {'"Ex", "Ey", "Ez"'};
  v.col_labels = {{'x','y','z'},''};
  v.field_name = {'Perpendicular Electric field'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file_old = 'mEDI';
  v.quant = 'edi';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% FGM B PP
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?BPP[1-4]$')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'FGM';
  v.frame = 'sc';
  v.sig = 'B';
  v.sen = '';
  if strcmp(vvs(1:2),'di')
    v.cs = {'ISR2'};
  else
    v.cs = {'GSE'};
  end
  v.units =  {'nT'};
  v.si_conv = {'1.0e-9>T'};
  v.size = 3;
  v.name = {'B'};
  v.labels = v.name;
  v.label_1 = {'"Bx", "By", "Bz"'};
  v.col_labels = {{'x','y','z'},''};
  v.field_name = {'Magnetic field PP'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mBPP';
  v.quant = 'b';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% B internal burst level 2 caa_export only no variable in .mat
  %% files
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^BB[1-4]$')==1
  v.data = 1;
  v.cl_id = vs(2);
  v.lev = 2;
  v.inst = 'FGM';
  v.frame = 'sc';
  v.sig = 'B';
  v.sen = '';
  v.cs = {'ISR2','ISR2','ISR2','na','na'};
  v.rep = {'scalar','scalar','scalar','scalar','scalar'};
  v.units =  {'nT','nT','nT','na','na'};
  v.si_conv = {'1.0e-9>T','1.0e-9>T','1.0e-9>T','',''};
  v.size = [ 1 1 1 1 1 ];
  v.name = {'B','B','B','P_Bitmask','P_Quality'};
  v.quant = 'b';
  v.labels = {'Bx','By','Bz','Bitmask','Quality'};
  v.field_name = {'Magnetic field','Magnetic field','Magnetic field','Bitmask','Quality'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT','FLOAT','FLOAT','INT','INT'};
  v.sigdig = [ 6 6 6 5 1 ];
  %	v.ent = {'Instrument'};
  %	v.prop = {'Probe_Potential'};
  %	v.fluc = {'Waveform'};
  v.com = '';
  v.file = 'cef export only';  % caa_export_new() cef export only
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% B
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^diBSC4kHz?')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'FGM';
  v.frame = 'sc';
  v.sig = 'B';
  v.sen = '';
  if strcmp(vvs(1:2),'di')
    v.cs = {'ISR2'};
  else
    v.cs = {'GSE'};
  end
  v.units =  {'nT'};
  v.si_conv = {'1.0e-9>T'};
  v.size = 3;
  v.name = {'B'};
  v.labels = v.name;
  v.label_1 = {'"Bx", "By", "Bz"'};
  v.col_labels = {{'x','y','z'},''};
  v.field_name = {'Magnetic field'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mEFWburst';
  v.quant = 'b';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% B
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^wBSC4kHz?')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'FGM';
  v.frame = 'sc';
  v.sig = 'B';
  v.sen = '';
  if strcmp(vvs(1:2),'di')
    v.cs = {'ISR2'};
  else
    v.cs = {'GSE'};
  end
  v.units =  {'nT'};
  v.si_conv = {'1.0e-9>T'};
  v.size = 3;
  v.name = {'B'};
  v.labels = v.name;
  v.label_1 = {'"Bx", "By", "Bz"'};
  v.col_labels = {{'x','y','z'},''};
  v.field_name = {'Magnetic field'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mEFWburstR';
  v.quant = 'b';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% FGM B full resolution and resampled
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?B(r|rs)?[1-4]$')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'FGM';
  v.frame = 'sc';
  v.sig = 'B';
  v.sen = '';
  if strcmp(vvs(1:2),'di')
    v.cs = {'ISR2'};
  else
    v.cs = {'GSE'};
  end
  v.units =  {'nT'};
  v.si_conv = {'1.0e-9>T'};
  v.size = 3;
  v.name = {'B'};
  v.labels = v.name;
  v.label_1 = {'"Bx", "By", "Bz"'};
  v.col_labels = {{'x','y','z'}};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  if regexp(vs,'^(di)?B(r|rs)[1-4]$')
    v.file = 'mBr';
    if regexp(vs,'^(di)?Brs[1-4]$')
      v.field_name = {'Magnetic field resampled to E (4 sec resolution)'};
      v.com = 'Resampled to E (4 sec resolution)';
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
  %% Spin axis orientation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^WHINAT[1-4]$')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'Whisper';
  v.frame = 'sc';
  v.sig = 'natural';
  v.sen = '';
  v.cs = {'sc'};
  v.units =  {'(V/m)^2/Hz'};
  v.si_conv = {''};
  v.size = 494;
  v.name = {'E'};
  v.labels = {'E'};
  v.field_name = {'E'};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.com = '';
  v.file = 'mWHI';
  v.quant = 'whinat';
  v.lev = 0;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% STAFF SC raw data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif any(regexp(vs,'^wBSC[1-4]$'))
  v.data = 1;
  v.cl_id = vs(5);
  v.inst = 'STAFF';
  v.frame = 'na';
  v.sig = 'B_SC';
  v.sen = 'Bx_By_Bz';
  v.cs = {'SC'};
  v.rep = {'vector'};
  v.units =  {'nT'};
  v.si_conv = {'1.0e-9>T'};
  v.size = 3;
  v.name = {'B_SC'};
  v.labels = {'B SC'};
  v.label_1 = {'"Bx", "By", "Bz"'};
  v.col_labels = {{'x','y','z'}};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.ent = {'Magnetic_Field'};
  v.prop = {'Vector'};
  v.fluc = {'Waveform'};
  v.file = 'mBSCR';
  v.field_name = {'Magnetic field'};
  v.com = '';
  v.lev = 0;
  v.quant = 'bsc';
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% STAFF SC despun data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif regexp(vs,'^(di)?BSC?[1-4]$')
  v.data = 1;
  v.cl_id = vs(end);
  v.inst = 'STAFF';
  v.frame = 'sc';
  v.sig = 'B_SC';
  v.sen = 'Bx_By_Bz';
  if strcmp(vvs(1:2),'di')
    v.cs = {'ISR2'};
  else
    v.cs = {'GSE'};
  end
  v.units =  {'nT'};
  v.si_conv = {'1.0e-9>T'};
  v.size = 3;
  v.name = {'B_SC'};
  v.labels = {'B SC'};
  v.label_1 = {'"Bx", "By", "Bz"'};
  v.col_labels = {{'x','y','z'}};
  v.ptype = {'Data'};
  v.valtype = {'FLOAT'};
  v.sigdig = 6;
  v.file = 'mBSC';
  v.field_name = {'Magnetic field'};
  v.com = '';
  v.quant = 'dibsc';
  v.lev = 1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% ADDITIONAL HELP IN PLOTTING, NOT SPECIFIC TO CLUSTER
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif strcmp(vs,'B') || strcmp(vs,'j') || strcmp(vs,'jz') || strcmp(vs,'jxB')
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
  %irf_log('fcal',['Variable name not recognized: ' vs]);
  varargout={[]};
  return
end

%% Construct the output
if nargout==1, varargout = {v};
else
  % Just print out everthing
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
