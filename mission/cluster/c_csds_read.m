function data=c_csds_read(data_path,start_time,dt,cl_id,quantity)
%C_CSDS_READ read CSDS data from ISDAT or disk
%
% data = c_csds_read(data_path,start_time,dt,cl_id,quantity)
%
% Input:
%	data_path - ISDAT database strings and directories containing CSDS
%	subdirectory with data separated by '|'. These data sources are
%	being tried in the same order. A typical example would be:
%		'db:10|/data/cluster'
%	start_epoch - start time (isdat epoch)
%	dt - time interval in sec
%	cl_id - SC#
%	quantity - one of the following:
%	CSDS PP:
%		'b' : B FGM
%		'edi' : E EDI
%		'ncis_p' : N CIS CODIF(p)
%		'ncis_h' : N CIS HIA
%		'tcis_ppar','tcis_pper' : T CIS CODIF(p) parallel, perpendicular
%		'tcis_hpar','tcis_hper' : T CIS HIA parallel, perpendicular
%		'vcis_p' : V CIS CODIF(p)
%		'vcis_h' : V CIS HIA
%	CSDS SP:
%		'slat' : spin axis latitude
%		'slong' : spin axis longitude
%
% See also C_GET, TOEPOCH
%
% $Id$

% ----------------------------------------------------------------------------
% SPDX-License-Identifier: Beerware
% "THE BEER-WARE LICENSE" (Revision 42):
% <yuri@irfu.se> wrote this file.  As long as you retain this notice you
% can do whatever you want with this stuff. If we meet some day, and you think
% this stuff is worth it, you can buy me a beer in return.   Yuri Khotyaintsev
% ----------------------------------------------------------------------------

narginchk(5,5)

cl_id_s = num2str(cl_id);
pp_infix='PP';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define request parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch (quantity)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSDS PP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'b'
    if (start_time > 1136073599.9) % After 2006-01-01, use the UP files
      r.inst  = 'FGM';
      r.var	= ['B_xyz_gse__C' cl_id_s '_UP_' r.inst];
      r.pr	= 'CSDS_PP';
      pp_infix='UP';
    else
      if ((start_time+dt) > 1136073600)
        warning('c_csds_read:PP_UP_transition',...
          '***** Interval spans both PP and UP files. Data truncated at 2006-01-01.')
      end
      r.inst  = 'FGM';
      r.var	= ['B_xyz_gse__C' cl_id_s '_PP_' r.inst];
      r.pr	= 'CSDS_PP';
    end
  case 'edi'
    r.inst  = 'EDI';
    r.var	= ['E_xyz_gse__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'ncis_p'
    r.inst  = 'CIS';
    r.var	= ['N_p__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'ncis_h'
    r.inst  = 'CIS';
    r.var	= ['N_HIA__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'tcis_ppar'
    r.inst  = 'CIS';
    r.var	= ['T_p_par__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'tcis_pper'
    r.inst  = 'CIS';
    r.var	= ['T_p_perp__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'tcis_hpar'
    r.inst  = 'CIS';
    r.var	= ['T_HIA_par__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'tcis_hper'
    r.inst  = 'CIS';
    r.var	= ['T_HIA_perp__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'vcis_p'
    r.inst  = 'CIS';
    r.var	= ['V_p_xyz_gse__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
  case 'vcis_h'
    r.inst  = 'CIS';
    r.var	= ['V_HIA_xyz_gse__C' cl_id_s '_PP_' r.inst];
    r.pr	= 'CSDS_PP';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CSDS SP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'slat'
    r.inst  = 'AUX';
    r.var	= ['sc_at' cl_id_s '_lat__CL_SP_' r.inst];
    r.pr	= 'CSDS_SP';
  case 'slong'
    r.inst  = 'AUX';
    r.var	= ['sc_at' cl_id_s '_long__CL_SP_' r.inst];
    r.pr	= 'CSDS_SP';
  otherwise
    error('caa:noSuchQuantity','Quantity ''%s'' is not recongized',quantity)
end

% variables
if strcmp(r.pr,'CSDS_SP')
  r.mem	= 'CL';
  r.file	= ['SP/AUX/CL_SP_' r.inst '_'];
  r.sen	= r.var;
elseif strcmp(r.pr,'CSDS_PP')
  r.mem	= ['C' cl_id_s];
  r.file	= ['PP/' r.inst '/C' cl_id_s '/C' cl_id_s '_' pp_infix '_' r.inst '_'];
  r.sen	= r.var;
else
  error('caa:noSuchQuantity','Project ''%s'' is not recongized',r.pr)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = [];
start_date_s = strrep(datestr(fromepoch(start_time),29),'-','');

p = tokenize(data_path,'|');

for i=1:length(p)
  
  if ~isempty(regexp(p{i},':\d{1,2}\>','once'))
    % use ISDAT
    useISDAT = 1;
  else
    if exist(p{i},'dir')
      % read from directory
      useISDAT = 0;
    else
      irf_log('dsrc',['Directory ' p{i} ' does not exist'])
      continue
    end
  end
  
  
  if useISDAT
    %irf_log('dsrc',,'Using ISDAT')
    lasterr('')
    try
      dbase = Mat_DbOpen(p{i});
    catch
      irf_log('dsrc','Cannot open ISDAT database')
      continue
    end
    if exist('dbase','var')
      [t, dat] = isGetDataLite( dbase, start_time, dt, ...
        r.pr, r.mem, r.inst, r.sen);
      Mat_DbClose(dbase)
      clear dbase
      
      if ~isempty(dat)
        % If dat has more the one column, we need to transpose it
        sz = size(dat);
        if sz( sz~=length(t) ) > 1, dat = dat'; end
        
        data = [double(t) double(dat)];
        return
      else
        % irf_log('dsrc','No data')
      end
    end
    
  else
    irf_log('dsrc','Using FILE')
    disp([p{i} '/CSDS/' r.file start_date_s '*']);
    data = irf_cdf_read([p{i} '/CSDS/' r.file start_date_s '*'], r.var,'latest');
    if ~isempty(data)
      data = irf_tlim(data,start_time + [0 dt]);
      return
    else
      % irf_log('dsrc','No data')
    end
  end
end

% if we are here means there was no data
irf_log('dsrc','No data')
