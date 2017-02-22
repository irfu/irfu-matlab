function varargout = make_pdist(file,varargin)
% Construct PDist skymap from file name + file path
%   PDist = mms.make_pdist(input)
%   [PDist PDistError] = mms.make_pdist(input)
%   PDistError = mms.make_pdist(input,'error')
%     input - can be a complete file path or a dataobj
%
%   db_info = datastore('mms_db');  
%   file  = [db_info.local_file_db_root '/mms?/fpi/brst/l2/des-dist/2015/10/16/mms?_fpi_brst_l2_des-dist_20151016103254_v2.1.0.cdf'];
%   c_eval('ePDist? = mms.make_pdist(file)',ic)
%
%   db_info = datastore('mms_db');  
%   file  = [db_info.local_file_db_root '/mms?/fpi/brst/l2/des-dist/2015/10/16/mms?_fpi_brst_l2_des-dist_20151016103254_v2.1.0.cdf'];
%   tmpDataObj = dataobj(file);
%   c_eval('ePDist? = mms.make_pdist(tmpDataObj)',ic)
 
loadError = 0;
loadDist = 1;
if ~isempty(varargin) && strcmp(varargin{1},'error')
  loadError = 1;
  loadDist = 0;
end
if nargout == 2
  loadError = 1;
  loadDist = 1;
end


if isa(file,'dataobj')
  tmpDataObj =  file;  
elseif isa(file,'char')
  tmpDataObj = dataobj(file);
else
  error('Input not recognized.')
end

if isempty(tmpDataObj); return; end

fileInfo = regexp(tmpDataObj.GlobalAttributes.Logical_source{1}, ...
  'mms(?<mmsId>[1-4])_fpi_(?<tmMode>(fast|brst))_(?<datalevel>\w*)_(?<detector>\w*)', 'names');

switch(lower(fileInfo.detector))
  case {'des'}
    species = 'electrons';
  case {'dis'}
    species = 'ions';
  otherwise
    irf.log('critical','File not recognized as particle distribution.');
    return;
end

%dataSetName = tmpDataObj.
%'mms1_des_dist_brst'
if loadError
  tmpError = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_disterr_' fileInfo.tmMode]);
  Error = tmpError.data;
end
if loadDist
  tmpDist = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_dist_' fileInfo.tmMode]);
  Dist = tmpDist.data;
else
  Dist = Error;
end


% Collect User data
ud = [];
ud.GlobalAttributes = tmpDist.GlobalAttributes;
ud.CATDESC          = tmpDist.CATDESC;
if isfield(tmpDist,'DISPLAY_TYPE'), ud.DISPLAY_TYPE     = tmpDist.DISPLAY_TYPE; end
ud.FIELDNAM         = tmpDist.FIELDNAM;
ud.VALIDMIN         = tmpDist.VALIDMIN;
ud.VALIDMAX         = tmpDist.VALIDMAX;
if isfield(tmpDist,'LABLAXIS'), ud.LABLAXIS = tmpDist.LABLAXIS; end
if isfield(tmpDist,'LABL_PTR_1'), ud.LABL_PTR_1 = tmpDist.LABL_PTR_1;
elseif isfield(tmpDist,'LABL_PTR_2'), ud.LABL_PTR_2 = tmpDist.LABL_PTR_2;
elseif isfield(tmpDist,'LABL_PTR_3'), ud.LABL_PTR_3 = tmpDist.LABL_PTR_3;
end


% Shift times to center of deltat- and deltat+ for l2 particle
% distributions and moments
if ~isempty(regexp(tmpDist.name,'^mms[1-4]_d[ei]s_','once'))
	if isfield(tmpDist.DEPEND_0,'DELTA_MINUS_VAR') && isfield(tmpDist.DEPEND_0,'DELTA_PLUS_VAR')
        if isfield(tmpDist.DEPEND_0.DELTA_MINUS_VAR,'data') && isfield(tmpDist.DEPEND_0.DELTA_PLUS_VAR,'data')
            irf.log('warning','Times shifted to center of dt-+. dt-+ are recalculated');
            flag_MINUS = 1e3;       flag_PLUS = 1e3;            % s --> ms
            if isfield(tmpDist.DEPEND_0.DELTA_MINUS_VAR, 'UNITS') && isfield(tmpDist.DEPEND_0.DELTA_PLUS_VAR, 'UNITS')
                if strcmp(tmpDist.DEPEND_0.DELTA_MINUS_VAR.UNITS, 's')
                    flag_MINUS = 1e3;           % s --> ms
                elseif strcmp(tmpDist.DEPEND_0.DELTA_MINUS_VAR.UNITS, 'ms')
                    flag_MINUS = 1;
                else
                    irf.log('warning','Epoch_minus_var units are not clear, assume s');
                    flag_MINUS = 1e3;       
                end
                if strcmp(tmpDist.DEPEND_0.DELTA_PLUS_VAR.UNITS, 's')
                    flag_PLUS = 1e3;           % s --> ms
                elseif strcmp(tmpDist.DEPEND_0.DELTA_PLUS_VAR.UNITS, 'ms')
                    flag_PLUS = 1;
                else
                    irf.log('warning','Epoch_plus_var units are not clear, assume s');
                    flag_PLUS = 1e3;
                end
            else
                irf.log('warning','Epoch_plus_var/Epoch_minus_var units are not clear, assume s');
            end                
            toffset = (int64(tmpDist.DEPEND_0.DELTA_PLUS_VAR.data*flag_PLUS)-int64(tmpDist.DEPEND_0.DELTA_MINUS_VAR.data*flag_MINUS))*1e6/2;
            tdiff = (int64(tmpDist.DEPEND_0.DELTA_PLUS_VAR.data*flag_PLUS)+int64(tmpDist.DEPEND_0.DELTA_MINUS_VAR.data*flag_MINUS))*1e6/2;    
            %toffset = (int64(tmpDist.DEPEND_0.DELTA_PLUS_VAR.data)-int64(tmpDist.DEPEND_0.DELTA_MINUS_VAR.data))*1e6/2;
            %tdiff = (int64(tmpDist.DEPEND_0.DELTA_PLUS_VAR.data)+int64(tmpDist.DEPEND_0.DELTA_MINUS_VAR.data))*1e6;
            tmpDist.DEPEND_0.DELTA_MINUS_VAR.data = tdiff;
            tmpDist.DEPEND_0.DELTA_PLUS_VAR.data = tdiff;
            tmpDist.DEPEND_0.data = tmpDist.DEPEND_0.data+toffset;
        end
	end
end

time = tmpDist.DEPEND_0.data;
if ~isa(time,'GenericTimeArray'), time = EpochTT(time); end
dt_minus = tmpDist.DEPEND_0.DELTA_MINUS_VAR.data;
dt_plus = tmpDist.DEPEND_0.DELTA_MINUS_VAR.data;

if is_version_geq(tmpDist.GlobalAttributes.Data_version{:}, '3.1.0')
  % FPI version 3.1.z files or newer
  energy = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_energy_' fileInfo.tmMode]);
  energy = energy.data;
  phi = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_phi_' fileInfo.tmMode]);
  phi = phi.data;
  theta = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_theta_' fileInfo.tmMode]);
  theta = theta.data;
  stepTable = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_steptable_parity_' fileInfo.tmMode]);
  stepTable = stepTable.data;

  % Make energytable from energy0, energy1 and energysteptable
  %energy = repmat(torow(energy),numel(stepTable),1);
  %energy(stepTable==1,:) = repmat(energy,sum(stepTable),1);

  % Construct PDist
  PD = PDist(time,Dist,'skymap',energy,phi,theta);
  PD.userData = ud;
  PD.name = tmpDist.name;
  %PD.units = tmpDist.units;
  PD.units = 's^3/cm^6';
  PD.siConversion = '1e12';
  PD.species = species;
  PD.ancillary.dt_minus = double(dt_minus)*1e-9;
  PD.ancillary.dt_plus = double(dt_plus)*1e-9;
  PD.ancillary.energy = energy;
  PD.ancillary.esteptable = stepTable;
else
  energy0 = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_energy0_' fileInfo.tmMode]);
  energy1 = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_energy1_' fileInfo.tmMode]);
  energy0 = energy0.data; energy1 = energy1.data;
  phi = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_phi_' fileInfo.tmMode]);
  theta = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_theta_' fileInfo.tmMode]);
  stepTable = get_variable(tmpDataObj,['mms' fileInfo.mmsId '_' fileInfo.detector '_steptable_parity_' fileInfo.tmMode]);
  phi = phi.data; theta = theta.data; stepTable = stepTable.data;

  % Make energytable from energy0, energy1 and energysteptable
  energy = repmat(torow(energy0),numel(stepTable),1);
  energy(stepTable==1,:) = repmat(energy1,sum(stepTable),1);

  % Construct PDist
  PD = PDist(time,Dist,'skymap',energy,phi,theta);
  PD.userData = ud;
  PD.name = tmpDist.name;
  %PD.units = tmpDist.units;
  PD.units = 's^3/cm^6';
  PD.siConversion = '1e12';
  PD.species = species;
  PD.ancillary.dt_minus = double(dt_minus)*1e-9;
  PD.ancillary.dt_plus = double(dt_plus)*1e-9;
  PD.ancillary.energy0 = energy0;
  PD.ancillary.energy1 = energy1;
  PD.ancillary.esteptable = stepTable;
end

if loadDist && ~loadError
  varargout = {PD};
elseif loadError && ~loadDist
  varargout = {PD};
else
  PDError = PD.clone(PD.time,Error);  PDError.name = tmpError.name;
  varargout = {PD,PDError};
end

end