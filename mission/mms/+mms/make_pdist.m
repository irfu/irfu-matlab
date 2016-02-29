function PD = make_pdist(file)
% Construct PDist skymap from file name + file path
%   ePDist = mms.make_pdist(FilePath)
%
%   db_info = datastore('mms_db');  
%   c_eval('ePDist? = mms.make_pdist([db_info.local_file_db_root ''/mms?/fpi/brst/l2/des-dist/2015/10/16/mms?_fpi_brst_l2_des-dist_20151016103254_v2.1.0.cdf''])',ic)

filePathAndName = file;
if strfind(filePathAndName,'des')
  vS = 'des';
  species = 'electrons';
elseif strfind(filePathAndName,'dis')
  vS = 'dis';
  species = 'ions';
end
fileParts = strsplit(filePathAndName,'/');
mmsId = fileParts{end}(4);

%disp('Loading electron distribution...')
tmpDataObj = dataobj(filePathAndName);
Dist = mms.variable2ts(get_variable(tmpDataObj,['mms' mmsId '_' vS '_dist_brst']));
energy0 = get_variable(tmpDataObj,['mms' mmsId '_' vS '_energy0_brst']);
energy1 = get_variable(tmpDataObj,['mms' mmsId '_' vS '_energy1_brst']);
phi = mms.variable2ts(get_variable(tmpDataObj,['mms' mmsId '_' vS '_phi_brst']));
theta = get_variable(tmpDataObj,['mms' mmsId '_' vS '_theta_brst']);
stepTable = mms.variable2ts(get_variable(tmpDataObj,['mms' mmsId '_' vS '_steptable_parity_brst']));
Epoch_plus_var = get_variable(tmpDataObj,'Epoch_plus_var');
Epoch_minus_var = get_variable(tmpDataObj,'Epoch_minus_var');

% Make energytable from energy0, energy1 and energysteptable
energy = repmat(torow(energy0.data),numel(stepTable.data),1);
energy(stepTable.data==1,:) = repmat(energy1.data,sum(stepTable.data),1);
% Shift time so that time stap is in middle of sweep ins tead of in the beginning
dt_shift = 0.5*(double(Epoch_plus_var.data)-double(Epoch_minus_var.data))*1e-3;
dt_minus = double(Epoch_minus_var.data)*1e-3-dt_shift;
dt_plus = double(Epoch_plus_var.data)*1e-3-dt_shift;
% Construct PDist
PD = PDist(Dist.time+dt_shift,Dist.data,'skymap',energy,phi.data,theta.data);
PD.userData = Dist.userData; 
PD.name = Dist.name; 
PD.units = Dist.units;
PD.units = 's^3/cm^6'; 
PD.siConversion = '1e12';
PD.species = species;
PD.ancillary.dt_minus = dt_minus; 
PD.ancillary.dt_plus = dt_plus; 
PD.ancillary.energy0 = energy0.data; 
PD.ancillary.energy1 = energy1.data;
PD.ancillary.energyStepTable = stepTable.data;